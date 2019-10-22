//-------------------------------------------------------------------------------------------
//
// Reaction.cpp
//
// Author: Struan Robertson
// Date:   23/Feb/2003
//
// This file contains the implementation of the Reaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "Reaction.h"
#include "ParseForPlugin.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{

  Reaction::Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id)
    :m_ppPersist(),
    m_TransitionState(NULL),
		m_ExcessReactant(NULL),
		m_pMoleculeManager(pMoleculeManager),
    m_pMicroRateCalculator(NULL),
    m_pTunnelingCalculator(NULL),
    m_FluxCellZPE(0.0),
    m_FluxGrainZPE(0.0),
    m_FluxCellOffset(0),
    m_CellFlux(),
    m_GrainFlux(),
    m_GrainKfmc(),
    m_MtxGrnKf(),
		m_ERConc(0.0),
		m_Env(Env),
    m_Flags(Flags),
    m_Name(id),
    m_reCalcMicroRateCoeffs(true),
    m_UsesProductProperties(true),
    m_GrnFluxFirstNonZeroIdx(0),
    m_EffGrainedFwdThreshold(0),
    m_EffGrainedRvsThreshold(0)
  {}

  Reaction::~Reaction(){}

  /*
  Reaction::Reaction(const Reaction& reaction) {
  // Copy constructor - define later SHR 23/Feb/2003
  }

  Reaction& Reaction::operator=(const Reaction& reaction) {
  // Assignment operator - define later SHR 23/Feb/2003

  return *this ;
  }
  */


  //
  // Locate molecule in molecular map.
  //
  Molecule* Reaction::GetMolRef(PersistPtr pp, const char* defaultType)
  {
    Molecule* pMol = NULL;

    if(!pp) return NULL;
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) return NULL;

    const char* reftxt = ppmol->XmlReadValue("ref");//using const char* in case NULL returned
    if(reftxt) // if got the name of the molecule
    {
      const char* typetxt = ppmol->XmlReadValue("me:type",optional);
      if(!typetxt)
        typetxt = ppmol->XmlReadValue("role");
      if(!typetxt && defaultType)
        typetxt=defaultType;
      if(typetxt){ // initialize molecule here with the specified type (need to know m_ppIOPtr)
        PersistPtr ppMolList = m_pMoleculeManager->get_PersistPtr();
        if(!ppMolList)
        {
          cerr << "No molecules have been specified." << endl;
          return NULL;
        }
        pMol = m_pMoleculeManager->addmol(string(reftxt), string(typetxt), getEnv(), getFlags());
        if(string(typetxt)==string("excessReactant"))
          m_ExcessReactant = pMol;
      }
    }

    if(!pMol) {
      cinfo << "Failed to get a molecular reference." << endl;
      return NULL;
    }

    return pMol;
  }

  //
  // Calculate grain averaged microcanonical rate coefficients.
  //
  bool Reaction::calcGrnAvrgMicroRateCoeffs() {
    if (m_reCalcMicroRateCoeffs){
      if (m_CellFlux.size()) m_CellFlux.clear();

      // Calculate microcanonical rate coefficients.
      if(!m_pMicroRateCalculator->calculateMicroCnlFlux(this))
        return false;

      // report Transition State Flux in cells to test output
      const int MaximumCell = getEnv().MaxCell;
      if (getFlags().cellFluxEnabled){
        ctest << "\nFlux(e) cells for " << getName() << ":\n{\n";
        for (int i = 0; i < MaximumCell; ++i){
          ctest << m_CellFlux[i] << endl;
        }
        ctest << "}\n";
      }

      // Calculate Grain-averaged microcanonical flux coefficients.
      if (!grnAvrgMicroFluxCoeffs())
        return false;

      // test grained microcanonical rate coefficients
      if (getFlags().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_ppPersist) )
        return false;
    }
    m_reCalcMicroRateCoeffs = false; // reset the flag
    return true;
  }

  //
  // Access microcanonical flux coefficients - cell values are averaged
  // to give grain values. 
  //
  bool Reaction::grnAvrgMicroFluxCoeffs() {

    // convert flux from cells to grains
    fluxCellToGrain();

    // Calculate forward and backward grained microcanonical rate coefficients
    calcGrainRateCoeffs();

    return true;
  }

  // set the bottom energy of m_CellFlux
  void Reaction::setCellFluxBottom(const double fluxBottomZPE){
    m_FluxCellZPE = fluxBottomZPE;
    m_FluxGrainZPE = fluxBottomZPE / getEnv().GrainSize ; //convert to grain
    m_FluxCellOffset = size_t(fmod(fluxBottomZPE, double(getEnv().GrainSize))/getEnv().CellSize) ;
  }

  // Calculate grain flux by summing over cells belong to each grain 
  // taking account of the cell offset against PES grid by altering
  // the range of first grain average.
  void Reaction::fluxCellToGrain() {

    const size_t maxGrn       = getEnv().MaxGrn;
    const size_t cellPerGrain = getEnv().cellPerGrain() ;
    const size_t cellOffset   = getFluxCellOffset();

    m_GrainFlux.clear();
    m_GrainFlux.resize(maxGrn, 0.0);

    for (size_t i(0), cIdx(0) ; i < maxGrn ; ++i) {
	  const size_t cellRange = (i == 0) ? cellPerGrain - cellOffset : cellPerGrain ; 
      for (size_t j(0) ; j < cellRange; ++j, ++cIdx) {
        m_GrainFlux[i] += m_CellFlux[cIdx];
      }
    }

    if (getFlags().grainFluxEnabled){
      ctest << "\nFlux(e) grains for " << getName() << ":\n{\n";
      for (size_t i(0) ; i < maxGrn; ++i){
        ctest << m_GrainFlux[i] << endl;
      }
      ctest << "}\n";
    }
  }
 
  void Reaction::calcFluxFirstNonZeroIdx(void) {
    double thresh = get_ThresholdEnergy();
    double RxnHeat = getHeatOfReaction();
    if(thresh<0.0)
      m_GrnFluxFirstNonZeroIdx = int(-thresh/m_Env.GrainSize);
    else if(thresh>0.0 && thresh<RxnHeat)
      m_GrnFluxFirstNonZeroIdx = int(RxnHeat - thresh)/m_Env.GrainSize;
    else
      m_GrnFluxFirstNonZeroIdx = 0;
  };

  // Read excess reactant concentration
  bool Reaction::ReadExcessReactantConcentration(PersistPtr ppReac){
    const char* pERConctxt = ppReac->XmlReadValue("me:excessReactantConc");
    if (!pERConctxt){
      cerr << "Concentration of excess reactant has not been specified.";
      return false;
    } else {
      stringstream s3(pERConctxt) ;
      s3 >> m_ERConc ;
    }
    return true;
  }

  // Read parameters requires to determine reaction heats and rates.
  bool Reaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    // Determine the method of MC rate coefficient calculation.
    // The name of the method may be in a text element e.g.<me:MCRCMethod>SimpleRRKM</me:MCRCMethod>
    // OR in a name or a xsi:type attribute e.g <me:MCRCMethod xsi:type ="me:MesmerILT">

    // Read the transition state (if present)
    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if(!ppTransitionState)
      ppTransitionState = ppReac->XmlMoveTo("transitionState") ;
    if (ppTransitionState)
    {
      Molecule* pTrans = GetMolRef(ppTransitionState,"transitionState");
      if(pTrans) {
        m_TransitionState = pTrans ;
        m_TransitionState->activateRole(string("transitionState"));
      }
    }

    // Determine the method of MC rate coefficient calculation.
    // The name of the method may be in a text element e.g.<me:MCRCMethod name="RRKM"/>
    // OR in a name or a xsi:type attribute e.g <me:MCRCMethod xsi:type ="me:MesmerILT">

    m_pMicroRateCalculator = ParseForPlugin<MicroRateCalculator>(this, "me:MCRCMethod", ppReac);

    // Determine the method of estimating tunneling coefficients. Note data may be in TS.
    m_pTunnelingCalculator = ParseForPlugin<TunnelingCalculator>(this, "me:tunneling", ppReac, optional);
    if(!m_pTunnelingCalculator)
      cinfo << "No tunneling method used for " << getName() << endl;

    if ((getReactionType() == ASSOCIATION || getReactionType() == IRREVERSIBLE_EXCHANGE 
      || getReactionType() == BIMOLECULAR_SINK || getReactionType() == PSEUDOISOMERIZATION
	    || getReactionType() == SECONDORDERASSOCIATION || getReactionType() == BIMOLECULAR_EXCHANGE)
      && (m_ERConc==0.0 || IsNan(m_ERConc)))
    {
      // If not already read in the MicroRateCalculator
      cinfo << "Not a unimolecular reaction: look for excess reactant concentration." << endl;
      if (!ReadExcessReactantConcentration(ppReac)) return false;
    }

    return true ;
  }


  void Reaction::setUsesProductProperties(bool b)
  {
    m_UsesProductProperties = b;

    //Ensure appropriate product properties have been read in..
    if(b)
      getHeatOfReaction();
  }

  string Reaction::getReactionString(reactionType type)
  {
    string s;
    int n;
    vector<Molecule*> reactants, products;
    if(type != productsOnly)
    {
      n = get_reactants(reactants);
      for(int i=0; i<n; ++i)
      {
        s += reactants[i]->getName();
        if(i<n-1)
          s += " + ";
      }
    }
    if(type==all)
      s += " => ";
    if(type==rev)
      s += " => ";
    if(type!=reactantsOnly)
    {
      n = get_products(products);
      for(int i=0; i<n; ++i)
      {
        s += products[i]->getName();
        if(i<n-1)
          s += " + ";
      }
    }
    return s;
  }
}//namespace
