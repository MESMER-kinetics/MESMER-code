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

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{

  Reaction::Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id)
    :m_ppPersist(),
    m_TransitionState(NULL),
    m_pMoleculeManager(pMoleculeManager),
    m_pMicroRateCalculator(NULL),
    m_pTunnelingCalculator(NULL),
    m_pCrossingCalculator(NULL),
    m_FluxCellZPE(0.0),
    m_FluxGrainZPE(0.0),
    m_FluxCellOffset(0),
    m_CellFlux(),
    m_GrainFlux(),
    m_GrainKfmc(),
    m_ERConc(0.0),
    m_fwdGrnCanonicalRate(0.0),
    m_rvsGrnCanonicalRate(0.0),
    m_fwdCellCanonicalRate(0.0),
    m_rvsCellCanonicalRate(0.0),
    m_Env(Env),
    m_Flags(Flags),
    m_Name(id),
    m_reCalcMicroRateCoeffs(true),
    m_kfwd(0.0),
    m_GrnFluxFirstNonZeroIdx(0),
    m_EffGrainedFwdThreshold(0),
    m_EffGrainedRvsThreshold(0),
    m_UsesProductProperties(true)
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
      if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this))
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

      // Calculate Grain-averaged microcanonical rate coefficients.
      if (!grnAvrgMicroRateCoeffs())
        return false;

      // test grained microcanonical rate coefficients
      if (getFlags().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_ppPersist) )
        return false;
    }
    m_reCalcMicroRateCoeffs = false; // reset the flag
    return true;
  }

  //
  // Access microcanonical rate coefficients - cell values are averaged
  // to give grain values. This code is similar to that in Molecule.cpp
  // and this averaging should be done there. SHR 19/Sep/2004.
  //
  bool Reaction::grnAvrgMicroRateCoeffs() {
    // This grain averaging of the microcanonical rate coefficients is
    // based on the view from the species that is
    // moving in the current reaction toward the opposite species.

    std::vector<double> shiftedCellFlux;
    shiftCellFlux(shiftedCellFlux);

    // convert flux from cells to grains
    fluxCellToGrain(shiftedCellFlux);

    // Calculate forward and backward grained microcanonical rate coefficients
    calcGrainRateCoeffs();

    return true;
  }

  // set the bottom energy of m_CellFlux
  void Reaction::setCellFluxBottom(const double fluxBottomZPE){
    m_FluxCellZPE = fluxBottomZPE;
    m_FluxGrainZPE = fluxBottomZPE / getEnv().GrainSize ; //convert to grain
    m_FluxCellOffset = int(fmod(fluxBottomZPE, getEnv().GrainSize));
  }

  // shift transition state cell flux
  void Reaction::shiftCellFlux(std::vector<double>& shiftedCellFlux){
    int cellOffset = getFluxCellOffset();
    const int MaximumCell  = getEnv().MaxCell;
    for(int i = 0; i < cellOffset; ++i){
      shiftedCellFlux.push_back(0.0);
    }
    for(int i = cellOffset, j = 0; i < MaximumCell; ++i, ++j){
      shiftedCellFlux.push_back(m_CellFlux[j]);
    }
  }

  // calculate flux in grains
  void Reaction::fluxCellToGrain(const std::vector<double>& shiftedCellFlux)
  {
    const int maxGrn = getEnv().MaxGrn;
    const int grnSiz = getEnv().GrainSize;

    // resize m_GrainFlux to maxGrn and initialize all members to zero
    m_GrainFlux.clear();
    m_GrainFlux.resize(maxGrn, 0.0);

    int cIdx = 0; // cell iterator

    for (int i = 0; i < maxGrn ; ++i) {
      for (int j = 0; j < grnSiz; ++j, ++cIdx) {
        m_GrainFlux[i] += shiftedCellFlux[cIdx];
      }
    }

    if (getFlags().grainFluxEnabled){
      ctest << "\nFlux(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < maxGrn; ++i){
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
    const char* pMCRCMethodtxt;
    PersistPtr pp = ppReac->XmlMoveTo("me:MCRCMethod");
    if(pp)
      pMCRCMethodtxt = pp->XmlReadValue("xsi:type",optional);
    if(!pMCRCMethodtxt)
      pMCRCMethodtxt = pp->XmlReadValue("name",optional);
    if(!pMCRCMethodtxt)
      pMCRCMethodtxt = ppReac->XmlReadValue("me:MCRCMethod", true, MicroRateCalculator::typeID()) ;
    if(pMCRCMethodtxt)
    {
      if(pMCRCMethodtxt[2]==':')
        pMCRCMethodtxt+=3; //Remove prefix "me:"
      m_pMicroRateCalculator = MicroRateCalculator::Find(pMCRCMethodtxt);
      if(!m_pMicroRateCalculator)
      //{
      //  cerr << "Unknown method " << pMCRCMethodtxt
      //    << " for the determination of Microcanonical rate coefficients in reaction "
      //    << getName() << endl;
        return false;
      //}
    }

    // Read the transition state (if present)
    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if (ppTransitionState)
    {
      Molecule* pTrans = GetMolRef(ppTransitionState,"transitionState");
      if(pTrans) m_TransitionState = pTrans;
    }

    //---------------------------------------------------------
    // Microcanonical rate constants methods dependent section.
      
      if (!m_pMicroRateCalculator->ReadParameters(this)) return false;
      
      // Determine the method of estimating crossing coefficients.
      const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling", optional) ;
      if(pTunnelingtxt)
      {
        m_pTunnelingCalculator = TunnelingCalculator::Find(pTunnelingtxt);
        if(!m_pTunnelingCalculator)
        {
          cerr << "Unknown method " << pTunnelingtxt
            << " for the determination of tunneling coefficients in reaction "
            << getName();
          return false;
        }
      }
      else{
        cinfo << "No tunneling method requested for " << getName() << endl;
      }

      // Determine the method of estimating crossing coefficients.
      const char* pCrossingtxt = ppReac->XmlReadValue("me:crossing", optional) ;
      if(pCrossingtxt)
      {
        m_pCrossingCalculator = CrossingCalculator::Find(pCrossingtxt);
        if(!m_pCrossingCalculator)
        {
          cerr << "Unknown method " << pCrossingtxt
            << " for the determination of crossing coefficients in reaction "
            << getName();
          return false;
        }
      }
      else{
        cinfo << "No crossing method requested for " << getName() << endl;
      }

    //
    //---------------------------------------------------------



    if (getReactionType() == ASSOCIATION || getReactionType() == IRREVERSIBLE_EXCHANGE || getReactionType() == BIMOLECULAR_SINK  || getReactionType() == PSEUDOISOMERIZATION){
      cinfo << "Not a unimolecular reaction: look for excess reactant concentration." << endl;
      if (!ReadExcessReactantConcentration(ppReac)) return false;
    }

    return true ;
  }

  //Returns true if the XML input contains ZPE for all products, but does not read them
  bool Reaction::ProductEnergiesSupplied() const
  {
    vector<Molecule*> prods;
    unsigned nprods = get_products(prods);
    for(unsigned i=0;i<nprods;++i)
      if(!(prods[i]->get_PersistentPointer())->XmlMoveToProperty("me:ZPE"))
        return false;

    return true;
  }

  void Reaction::setUsesProductProperties(bool b)
  {
    m_UsesProductProperties = b;

    //Ensure appropriate product properties have been read in..
    if(b)
      getHeatOfReaction();
  }

}//namespace
