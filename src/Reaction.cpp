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
  Reaction::Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
    m_Env(Env),
    m_Name(id),
    restartCalc(true),
    m_PreExp(0.0),
    m_NInf(0.0),
    m_ERConc(1.0),
    m_kfwd(0.0),
    m_HeatOfReaction(0.0),
    m_pMoleculeManager(pMoleculeManager),
    m_pMicroRateCalculator(NULL),
    m_pTunnelingCalculator(NULL),
    m_ppPersist(),
    m_srct(NULL),
    m_rct1(NULL),
    m_rct2(NULL),
    m_pdt1(NULL),
    m_pdt2(NULL),
    m_TransitionState(NULL),
    m_CellKfmc(),
    m_CellKbmc(),
    m_GrainKfmc(),
    m_GrainKbmc()
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
// Read the Molecular data from input stream.
//
bool Reaction::InitializeReaction(PersistPtr ppReac)
{
    m_ppPersist = ppReac;

    Molecule* pMol1(NULL) ;
    Molecule* pMol2(NULL) ;

    //Read reactants
    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
    pMol1 = GetMolRef(ppReactant1);
    if(!pMol1) return false;
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    if(ppReactant2)
    {
      pMol2 = GetMolRef(ppReactant2);
      if(!pMol2){
        {string errorMsg = "Reactant 2 defined in Reaction " + m_Name + " is incomplete.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);}
        return false;
      }
    }

    //Put the CollidingMolecule into m_rct1, even if it is second in datafile
    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);
    }
    else{
      pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
      if(!pColMol){
        meErrorLog.ThrowError(__FUNCTION__, string("At least one of the reactants has to be a colliding molecule"), obError);
        return false;
      }
      m_rct2 = dynamic_cast<ModelledMolecule*>(pMol1);
    }
    m_rct1 = pColMol;

    if (m_rct1 && m_rct2){ // the reactant side has two molecules
      {stringstream errorMsg;
      errorMsg << "Reaction " << m_Name << " has two reactants. ";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

      // check whether there is any SuperMolecule in m_molmap contains pMol1 & pMol2
      string id; //shoud not set any name for it.
      SuperMolecule* pSupMol = NULL;
      while(m_pMoleculeManager->GetNextMolecule(id, pSupMol)){ // get next SuperMolecule
        // if found a SuperMolecule
        ModelledMolecule*  rm1 = pSupMol->getMember1();
        ModelledMolecule*  rm2 = pSupMol->getMember2();
        if (!rm1 && !rm2){// there is no data inside, occupy it!
          pSupMol->setMembers(m_rct1, m_rct2);
          m_srct = pSupMol;
          {stringstream errorMsg;
          errorMsg << "Set members of the SuperMolecule: " << m_srct->getName();
          meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}
          break;
        }
      }
      if (!pSupMol){
        meErrorLog.ThrowError(__FUNCTION__, string("No SuperMolecule was found."), obInfo);
        // there will always at least one SuperMolecule in m_molmap, check the end of addmol()
        // in MoleculeManager.cpp.
        /* need to create one (mark _2007_12_10__17_10_18_)
        write a SuperMolecule creator that acquire a position in the XML
        */
      }
    }
    else{
      string errorMsg = "Reaction " + m_Name + " has only one reactant";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg, obInfo);
    }

    //Read products ... if any.
    pMol1=pMol2=NULL;
    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);

      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      pMol2 = GetMolRef(ppProduct2);

      //Put the Colliding Molecule into m_pdt1, even if it is second in datafile
      pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
      if(pColMol)
      {
        m_pdt1 = pColMol ;
        m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2);
      }
      else
      {
        pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
        if(!pColMol) { // both products are not CollidingMolecule -> dissociation reaction
          m_pdt1 = dynamic_cast<CollidingMolecule*>(pMol1);
          m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2); // do we need to check the existence of m_pdt2 ?
        } else {
          m_pdt1 = pColMol ;
          m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol1);
        }
      }
    }

    // Read the transition state (if present)
    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if (ppTransitionState)
    {
      TransitionState* pTrans = dynamic_cast<TransitionState*>(GetMolRef(ppTransitionState));
      if(pTrans) m_TransitionState = pTrans;
    }

    // Read heat of reaction and rate parameters.

    return ReadRateCoeffParameters(ppReac) ;
}


//
// Locate molecule in molecular map.
//
Molecule* Reaction::GetMolRef(PersistPtr pp)
{
  Molecule* pMol = NULL;

  if(!pp) return NULL;
  PersistPtr ppmol = pp->XmlMoveTo("molecule");
  if(!ppmol) return NULL;

  string pRef = ppmol->XmlReadValue("ref");
  if(pRef.size()){ // if got the name of the molecule
    string pType = ppmol->XmlReadValue("me:type");
    if(pType.size()){ // initialize molecule here with the specified type (need to know m_ppIOPtr)
      PersistPtr ppMolList = m_pMoleculeManager->get_PersistPtr();
      if(!ppMolList)
      {
        meErrorLog.ThrowError(__FUNCTION__, string("No molecules have been specified"), obWarning);
        return false;
      }
      pMol = m_pMoleculeManager->addmol(pRef, pType, ppMolList, getEnv());
    }
  }

  if(!pMol) {
    meErrorLog.ThrowError(__FUNCTION__, string("Cannot find molecule: "), obInfo);
    return NULL;
  }

  return pMol;
}

// Read parameters requires to determine reaction heats and rates.
bool Reaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    // Read the heat of reaction (if present).
    const char* pHeatRxntxt = ppReac->XmlReadValue("me:HeatOfReaction",false);
    if (pHeatRxntxt){
      stringstream s1(pHeatRxntxt);
      s1 >> m_HeatOfReaction ;
    } else { // Calculate heat of reaction.
      double zpe_pdt1 = m_pdt1 ? m_pdt1->get_zpe() : 0.;
      double zpe_pdt2 = m_pdt2 ? m_pdt2->get_zpe() : 0.;
      double zpe_rct1 = m_rct1 ? m_rct1->get_zpe() : 0.;
      double zpe_rct2 = m_rct2 ? m_rct2->get_zpe() : 0.;
      m_HeatOfReaction = zpe_pdt1 + zpe_pdt2 - zpe_rct1 - zpe_rct2;
    }

    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential",false);
    if (pPreExptxt)
    {
      stringstream s2(pPreExptxt); s2 >> m_PreExp ;
    }
    const char* pNInftxt   = ppReac->XmlReadValue("me:nInfinity",false);
    if (pNInftxt)
    {
      stringstream s3(pNInftxt); s3 >> m_NInf ;
    }
    const char* pERConctxt   = ppReac->XmlReadValue("me:excessReactantConc",false);
    if (pERConctxt)
    {
      stringstream s3(pERConctxt); s3 >> m_ERConc ;
    }

    // Determine the method of MC rate coefficient calculation.
    const char* pMCRCMethodtxt = ppReac->XmlReadValue("me:MCRCMethod") ;
    if(pMCRCMethodtxt)
    {
      m_pMicroRateCalculator = MicroRateCalculator::Find(pMCRCMethodtxt);
      if(!m_pMicroRateCalculator)
      {
        cerr << "Unknown method " << pMCRCMethodtxt
             << " for the determination of Microcanonical rate coefficients in reaction "
             << m_Name;
        return false;
      }
    }

    // Determine the method of estimating tunneling effect.
    const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling") ;
    if(pTunnelingtxt)
    {
      m_pTunnelingCalculator = TunnelingCalculator::Find(pTunnelingtxt);
      if(!m_pTunnelingCalculator)
      {
        cerr << "Unknown method " << pTunnelingtxt
             << " for the determination of tunneling coefficients in reaction "
             << m_Name;
        return false;
      }
    }

  return true ;
}


  //
  // Calculate reaction equilibrium constant.
  //
  double Reaction::calcEquilibriumConstant() {
    // equilibrium constant:

    double Keq(0.0) ;

    // Get Canonical partition functions.

    double Qrcts = 1.0;
    if (m_rct2)
      Qrcts = m_srct->grnCanPrtnFn();
    else
      Qrcts = m_rct1->grnCanPrtnFn() ;

    double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

    double mass_rct1 = m_rct1->getMass();
    double mass_rct2 = (m_rct2)? m_rct2->getMass() : 0.0;
    double mass_srct = mass_rct1 + mass_rct2;

    // Calculate the equilibrium constant.
    double beta = getEnv().beta ;

    Keq = Qpdt1/Qrcts ;
    if(debugFlag) ctest << "Keq = " << Keq << endl;

    /* Electronic degeneracies were already accounted for in DOS calculations */
    // Heat of reaction
    Keq *= exp(-beta * m_HeatOfReaction * kJPerMolInRC) ;
    if(debugFlag) ctest << "Keq = " << Keq << endl;

    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tau = 2.0593e19 * pow(2. * M_PI ,1.5);
    if (m_rct2)
      Keq /= (tp_C * pow(mass_rct1 * mass_rct2 / (mass_srct * beta), 1.5));

    if(debugFlag) ctest << "Keq = " << Keq << endl;
    return Keq ;
  }

//
// Access microcanoincal rate coefficients.
//
void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc) {
    calcGrnAvrgMicroRateCoeffs();
    kmc = m_GrainKfmc ;
}

double Reaction::get_ActivationEnergy(void) {
  if (!m_TransitionState) {
    cinfo << "No TransitionState for " << getName() << ", activation energy = 0.";
    return 0.0;
  }
  double zpeReactants = m_rct2 ? m_rct1->get_zpe() + m_rct2->get_zpe() : m_rct1->get_zpe();
  double AE = m_TransitionState->get_zpe() - zpeReactants;
  if(IsNan(AE)){
    cerr << "To use ILT for reaction " << getName() << " the ZPE of the transition state needs to be set.";
    exit(1);
  }
  return AE;
} ;


//
// Calculate grain averaged microcanonical rate coefficients.
//
bool Reaction::calcGrnAvrgMicroRateCoeffs() {
  if (restartCalc){
    if (m_CellKfmc.size()) m_CellKfmc.clear();
    if (m_CellKbmc.size()) m_CellKbmc.clear();
    restartCalc = false; // reset the flag
    // Calculate microcanonical rate coefficients.
    if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this) ||
      (getEnv().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_ppPersist)))
      return false;
    // Calculate Grain-averaged microcanonical rate coefficients.
    return grnAvrgMicroRateCoeffs() ;
  }
}

void Reaction::grnAvg(const int _MG, const int _gsz, const std::vector<double> &DOS, const std::vector<double> &CellRC, std::vector<double> &GrainRC){
  int idx1(0), idx2(0);
  for (int i = 0; i < _MG ; ++i ) {
    int idx3(idx1);

    // Calculate the number of states in a grain.
    double gNOS = .0 ;
    for (int j = 0 ; j < _gsz ; ++j, ++idx1 ) gNOS += DOS[idx1] ;

    // Calculate average energy of the grain if it contains sum states.
    if ( gNOS > 0.0 ) {
      double gSE = .0;
      for (int j= 0 ; j < _gsz ; ++j, ++idx3 ) gSE += CellRC[idx3] * DOS[idx3] ;
      GrainRC[idx2] = gSE/gNOS ;
    }
    idx2++ ;
  }

  // Issue warning if number of grains produced is less that requested.
  if ( idx2 != _MG ) {
    cerr << "Number of grains produced is not equal to that is requested" << endl
             << "Number of grains requested: " << _MG << endl
             << "Number of grains produced : " << idx2 << " in " << getName();
  }
  else{
      //      cinfo << "Number of grains requested: " << MaximumGrain << endl
      //            << "Number of grains produced : " << idx2 << " in " << getName() << endl;
  }
}

//
// Access microcanonical rate coefficients - cell values are averaged
// to give grain values. This code is similar to that in Molecule.cpp
// and this averaging should be done there. SHR 19/Sep/2004.
//
bool Reaction::grnAvrgMicroRateCoeffs() {

    int MaximumGrain = getEnv().MaxGrn;
    int grnSize = static_cast<int>(getEnv().GrainSize);

    // Always calculate forward grain microcanonical rate coefficients.
    m_GrainKfmc.resize(MaximumGrain, 0.);
    if (m_pdt1)
      m_GrainKbmc.resize(MaximumGrain, 0.);

    // Extract density of states of equilibrium molecule.
    std::vector<double> rctCellDOS;
    std::vector<double> pdtCellDOS;
    if (m_rct2)
      m_srct->getCellDensityOfStates(rctCellDOS);
    else
      m_rct1->getCellDensityOfStates(rctCellDOS);

    if (m_pdt1)
      m_pdt1->getCellDensityOfStates(pdtCellDOS);

    // Check if there are enough cells.
    if (grnSize < 1) {
      cerr << "The requested grain size is invalid.";
      exit(1) ;
    }

    grnAvg(MaximumGrain, grnSize, rctCellDOS, m_CellKfmc, m_GrainKfmc);
    if (m_pdt1)
      grnAvg(MaximumGrain, grnSize, pdtCellDOS, m_CellKbmc, m_GrainKbmc);

    if (getEnv().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }

    if (getEnv().kbEGrainsEnabled && m_pdt1){
      ctest << "\nk_b(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKbmc[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
}

//
// Add microcanonical terms to collision operator
//
void Reaction::AddMicroRates(dMatrix *CollOptr,
                               isomerMap &isomermap,
                               const double rMeanOmega)
{
    // Calculate Microcanonical rate coefficients.
    calcGrnAvrgMicroRateCoeffs() ;

    // Add microcanonical rates to the collision operator.
    AddReactionTerms(CollOptr, isomermap, rMeanOmega) ;
}

//
// DetailedBalance
//
void Reaction::detailedBalance(const int dir){
  // if dir == -1 then do backward -> forward conversion

  int MaximumCell = getEnv().MaxCell;

  if (dir == -1) {
    m_CellKfmc.resize(MaximumCell, 0.0);

    for (int i = 0; i < MaximumCell; ++i)
      m_CellKfmc[i] = m_CellKbmc[i] * m_srct->m_cellDOS[i] / m_pdt1->m_cellDOS[i];

    if (getEnv().kfECellsEnabled){
      ctest << "k_f(e) cells:\n{\n";
      for (int i = 0; i < MaximumCell; ++i)
          ctest << m_CellKfmc[i] << endl;
      ctest << "}" << endl;
    }
  }
  else{
    m_CellKbmc.resize(MaximumCell, 0.0);

    for (int i = 0; i < MaximumCell; ++i)
      m_CellKbmc[i] = m_CellKfmc[i] * m_pdt1->m_cellDOS[i] / m_srct->m_cellDOS[i];

    if (getEnv().kbECellsEnabled){
      ctest << "k_b(e) cells:\n{\n";
      for (int i = 0; i < MaximumCell; ++i)
          ctest << m_CellKbmc[i] << endl;
      ctest << "}" << endl;
    }
  }

}


}//namespace
