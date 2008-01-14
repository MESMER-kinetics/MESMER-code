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
m_pMoleculeManager(pMoleculeManager),
m_srct(NULL),
m_rct1(NULL),
m_rct2(NULL),
m_pdt1(NULL),
m_pdt2(NULL),
m_TransitionState(NULL),
m_HeatOfReaction(0.0),
m_kfwd(0.0),
m_CellKfmc(),
m_GrainKfmc(),
m_ActEne(0.0),
m_PreExp(0.0),
m_NInf(0.0),
m_ERConc(0.0),
m_ppPersist(),
m_pMicroRateCalculator(NULL)
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

    //Read in Kinetic rate parameters, if present
    m_ActEne = std::numeric_limits<double>::quiet_NaN();//means not set

    const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy",false);
    if (pActEnetxt)
    {
      stringstream s1(pActEnetxt); s1 >> m_ActEne ;
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
        {stringstream errorMsg;
        errorMsg << "Unknown method " << pMCRCMethodtxt
                 << " for the determination of Microcanonical rate coefficients in reaction "
                 << m_Name;
        meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);}
        return false;
      }
    }

	return true ;
}


  //
  // Calculate reaction equilibrium constant.
  //
  double Reaction::calcEquilibriumConstant() {

    double Keq(0.0) ;

    // Get Canonical partition functions.

    double Qrct1 = m_rct1->grnCanPrtnFn() ;
    double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

    double Qrct2 = (m_rct2)? m_rct2->grnCanPrtnFn() : 1.0 ;
    double Qpdt2 = (m_pdt2)? m_pdt2->grnCanPrtnFn() : 1.0 ;

    // Calculate the equilibrium constant.
    if(0){stringstream errorMsg;
    errorMsg << "Qrct1 = " << Qrct1 << ", Qpdt1 = " << Qpdt1 << ", Qrct2 = " << Qrct2 << ", Qpdt2 = " << Qpdt2;
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

    double beta = getEnv().beta ;

    Keq = (Qpdt1*Qpdt2)/(Qrct1*Qrct2) ;     // Whatever reaction.
    Keq *= exp(-beta*m_HeatOfReaction) ;

    // check in ASSOCIATION reaction, if the partition function of the SuperMolecule is equal to the product of the
    // partitions of the two reactants.
    /*    if (m_rct1 && m_rct2 && m_pdt1 && !m_pdt2){
    double Qrcts = m_srct->grnCanPrtnFn();
    if (Qrcts != Qrct1){
    double diff = Qrcts - Qrct1;
    stringstream errorMsg;
    errorMsg << "Partition function of the SuperMolecule is not consistent with the product of partition functions of the reactants. Diff = " << diff;
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    }
    } */

    return Keq ;
}

//
// Access microcanoincal rate coefficients.
//
void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc) {
    calcGrnAvrgMicroRateCoeffs();
    kmc = m_GrainKfmc ;
}

//
// Calculate grain averaged microcanonical rate coefficients.
//
bool Reaction::calcGrnAvrgMicroRateCoeffs() {

    // Calculate microcanonical rate coefficients.
    if (m_CellKfmc.size()==0)
    {
      if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this, m_CellKfmc, getEnv()) ||
        (getEnv().microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_CellKfmc, m_ppPersist, getEnv())))
        return false;
    }

    // Calculate Grain averages of microcanonical rate coefficients.
    if (m_GrainKfmc.size()==0)
      return grnAvrgMicroRateCoeffs() ;
    return true;
}

//
// Access microcanonical rate coefficients - cell values are averaged
// to give grain values. This code is similar to that in Molecule.cpp
// and this averaging should be done there. SHR 19/Sep/2004.
//
bool Reaction::grnAvrgMicroRateCoeffs() {

    int MaximumGrain = getEnv().MaxGrn;
    double currentGrainSize = getEnv().GrainSize;
    m_GrainKfmc.resize(MaximumGrain, 0.);

    // Extract density of states of equilibrium molecule.

    vector<double> cellDOS(getEnv().MaxCell,0.0) ;
    m_rct1->getCellDensityOfStates(cellDOS) ;

    // Check that there are enough cells.

    if (currentGrainSize < 1) {
        meErrorLog.ThrowError(__FUNCTION__, string("Not enought Cells to produce requested number of Grains."), obError);
      exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i = 0; i < MaximumGrain ; ++i ) {

      int idx3 = idx1 ;

      // Calculate the number of states in a grain.

      double gNOS = .0 ;
      for (int j = 0 ; j < currentGrainSize ; ++j, ++idx1 )
        gNOS += cellDOS[idx1] ;

      // Calculate average energy of the grain if it contains sum states.
        // need to think about how to deal with DOS/ENE of atoms. (there should 
		// have no rovibrational DOS/ENE for atoms)
      if ( gNOS > 0.0 ) {

        double gSE = .0;
        for (int j= 0 ; j < currentGrainSize ; ++j, ++idx3 )
          gSE += m_CellKfmc[idx3] * cellDOS[idx3] ;

        m_GrainKfmc[idx2] = gSE/gNOS ;
      }
      idx2++ ;
    }

    // Issue warning if number of grains produced is less that requested.
    if ( idx2 != MaximumGrain ) {
      {stringstream errorMsg;
      errorMsg << "Number of grains produced is not equal to that requested" << endl
               << "Number of grains requested: " << MaximumGrain << endl
               << "Number of grains produced : " << idx2 << " in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);}
    }
    else{
        //      {stringstream errorMsg;
        //      errorMsg << "Number of grains requested: " << MaximumGrain << endl
        //               << "Number of grains produced : " << idx2 << " in " << getName();
        //      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}
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


}//namespace
