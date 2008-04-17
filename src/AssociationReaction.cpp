//-------------------------------------------------------------------------------------------
//
// AssociationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the AssociationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "AssociationReaction.h"
#include <math.h>

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular for association reaction data from input stream.
  // Note: the convention adopted here is that there are two reactants
  // and one product (adduct).
  //
  bool AssociationReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.

    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      string errorMsg = "Association reaction " + getName() + " has no reactants.";
      meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
      return false;
    }
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    Molecule* pMol2 = GetMolRef(ppReactant2);
    if(!pMol2)
    {
      string errorMsg = "Association reaction " + getName() + " does not have two reactants.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
      return false;
    }

    // Put the CollidingMolecule into m_rct1, even if it is second in datafile

    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);
    } else {
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
      errorMsg << "Reaction " << getName() << " has two reactants. ";
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
      }
    }
    else{
      string errorMsg = "Reaction " + getName() + " has only one reactant";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg, obInfo);
    }

    //Read product details.

    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    pMol1 = GetMolRef(ppProduct1);
    if (!pMol1) {
      string errorMsg = "Association reaction " + getName() + " has no product.";
      meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
      return false;
    }
    PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
    if(ppProduct2)
    {
      string errorMsg = "Association reaction " + getName() + " has more than one product.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
      return false;
    }

    // Save product as CollidingMolecule.

    pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_pdt1 = pColMol;
    } else {
      meErrorLog.ThrowError(__FUNCTION__, string("Isomer product must be a colliding molecule"), obError);
      return false;
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
  // Add (REVERSIBLE) association reaction terms to collision matrix.
  //
  void AssociationReaction::AddReactionTerms(dMatrix      *CollOptr,
    isomerMap    &isomermap,
    const double rMeanOmega)
  {
    // Locate isomers in system matrix.
    const int pdtLoc =      isomermap[dynamic_cast<CollidingMolecule*>(m_pdt1)] ;
    const int jj     = (*m_sourceMap)[dynamic_cast<SuperMolecule    *>(m_srct)] ;

    // Get equilibrium constant.
    const double Keq = calcEquilibriumConstant() ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn ;
    vector<double> adductBoltz ;
    m_pdt1->grnBoltzDist(adductBoltz) ;

    const int colloptrsize = dynamic_cast<CollidingMolecule*>(m_pdt1)->get_colloptrsize() ; // collision operator size of the isomer

    double DissRateCoeff(0.0) ;

    // Multiply equilibrum constant by concentration of excess reactant.
    // concentration of excess reactant should be in molec/cm3. This gives
    // a dimensionless pseudo-isomerization equilibrium constant.

    for ( int i = 0; i < colloptrsize; ++i) {
      int ll = i;
      int ii(pdtLoc + i) ;

      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKbmc[ll] ;                               // Backward loss of the adduct
      (*CollOptr)[jj][ii]  = rMeanOmega * m_GrainKbmc[ll] * sqrt(adductBoltz[ll] * Keq) ; // Reactive gain of the source
      (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                                        // Reactive gain
      DissRateCoeff       += m_GrainKbmc[ll] * adductBoltz[ll] ;
    }
    (*CollOptr)[jj][jj] -= (rMeanOmega * DissRateCoeff * Keq);       // Backward loss reaction from detailed balance.
    
  }

  //
  // Calculate reaction equilibrium constant for the general reaction
  //        A + B  <===> C
  //
  double AssociationReaction::calcEquilibriumConstant() {

    // equilibrium constant:
    long double Keq(0.0) ;

    // Get Canonical partition functions.
    const long double Qrct1 = m_rct1->grnCanPrtnFn();
    const long double Qrct2 = m_rct2->grnCanPrtnFn();
    const long double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

    const long double mass_rct1 = m_rct1->getMass();
    const long double mass_rct2 = m_rct2->getMass();
    const long double mass_srct = mass_rct1 + mass_rct2;

    // Calculate the equilibrium constant.
    const long double beta = getEnv().beta ;

    Keq = Qpdt1 / (Qrct1 * Qrct2);

    /* Electronic degeneracies were already accounted for in DOS calculations */
    // Heat of reaction
    double HeatOfReaction = getHeatOfReaction() ;
    long double _expon = -beta * HeatOfReaction * kJPerMolInRC;
    Keq *= exp(_expon) ;

    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tau = 2.0593e19 * pow(2. * M_PI ,1.5);
    Keq /= (tp_C * pow(mass_rct1 * mass_rct2 / (mass_srct * beta), 1.5l));

    const double excess = getExcessReactantConc();
    Keq *= excess ;

    return (double) Keq ;

  }

  // Read parameters requires to determine reaction heats and rates.
  bool AssociationReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    // Read the heat of reaction (if present).
    const char* pHeatRxntxt = ppReac->XmlReadValue("me:HeatOfReaction",false);
    if (pHeatRxntxt){
      stringstream s1(pHeatRxntxt);
      double value = 0.0; s1 >> value ; setHeatOfReaction(value);
    } else { // Calculate heat of reaction.
      double zpe_pdt1 = m_pdt1->get_zpe();
      double zpe_rct1 = m_rct1->get_zpe();
      double zpe_rct2 = m_rct2->get_zpe();
      setHeatOfReaction(zpe_pdt1 - zpe_rct1 - zpe_rct2);
    }

    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential",false);
    if (pPreExptxt)
    {
      double value = 0.0;
      stringstream s2(pPreExptxt); s2 >> value ; set_PreExp(value);
    }
    const char* pNInftxt   = ppReac->XmlReadValue("me:nInfinity",false);
    if (pNInftxt)
    {
      double value = 0.0; stringstream s3(pNInftxt); s3 >> value ; set_NInf(value);
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
          << getName();
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
          << getName();
        return false;
      }
    }

    return true ;
  }

}//namespace
