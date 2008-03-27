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
        /* need to create one (mark _2007_12_10__17_10_18_)
        write a SuperMolecule creator that acquire a position in the XML
        */
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
    double Keq = calcEquilibriumConstant() ;
    const double excess = getExcessReactantConc();

    vector<double> rctDOS;
    vector<double> pdtDOS;
    m_srct->grnDensityOfStates(rctDOS) ;
    m_pdt1->grnDensityOfStates(pdtDOS) ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn ;
    vector<double> adductBoltz ;
    m_pdt1->grnBoltzDist(adductBoltz) ;

    const int colloptrsize = dynamic_cast<CollidingMolecule*>(m_pdt1)->get_colloptrsize() ; // collision operator size of the isomer

    const int grnWellBot = MaximumGrain - colloptrsize; // well bottom of the isomer

    double DissRateCoeff(0.0) ;
    if (debugFlag) ctest << "colloptrsize = " << colloptrsize << ", grnWellBot = " << grnWellBot << endl;
    if (debugFlag) ctest << "Keq * excess = " << Keq << " * " << excess << " = " << Keq * excess << endl;

    // Multiply equilibrum constant by concentration of excess reactant.
    // concentration of excess reactant should be in molec/cm3. This gives
    // a dimensionless pseudo-isomerization equilibrium constant.
    Keq *= excess ;

    if (collisionOperatorCheck){
      ctest << "\nSystem collision operator check before adding " << getName() << " microrates:\n";
      (*CollOptr).showFinalBits(8);
    }

    for ( int i = 0; i < colloptrsize; ++i) {
      int ll = i + grnWellBot ;
      int ii(pdtLoc + i) ;

      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKbmc[ll] ;                           // Backward loss of the adduct
      if (i < colloptrsize && i > (colloptrsize - 8) && debugFlag) {
        ctest << rMeanOmega << " * " << "m_GrainKbmc[" << ll
        << "](" << m_GrainKbmc[ll] << ") = " << rMeanOmega * m_GrainKbmc[ll] << endl;
      }
      (*CollOptr)[jj][ii]  = rMeanOmega * m_GrainKbmc[ll] * sqrt(adductBoltz[ll] * Keq) ; // Reactive gain of the source
      (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                                    // Reactive gain
      DissRateCoeff       += rMeanOmega * m_GrainKbmc[ll] * adductBoltz[ll] ;
    }

    if (debugFlag) ctest << "DissRateCoeff = " << DissRateCoeff << endl;
    (*CollOptr)[jj][jj] -= (DissRateCoeff * Keq);       // Backward loss reaction from detailed balance.

    if (collisionOperatorCheck){
      ctest << "\nSystem collision operator check after adding " << getName() << " microrates:\n";
      (*CollOptr).showFinalBits(8);
    }

  }

  //
  // Calculate reaction equilibrium constant for the general reaction
  //        A + B  <===> C
  //
  double AssociationReaction::calcEquilibriumConstant() {

    // equilibrium constant:
    long double Keq(0.0) ;

    // Get Canonical partition functions.
    long double Qrcts = m_srct->grnCanPrtnFn();
    long double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

    long double mass_rct1 = m_rct1->getMass();
    long double mass_rct2 = m_rct2->getMass();
    long double mass_srct = mass_rct1 + mass_rct2;

    // Calculate the equilibrium constant.
    long double beta = getEnv().beta ;

    Keq = Qrcts/Qpdt1 ;
    if(debugFlag) ctest << "Keq = " << Keq << endl;

    /* Electronic degeneracies were already accounted for in DOS calculations */
    // Heat of reaction
    long double _expon = -beta * getHeatOfReaction() * kJPerMolInRC;
    Keq /= exp(_expon) ;
    if(debugFlag) ctest << "Keq = " << Keq 
      << ", getHeatOfReaction() = " << getHeatOfReaction() 
      << ", beta = " << beta 
      << ", _expon = " << _expon 
      << ", dh = " << getHeatOfReaction() * kJPerMolInRC 
      << ", kJPerMolInRC = " << kJPerMolInRC << endl;

    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tau = 2.0593e19 * pow(2. * M_PI ,1.5);
    Keq *= (tp_C * pow(mass_rct1 * mass_rct2 / (mass_srct * beta), 1.5l));
    if(debugFlag) ctest << "Keq = " << Keq << endl;

    Keq = 1./Keq;
    if(debugFlag) ctest << "Keq = " << Keq << endl;
    return (double) Keq ;

  }

}//namespace
