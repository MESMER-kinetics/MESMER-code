//-------------------------------------------------------------------------------------------
//
// ExchangeReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the ExchangeReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "ExchangeReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool ExchangeReaction::InitializeReaction(PersistPtr ppReac)
  {
        m_ppPersist = ppReac;

        // Read reactant details.

        PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
        Molecule* pMol1 = GetMolRef(ppReactant1);
        if(!pMol1){
            cerr << "Cannot find 1st reactant molecule definition for association reaction " << getName() << ".";
            return false;
        }
        PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        Molecule* pMol2 = GetMolRef(ppReactant2);
        if(!pMol2)
        {
          cerr << "Cannot find 2nd reactant molecule definition for association reaction " << getName() << ".";
          return false;
        }

        // if deficientReactantLocation=true, then pMol1 (the first rct
        // in the XML input) is the deficient reactant (m_rct1)

        ModelledMolecule* tmp_rct1 = dynamic_cast<ModelledMolecule*>(pMol1);
        ModelledMolecule* tmp_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);

        if(deficientReactantLocation){
          m_rct1 = tmp_rct1;
          m_rct2 = tmp_rct2;
        }
        else {
          m_rct1 = tmp_rct1;
          m_rct2 = tmp_rct2;
        }

        if(!m_rct1){
          cerr << "the deficient reactant in the association reaction is undefined" << endl;
          return false;
        }
        if(!m_rct2){
          cerr << "the excess reactant in the association reaction is undefined" << endl;
          return false;
        }

    // Read product details. Save them as type Molecule.

    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);
      if (pMol1) {
        m_pdt1 = dynamic_cast<ModelledMolecule*>(pMol1) ;
      } else {
        cerr << "Exchange reaction" << getName() << " has no products defined.";
      }

      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      if (ppProduct2) {
        pMol2 = GetMolRef(ppProduct1);
        if (pMol2) {
          m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2) ;
        } else {
          cerr << "Exchange reaction " << getName() << " has only one product defined.";
        }
      }
    }

    // Read the transition state (if present).

    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if (ppTransitionState)
    {
      TransitionState* pTrans = dynamic_cast<TransitionState*>(GetMolRef(ppTransitionState));
      if(pTrans)
        m_TransitionState = pTrans;
    }

    // Read heat of reaction and rate parameters.

    return ReadRateCoeffParameters(ppReac) ;
  }

  //
  // Calculate reaction equilibrium constant.
  //
  double ExchangeReaction::calcEquilibriumConstant() {

    double Keq(0.0) ;

    // Get Canonical partition functions.

    const double Qrct1 = m_rct1->rovibronicGrnCanPrtnFn() ;
    const double Qrct2 = m_rct2->rovibronicGrnCanPrtnFn() ;
    const double Qpdt1 = m_pdt1->rovibronicGrnCanPrtnFn() ;
    const double Qpdt2 = m_pdt2->rovibronicGrnCanPrtnFn() ;

    const double mass_rct1 = m_rct1->getMass() ;
    const double mass_rct2 = m_rct2->getMass() ;
    const double mass_pdt1 = m_pdt1->getMass() ;
    const double mass_pdt2 = m_pdt2->getMass() ;

    // Calculate the equilibrium constant.
    double beta = getEnv().beta ;

    Keq = (Qpdt1 * Qpdt2)/(Qrct1 * Qrct2) ;

    // Electronic degeneracies were already accounted for in DOS calculations.
    // Heat of reaction
    double HeatOfReaction = getHeatOfReaction();
    Keq *= exp(-beta * HeatOfReaction) ;

    // Translational contribution

    Keq *= (pow(mass_pdt1 * mass_pdt2 / (mass_rct1 * mass_rct2), 1.5));

    return Keq ;
  }

  //
  // Add exchange reaction terms to collision matrix.
  //
  void ExchangeReaction::AddReactionTerms(qdMatrix      *CollOptr,
    isomerMap    &isomermap,
    const double rMeanOmega)
  {

  }

  // Read parameters requires to determine reaction heats and rates.
  bool ExchangeReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    return true ;
  }

  //
  // Calculate grained forward and reverse k(E)s from trainsition state flux
  //
  void ExchangeReaction::calcGrainRateCoeffs(){
    if (getEnv().testRateConstantEnabled)
      testRateConstant();
  }
  
  // Test k(T)
  void ExchangeReaction::testRateConstant() {
  
  }

}//namespace
