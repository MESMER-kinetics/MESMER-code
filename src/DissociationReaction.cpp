//-------------------------------------------------------------------------------------------
//
// DissociationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the DissociationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "DissociationReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool DissociationReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.

    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      string errorMsg = "Dissociation reaction " + getName() + " has no reactant.";
      meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
      return false;
    }        
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    if(ppReactant2)
    {
      string errorMsg = "Dissociation reaction " + getName() + " has more than one reactant.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
      return false;
    }

    // Save reactant as CollidingMolecule.

    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_rct1 = pColMol;
    } else {
      meErrorLog.ThrowError(__FUNCTION__, string("Dissociating reactant must be a colliding molecule"), obError);
      return false;
    }

    // Read product details. The detail of products may be absent or, may be needed
    // to calculate the microcanonical rates. If there are products, save them as 
    // type Molecule.

    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);
      if (pMol1) {
        m_pdt1 = dynamic_cast<CollidingMolecule*>(pMol1) ;
      } else {
        string errorMsg = "Dissociation reaction" + getName() + " has no products defined.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg, obWarning);
      }

      Molecule* pMol2 = NULL ;
      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      if (ppProduct2) {
        pMol2 = GetMolRef(ppProduct1);
        if (pMol2) {
          m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2) ;
        } else {
          string errorMsg = "Dissociation reaction " + getName() + " has only one product defined.";
          meErrorLog.ThrowError(__FUNCTION__, errorMsg, obWarning);
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
  double DissociationReaction::calcEquilibriumConstant() {

    double Keq(0.0) ;

    return Keq ;
  }

  //
  // Add dissociation reaction terms to collision matrix.
  //
  void DissociationReaction::AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_rct1)] ;
    const int colloptrsize = dynamic_cast<CollidingMolecule*>(m_rct1)->get_colloptrsize();

    if (collisionOperatorCheck){
      ctest << "\nSystem collision operator check before adding " << getName() << " microrates:\n";
      (*CollOptr).showFinalBits(8);
    }

    for ( int i = 0 ; i < colloptrsize; ++i ) {
      int ii(rctLocation + i) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[i] ;   // Forward loss reaction.
    }

    if (collisionOperatorCheck){
      ctest << "\nSystem collision operator check after adding " << getName() << " microrates:\n";
      (*CollOptr).showFinalBits(8);
    }

  }

  // Read parameters requires to determine reaction heats and rates.
  bool DissociationReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    // Read the heat of reaction (if present).
    const char* pHeatRxntxt = ppReac->XmlReadValue("me:HeatOfReaction",false);
    if (pHeatRxntxt){
      stringstream s1(pHeatRxntxt);
      double value = 0.0; s1 >> value ; setHeatOfReaction(value);
    } else { // Calculate heat of reaction.
      double zpe_pdt1 = m_pdt1->get_zpe();
      double zpe_pdt2 = m_pdt2->get_zpe();
      double zpe_rct1 = m_rct1->get_zpe();
      setHeatOfReaction(zpe_pdt1 + zpe_pdt2 - zpe_rct1);
    }

    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential",false);
    if (pPreExptxt)
    {
      double value = 0.0; stringstream s2(pPreExptxt); s2 >> value  ; set_PreExp(value);
    }
    const char* pNInftxt   = ppReac->XmlReadValue("me:nInfinity",false);
    if (pNInftxt)
    {
      double value = 0.0; stringstream s3(pNInftxt); s3 >> value ; set_NInf(value);
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

  //
  // Access microcanonical rate coefficients - cell values are averaged
  // to give grain values. This code is similar to that in Molecule.cpp
  // and this averaging should be done there. SHR 19/Sep/2004.
  //
  bool DissociationReaction::grnAvrgMicroRateCoeffs() {
    // This grain averaging of the microcanonical rate coefficients is 
    // based on the view from the species that is
    // moving in the current reaction toward the opposite species.

    int MaximumGrain = getEnv().MaxGrn;
    int grnSize = static_cast<int>(getEnv().GrainSize);
    // Check if there are enough cells.
    if (grnSize < 1) {
      cerr << "The requested grain size is invalid.";
      exit(1) ;
    }

    // Extract density of states of equilibrium molecule.
    std::vector<double> rctCellDOS;

    m_GrainKfmc.resize(MaximumGrain, 0.);

    m_rct1->getCellDensityOfStates(rctCellDOS);
    rateConstantGrnAvg(MaximumGrain, grnSize, rctCellDOS, m_CellKfmc, m_GrainKfmc);

    //no detailed balance required
    //code for printing dissociation k(E)s follows

    if (getEnv().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }

}//namespace
