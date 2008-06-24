//-------------------------------------------------------------------------------------------
//
// IsomerizationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IsomerizationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IsomerizationReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool IsomerizationReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      cerr << "Cannot get reactant definition for Isomerization reaction " << getName() << ".";
      return false;
    }

    // Save reactant as CollidingMolecule.
    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_rct1 = pColMol;
    } else {
      cerr << "Isomer reactant must be a colliding molecule";
      return false;
    }

    //Read product details.
    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    pMol1 = GetMolRef(ppProduct1);
    if (!pMol1) {
      cerr << "Cannot get product definition for Isomerization reaction " << getName() << ".";
      return false;
    }

    // Save product as CollidingMolecule.
    pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_pdt1 = pColMol;
    } else {
      cerr << "Isomer product must be a colliding molecule";
      return false;
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
  double IsomerizationReaction::calcEquilibriumConstant() {

    double Keq(0.0) ;

    // Get Canonical partition functions.
    double Qrct1 = m_rct1->rovibronicGrnCanPrtnFn() ;
    double Qpdt1 = m_pdt1->rovibronicGrnCanPrtnFn() ;

    double beta = getEnv().beta ;

    double HeatOfReaction = getHeatOfReaction();
    Keq = (Qpdt1 / Qrct1)*exp(-beta * HeatOfReaction) ;

    return Keq ;
  }

  //
  // Add isomer reaction terms to collision matrix.
  //
  void IsomerizationReaction::AddReactionTerms(dMatrix         *CollOptr,
    isomerMap       &isomermap,
    const double    rMeanOmega)
  {
    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    vector<double> pdtDOS;
    m_rct1->getGrainDensityOfStates(rctDOS) ;
    m_pdt1->getGrainDensityOfStates(pdtDOS) ;

    // Locate isomers in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    const int pdtLocation = isomermap[m_pdt1] ;

    const int TSgrainZPE  = getTSFluxGrnZPE();
    const int rctGrainZPE = m_rct1->get_grnZpe();
    const int pdtGrainZPE = m_pdt1->get_grnZpe();
    const int MaximumGrain = getEnv().MaxGrn;

    for ( int i = 0 ; i < MaximumGrain - TSgrainZPE ; ++i ) {
      int ll = i + TSgrainZPE - rctGrainZPE;
      int mm = i + TSgrainZPE - pdtGrainZPE;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation + mm) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainTSFlux[i] / rctDOS[ll];                     // Forward loss reaction.
      (*CollOptr)[jj][jj] -= rMeanOmega * m_GrainTSFlux[i] / pdtDOS[mm] ;                    // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = rMeanOmega * m_GrainTSFlux[i] / sqrt(rctDOS[ll] * pdtDOS[mm]) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                           // Reactive gain.
    }
  }

  // Read parameters requires to determine reaction heats and rates.
  bool IsomerizationReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    // Read the heat of reaction (if present).
    const char* pHeatRxntxt = ppReac->XmlReadValue("me:HeatOfReaction",false);
    if (pHeatRxntxt){
      stringstream s1(pHeatRxntxt);
      double value = 0.0; s1 >> value ; setHeatOfReaction(value);
    } else { // Calculate heat of reaction.
      setHeatOfReaction(get_relative_pdtZPE(), get_relative_rctZPE());
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
  // Calculate grained forward and reverse k(E)s from trainsition state flux
  //
  void IsomerizationReaction::calcGrainRateCoeffs(){
    vector<double> rctGrainDOS;
    vector<double> pdtGrainDOS;
    m_rct1->getGrainDensityOfStates(rctGrainDOS) ;
    m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;

    const int TSgrainZPE  = getTSFluxGrnZPE();
    const int rctGrainZPE = m_rct1->get_grnZpe();
    const int pdtGrainZPE = m_pdt1->get_grnZpe();

    const int MaximumGrain = getEnv().MaxGrn;
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain , 0.0);
    m_GrainKbmc.clear();
    m_GrainKbmc.resize(MaximumGrain , 0.0);

    // For AssociationReaction, TSFlux is calculated from ILT
    for (int i = TSgrainZPE - pdtGrainZPE, j = 0; i < MaximumGrain; ++i, ++j){
      m_GrainKbmc[i] = m_GrainTSFlux[j] / pdtGrainDOS[i];
    }
    for (int i = TSgrainZPE - rctGrainZPE, j = 0; i < MaximumGrain; ++i, ++j){
      m_GrainKfmc[i] = m_GrainTSFlux[j] / rctGrainDOS[i];
    }

    // the code that follows is for printing of the f & r k(E)s
    if (getEnv().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getEnv().kbEGrainsEnabled){
      ctest << "\nk_b(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKbmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getEnv().testRateConstantEnabled)
      testRateConstant();
  }

  // Test k(T)
  void IsomerizationReaction::testRateConstant() const {
  
  }
}//namespace
