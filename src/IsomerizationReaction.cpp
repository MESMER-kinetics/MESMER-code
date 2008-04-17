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
      string errorMsg = "Isomerizaton reaction " + getName() + " has no reactant.";
      meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
      return false;
    }
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    if(ppReactant2)
    {
      string errorMsg = "Isomerizaton reaction " + getName() + " has more than one reactant.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
      return false;
    }

    // Save reactant as CollidingMolecule.

    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_rct1 = pColMol;
    } else {
      meErrorLog.ThrowError(__FUNCTION__, string("Isomer reactant must be a colliding molecule"), obError);
      return false;
    }

    //Read product details.

    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    pMol1 = GetMolRef(ppProduct1);
    if (!pMol1) {
      string errorMsg = "Isomerizaton reaction " + getName() + " has no product.";
      meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
      return false;
    }
    PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
    if(ppProduct2)
    {
      string errorMsg = "Isomerizaton reaction " + getName() + " has more than one product.";
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

    double Qrct1 = m_rct1->grnCanPrtnFn() ;
    double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

    // Calculate the equilibrium constant.
    {stringstream errorMsg;
    errorMsg << "Qrct1 = " << Qrct1 << ", Qpdt1 = " << Qpdt1 ;
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

    double beta = getEnv().beta ;

    Keq = (Qpdt1 / Qrct1)*exp(-beta * getHeatOfReaction()) ;

    return Keq ;
  }

  //
  // Add isomer reaction terms to collision matrix.
  //
  void IsomerizationReaction::AddReactionTerms(dMatrix         *CollOptr,
    isomerMap       &isomermap,
    const double    rMeanOmega)
  {
    // Locate isomers in system matrix.
    const int rctLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_rct1)] ;
    const int pdtLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_pdt1)] ;

    // Get densities of states for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn;
    vector<double> rctDOS;
    vector<double> pdtDOS;
    m_rct1->getGrainDensityOfStates(rctDOS) ;
    m_pdt1->getGrainDensityOfStates(pdtDOS) ;

    const int idx = m_pdt1->get_grnZpe() - m_rct1->get_grnZpe() ;
    for ( int i = max(0,-idx) ; i < min(MaximumGrain,(MaximumGrain-idx)) ; ++i ) {
      int ll = i + idx ;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation + i) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[ll] ;                            // Forward loss reaction.
      (*CollOptr)[jj][jj] -= rMeanOmega * m_GrainKfmc[ll]*rctDOS[ll]/pdtDOS[i] ;       // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(rctDOS[ll]/pdtDOS[i]) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                     // Reactive gain.
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
      double zpe_pdt1 = m_pdt1->get_zpe();
      double zpe_rct1 = m_rct1->get_zpe();
      setHeatOfReaction(zpe_pdt1 - zpe_rct1);
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
  bool IsomerizationReaction::grnAvrgMicroRateCoeffs() {
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
    m_GrainKbmc.resize(MaximumGrain, 0.);

    m_rct1->getCellDensityOfStates(rctCellDOS);
    rateConstantGrnAvg(MaximumGrain, grnSize, rctCellDOS, m_CellKfmc, m_GrainKfmc);
    grainRateCoeffDetailedBalance(1);

    // the following is for printing of f & r k(E)s

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

}//namespace
