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

  // Is reaction equilibrating and therefore contributes
  // to the calculation of equilibrium fractions.
  bool IsomerizationReaction::isEquilibratingReaction(double &Keq, ModelledMolecule **rct, ModelledMolecule **pdt) {

    Keq = calcEquilibriumConstant() ;

    *rct = m_rct1 ;
    *pdt = m_pdt1 ;

    return true ;
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
  void IsomerizationReaction::AddReactionTerms(qdMatrix         *CollOptr,
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

    const int colloptrsize = m_pdt1->get_colloptrsize();
    const int forwardThreshE = get_effectiveForwardTSFluxGrnZPE();
    const int reverseThreshE = get_effectiveReverseTSFluxGrnZPE();
    const int fluxStartIdx = get_TSFluxStartIdx();

    for ( int i=fluxStartIdx, j = reverseThreshE, k=0; j < colloptrsize; ++i, ++j, ++k) {
      int ll = k + forwardThreshE;
      int mm = k + reverseThreshE;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation + mm) ;
      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * m_GrainTSFlux[i] / rctDOS[ll]);                     // Forward loss reaction.
      (*CollOptr)[jj][jj] -= qd_real(rMeanOmega * m_GrainTSFlux[i] / pdtDOS[mm]) ;                    // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = qd_real(rMeanOmega * m_GrainTSFlux[i] / sqrt(rctDOS[ll] * pdtDOS[mm])) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                           // Reactive gain.
    }
 
}

// Read parameters requires to determine reaction heats and rates.
bool IsomerizationReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy", false);
    if (pActEnetxt)
    {
      PersistPtr ppActEne = ppReac->XmlMoveTo("me:activationEnergy") ;
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput;
      if (unitsTxt){
        unitsInput = unitsTxt;
      }
      else{
        unitsInput = "kJ/mol";
      }
      const char* pLowertxt = ppActEne->XmlReadValue("lower", false);
      const char* pUppertxt = ppActEne->XmlReadValue("upper", false);
      const char* pStepStxt = ppActEne->XmlReadValue("stepsize", false);
      double value(getConvertedEnergy(unitsInput, tmpvalue));
      if (pLowertxt && pUppertxt){
        double tmpvalueL(0.0), tmpvalueU(0.0), tmpstepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> tmpvalueL; s4 >> tmpvalueU; s5 >> tmpstepsize;
        double valueL(getConvertedEnergy(unitsInput, tmpvalueL));
        double valueU(getConvertedEnergy(unitsInput, tmpvalueL));
        double stepsize(getConvertedEnergy(unitsInput, tmpstepsize));
        set_EInf(valueL, valueU, stepsize);
      }
      else{
        set_EInf(value);
      }
    }

    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential", false);
    if (pPreExptxt)
    {
      PersistPtr ppPreExponential = ppReac->XmlMoveTo("me:preExponential") ;
      double value = 0.0;
      stringstream s2(pPreExptxt); s2 >> value ;
      const char* pLowertxt = ppPreExponential->XmlReadValue("lower", false);
      const char* pUppertxt = ppPreExponential->XmlReadValue("upper", false);
      const char* pStepStxt = ppPreExponential->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        set_PreExp(valueL, valueU, stepsize);
      }
      else{
        set_PreExp(value);
      }
    }

    const char* pNInftxt = ppReac->XmlReadValue("me:nInfinity", false);
    if (pNInftxt)
    {
      PersistPtr ppNInf = ppReac->XmlMoveTo("me:nInfinity") ;
      double value = 0.0;
      stringstream s2(pNInftxt); s2 >> value ;
      const char* pLowertxt = ppNInf->XmlReadValue("lower", false);
      const char* pUppertxt = ppNInf->XmlReadValue("upper", false);
      const char* pStepStxt = ppNInf->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        set_NInf(valueL, valueU, stepsize);
      }
      else{
        set_NInf(value);
      }
    }

    const char* pTInftxt = ppReac->XmlReadValue("me:TInfinity");  // read XML designation
    if (!pTInftxt){           // if there is no Tinfinity
      cinfo << "Tinfinity has not been specified; set to the default value of 298 K" << endl;
    }
    else{                     // if there is a Tinfinity
      double TInf(298.0);     // set Tinf to 298
      stringstream s3(pTInftxt);    // initialize stringstream
      s3 >> TInf;             // use the stringstream to get Tinfinty from the input file
      if(TInf <= 0){          // if Tinf is <= 0, set to 298
        cinfo << "Tinfinity is less than or equal to 0; set to the default value of 298 K";
      }
      else
        set_TInf(TInf);         // else set Tinf to the value in the input
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
// Calculate grained forward and reverse k(E)s from transition state flux
//
void IsomerizationReaction::calcGrainRateCoeffs(){
  vector<double> rctGrainDOS;
  vector<double> pdtGrainDOS;
  m_rct1->getGrainDensityOfStates(rctGrainDOS) ;
  m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;

  calculateEffectiveGrainedThreshEn();
  const int forwardTE = get_effectiveForwardTSFluxGrnZPE();
  int reverseTE = get_effectiveReverseTSFluxGrnZPE();
  calculateTSfluxStartIdx();
  const int fluxStartIdx = get_TSFluxStartIdx();

  const int MaximumGrain = (getEnv().MaxGrn-fluxStartIdx);
  m_GrainKfmc.clear();
  m_GrainKfmc.resize(MaximumGrain , 0.0);
  m_GrainKbmc.clear();
  m_GrainKbmc.resize(MaximumGrain , 0.0);

  for (int i = reverseTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
    m_GrainKbmc[i] = m_GrainTSFlux[j] / pdtGrainDOS[i];
  }
  for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
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
void IsomerizationReaction::testRateConstant() {

  double k_forward(0.0), k_backward(0.0);
  vector<double> rctGrainDOS, rctGrainEne, pdtGrainDOS, pdtGrainEne ;
  m_rct1->getGrainDensityOfStates(rctGrainDOS);
  m_pdt1->getGrainDensityOfStates(pdtGrainDOS);
  m_rct1->getGrainEnergies(rctGrainEne);
  m_pdt1->getGrainEnergies(pdtGrainEne);
  const int MaximumGrain = (getEnv().MaxGrn-get_TSFluxStartIdx());
  const double beta = getEnv().beta;
  const double temperature = 1. / (boltzmann_RCpK * beta);

  for(int i(0); i < MaximumGrain; ++i){
    k_forward += m_GrainKfmc[i] * exp( log(rctGrainDOS[i]) - beta * rctGrainEne[i]);
    k_backward += m_GrainKbmc[i] * exp( log(pdtGrainDOS[i]) - beta * pdtGrainEne[i]);
  }

  const double rctprtfn = canonicalPartitionFunction(rctGrainDOS, rctGrainEne, beta);
  const double pdtprtfn = canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, beta);
  k_forward /= rctprtfn;
  k_backward /= pdtprtfn;
  set_forwardCanonicalRateCoefficient(k_forward);
  set_backwardCanonicalRateCoefficient(k_backward);

  ctest << endl << "Canonical pseudo first order forward rate constant of isomerization reaction " 
    << getName() << " = " << get_forwardCanonicalRateCoefficient() << " s-1 (" << temperature << " K)" << endl;
  ctest << "Canonical pseudo first order backward rate constant of isomerization reaction " 
    << getName() << " = " << get_backwardCanonicalRateCoefficient() << " s-1 (" << temperature << " K)" << endl;
}

void IsomerizationReaction::calculateEffectiveGrainedThreshEn(void){  // see the comments in
  double thresh = get_ThresholdEnergy();    // calculateEffectiveGrainedThreshEn under AssociationReaction.cpp
  double RxnHeat = getHeatOfReaction();
  int TS_en = this->getTSFluxGrnZPE();
  int pdt_en = m_pdt1->get_grnZpe();
  int rct_en = m_rct1->get_grnZpe();
  int GrainedRxnHeat = pdt_en - rct_en;
  if(thresh<0.0){
    set_effectiveForwardTSFluxGrnZPE(0);
    set_effectiveReverseTSFluxGrnZPE(-GrainedRxnHeat);
  }
  else if(thresh>0.0 && thresh<RxnHeat){
    set_effectiveForwardTSFluxGrnZPE(GrainedRxnHeat);
    set_effectiveReverseTSFluxGrnZPE(0);
  }
  else{
    set_effectiveForwardTSFluxGrnZPE(TS_en-rct_en);
    set_effectiveReverseTSFluxGrnZPE(TS_en-pdt_en);
  }
}

}//namespace
