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
      m_rct1 = tmp_rct2;
      m_rct2 = tmp_rct1;
    }

    if(!m_rct1){
      cerr << "the deficient reactant in the association reaction is undefined" << endl;
      return false;
    }
    if(!m_rct2){
      cerr << "the excess reactant in the association reaction is undefined" << endl;
      return false;
    }

    //Read product details.

    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    pMol1 = GetMolRef(ppProduct1);
    if (!pMol1) {
      cerr << "Cannot find product molecule definition for association reaction " << getName() << ".";
      return false;
    }

    // Save product as CollidingMolecule.

    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol){
      m_pdt1 = pColMol;
    } else {
      cerr << "Isomer product must be a colliding molecule";
      return false;
    }

    // Read the transition state (if present)
    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if (ppTransitionState)
    {
      TransitionState* pTrans = dynamic_cast<TransitionState*>(GetMolRef(ppTransitionState,"transitionState"));
      if(pTrans) m_TransitionState = pTrans;
    }

    // Read heat of reaction and rate parameters.

    return ReadRateCoeffParameters(ppReac) ;
  }

  //
  // Add (REVERSIBLE) association reaction terms to collision matrix.
  //
  void AssociationReaction::AddReactionTerms(qdMatrix      *CollOptr,
    isomerMap    &isomermap,
    const double rMeanOmega)
  {
    // Get densities of states of the adduct for detailed balance.
    vector<double> pdtDOS;
    m_pdt1->getGrainDensityOfStates(pdtDOS) ;

    // Locate isomers in system matrix.
    const int pdtLoc =      isomermap[m_pdt1] ;
    const int jj     = (*m_sourceMap)[get_pseudoIsomer()] ;

    // Get equilibrium constant.
    const double Keq = calcEquilibriumConstant() ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn ;
    vector<double> adductPopFrac ; // Population fraction of the adduct

    m_pdt1->normalizedBoltzmannDistribution(adductPopFrac, MaximumGrain) ;
    qd_real DissRateCoeff(0.0) ;

    const int colloptrsize = m_pdt1->get_colloptrsize();
    const int forwardThreshE = get_effectiveForwardTSFluxGrnZPE();
    const int reverseThreshE = get_effectiveReverseTSFluxGrnZPE();
    const int fluxStartIdx = get_TSFluxStartIdx();

    for ( int i = reverseThreshE, j = fluxStartIdx; i < colloptrsize; ++i, ++j) {
      int ii(pdtLoc + i) ;

      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * m_GrainTSFlux[j] / pdtDOS[i]);                                // Loss of the adduct to the source
      (*CollOptr)[jj][ii]  = qd_real(rMeanOmega * m_GrainTSFlux[j] * sqrt(adductPopFrac[i] * Keq) / pdtDOS[i]);// Reactive gain of the source
      (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                                                      // Reactive gain (symmetrization)
      DissRateCoeff       += qd_real(m_GrainTSFlux[j] * adductPopFrac[i] / pdtDOS[i]);
    }
    (*CollOptr)[jj][jj] -= qd_real(rMeanOmega * DissRateCoeff * Keq);       // Loss of the source from detailed balance.
  }

  //
  // Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
  //
  double AssociationReaction::rctsRovibronicGrnCanPrtnFn() {
    vector<double> rctGrainDOS;
    vector<double> rctGrainEne;
    calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

    // Calculate the rovibronic partition function based on the grain DOS
    // The following catches the case where the molecule is a single atom
    double CanPrtnFn = max(canonicalPartitionFunction(rctGrainDOS, rctGrainEne, getEnv().beta), 1.0) ;
    if (CanPrtnFn == 1.0){
      // Electronic partition function for atom is accounted here.
      CanPrtnFn = double(m_rct1->getSpinMultiplicity() * m_rct2->getSpinMultiplicity()) ;
    }

    return CanPrtnFn ;
  }

  // Is reaction equilibrating and therefore contributes
  // to the calculation of equilibrium fractions.
  bool AssociationReaction::isEquilibratingReaction(double &Keq, ModelledMolecule **rct, ModelledMolecule **pdt) {

    Keq = calcEquilibriumConstant() ;

    *rct = m_rct1 ;
    *pdt = m_pdt1 ;

    return true ;
  }

  //
  // Calculate reaction equilibrium constant for the general reaction
  //        A + B  <===> C
  //
  double AssociationReaction::calcEquilibriumConstant() {

    // equilibrium constant:
    double Keq(0.0) ;
    const double beta = getEnv().beta ;

    // partition function for each reactant
    double Qrcts = rctsRovibronicGrnCanPrtnFn();

    // rovibronic partition function for reactants multiplied by translation contribution
    Qrcts *= translationalContribution(m_rct1->getMass(), m_rct2->getMass(), beta);

    // rovibronic partition function for product
    const double Qpdt1 = m_pdt1->rovibronicGrnCanPrtnFn() ;

    Keq = Qpdt1 / Qrcts;

    // Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells
    const double HeatOfReaction = getHeatOfReaction() ;
    const double _expon = -beta * HeatOfReaction;
    Keq *= exp(_expon) ;

    Keq *= m_ERConc ;
    //
    // K_eq = ( [C]/[A][B] ) * [A] = [C]/[B]
    //
    // where [A] is the reactant what is in excess (seen as constant).
    // Therefore, the K_eq here is essentially the pseudo-first-order equilibrium constant.

    return Keq ;

  }

  // Read parameters requires to determine reaction heats and rates.
  bool AssociationReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

    const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy", false);
    if (pActEnetxt)
    {
      PersistPtr ppActEne = ppReac->XmlMoveTo("me:activationEnergy") ;
      double value = 0.0;
      stringstream s2(pActEnetxt); s2 >> value ;
      const char* pLowertxt = ppActEne->XmlReadValue("lower", false);
      const char* pUppertxt = ppActEne->XmlReadValue("upper", false);
      const char* pStepStxt = ppActEne->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        set_ThresholdEnergy(value, valueL, valueU, stepsize);
      }
      else{
        set_ThresholdEnergy(value);
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

    const char* pERConctxt = ppReac->XmlReadValue("me:excessReactantConc",false);
    if (!pERConctxt){
      cerr << "Concentration of excess reactant has not been specified.";
      return false;
    } else {
      stringstream s3(pERConctxt) ;
      s3 >> m_ERConc ;
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
    else{
      cinfo << "No tunneling method was found for " << getName() << endl;
    }

    return true ;
  }

  //
  // Calculate grained forward and reverse k(E)s from transition state flux
  //
  void AssociationReaction::calcGrainRateCoeffs(){

    vector<double> rctGrainDOS;
    vector<double> rctGrainEne;
    vector<double> pdtGrainDOS;
    m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;
    calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

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

  void AssociationReaction::testRateConstant() {
    vector<double> pdtGrainDOS;
    vector<double> pdtGrainEne;

    m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;
    m_pdt1->getGrainEnergies(pdtGrainEne);

    const int MaximumGrain = (getEnv().MaxGrn-get_TSFluxStartIdx());
    const double beta = getEnv().beta;
    for (int i = 0; i < MaximumGrain; ++i){
      m_backwardCanonicalRate += m_GrainKbmc[i] * exp( log(pdtGrainDOS[i]) - beta * pdtGrainEne[i] ) ;
    }
    const double prtfn = canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, beta);
    m_backwardCanonicalRate /= prtfn;
    double Keq = calcEquilibriumConstant();
    m_forwardCanonicalRate = m_backwardCanonicalRate * Keq;
    const double temperature = 1. / (boltzmann_RCpK * beta);
    ctest << endl << "Canonical pseudo first order rate constant of association reaction "
      << getName() << " = " << m_forwardCanonicalRate << " s-1 (" << temperature << " K)" << endl;
    ctest << "Canonical bimolecular rate constant of association reaction "
      << getName() << " = " << m_forwardCanonicalRate/m_ERConc << " cm^3/mol/s (" << temperature << " K)" << endl;
    ctest << "Canonical first order rate constant for the reverse of reaction "
      << getName() << " = " << m_backwardCanonicalRate << " s-1 (" << temperature << " K)" << endl;
  }

  //
  // Calculate the rovibrational density of states of reactants.
  //
  bool AssociationReaction::calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)    // Calculate rovibrational reactant DOS
  {
    std::vector<double> rctsCellDOS;
    getRctsCellDensityOfStates(rctsCellDOS);

    std::vector<double> shiftedCellDOS;
    std::vector<double> shiftedCellEne;
    const int MaximumCell = getEnv().MaxCell;
    const int cellOffset = get_pseudoIsomer()->get_cellOffset();
    std::vector<double> rctsCellEne;
    getCellEnergies(MaximumCell, rctsCellEne);
    shiftCells(MaximumCell, cellOffset, rctsCellDOS, rctsCellEne, shiftedCellDOS, shiftedCellEne);

    string catName = m_rct1->getName() + " + " + m_rct2->getName();
    calcGrainAverages(getEnv().MaxGrn, getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, grainDOS, grainEne, catName);

    return true;
  }

  //
  // Get reactants cell density of states.
  //
  void AssociationReaction::getRctsCellDensityOfStates(vector<double> &cellDOS) {
    get_rctsDensityOfStatesCalculator()->countDimerCellDOS(m_rct1, m_rct2, cellDOS);
  }

  const int AssociationReaction::get_rctsGrnZpe(){
    double grnZpe = (m_rct1->get_zpe()+m_rct2->get_zpe()-getEnv().EMin) / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void AssociationReaction::calculateEffectiveGrainedThreshEn(void){  // calculate the effective forward and reverse
    double thresh = get_ThresholdEnergy();  // threshold energy for an association reaction
    double RxnHeat = getHeatOfReaction();
    int TS_en = this->getTSFluxGrnZPE();
    int pdt_en = m_pdt1->get_grnZpe();
    int rct_en = get_rctsGrnZpe();
    int GrainedRxnHeat = pdt_en - rct_en;
    if(thresh<0.0){                          // if the forward threshold energy is negative
      set_effectiveForwardTSFluxGrnZPE(0);   // forward grained flux threshold energy = 0
      set_effectiveReverseTSFluxGrnZPE(-GrainedRxnHeat);  // reverse grained flux threshold energy = heat of reaction
    }
    else if(thresh>0.0 && thresh<RxnHeat){   // if the reverse threshold energy is negative
      set_effectiveForwardTSFluxGrnZPE(GrainedRxnHeat);  // forward grained flux threshold energy = heat of reaction
      set_effectiveReverseTSFluxGrnZPE(0);   // reverse grained flux threshold energy = 0
    }
    else{
      set_effectiveForwardTSFluxGrnZPE(TS_en-rct_en);  // forward grained flux threshold energy = TS energy - rct energy
      set_effectiveReverseTSFluxGrnZPE(TS_en-pdt_en);  // reverse grained flux threshold energy = TS energy - pdt energy
    }
  }

}//namespace
