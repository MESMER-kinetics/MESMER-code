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
  // One fact to know is that whatever happens in the reaction operator in this routine is in the unit of 
  // "per collision". In addition, the expression of every entry has to be first similarly transformed to the 
  // symmetrized matrix corrdinates using detailed balance just like the way of constructing collision operator.
  // The flux has to be divided by omega before putting into the entry because flux is calculated in the unit of
  // second, we need to convert the flux into the unit of per collision.
  // 
  bool AssociationReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if(!ppReactantList)
      ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
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

    Molecule* tmp_rct1 = pMol1;
    Molecule* tmp_rct2 = pMol2;

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
    PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
    if(!ppProductList)
      ppProductList=ppReac; //Be forgiving; we can get by without a productList element

    PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
    pMol1 = GetMolRef(ppProduct1);
    if (!pMol1) {
      cerr << "Cannot find product molecule definition for association reaction " << getName() << ".";
      return false;
    }

    // Save product as Molecule.

    Molecule* pColMol = pMol1;
    if(pColMol){
      m_pdt1 = pColMol;
    } else {
      cerr << "Isomer product must be a colliding molecule";
      return false;
    }

    // Read heat of reaction and rate parameters.
    return ReadRateCoeffParameters(ppReac) ;

  }

  void AssociationReaction::AddReactionTermsWithReservoirState(qdMatrix      *CollOptr,
    molMapType    &isomermap,
    const double rMeanOmega)
  {
    // Get densities of states of the adduct for detailed balance.
    vector<double> pdtDOS;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;

    // Locate isomers in system matrix.
    const int pdtLoc =      isomermap[m_pdt1] ;
    const int jj     = (*m_sourceMap)[get_pseudoIsomer()] ;

    // Get equilibrium constant.
    const double Keq = calcEquilibriumConstant() ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn ;
    vector<double> adductPopFrac ; // Population fraction of the adduct
    m_pdt1->getColl().normalizedGrnBoltzmannDistribution(adductPopFrac, MaximumGrain) ;

    qd_real DissRateCoeff(0.0) ;

    const int pNGG(m_pdt1->getColl().getNumberOfGroupedGrains());
    const int pShiftedGrains(pNGG == 0 ? 0 : pNGG - 1);
    const int colloptrsize = m_pdt1->getColl().get_colloptrsize();
    //const int forwardThreshE = get_EffGrnFwdThreshold();
    const int reverseThreshE = get_EffGrnRvsThreshold();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    // Note: reverseThreshE will always be greater than pNGG here

    for ( int i = reverseThreshE, j = fluxStartIdx; i < colloptrsize; ++i, ++j) {
      int ii(pdtLoc + i - pShiftedGrains) ;

      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * m_GrainFlux[j] / pdtDOS[i]);                                // Loss of the adduct to the source
      (*CollOptr)[jj][ii]  = qd_real(rMeanOmega * m_GrainFlux[j] * sqrt(adductPopFrac[i] * Keq) / pdtDOS[i]);// Reactive gain of the source
      (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                                                      // Reactive gain (symmetrization)
      DissRateCoeff       += qd_real(m_GrainFlux[j] * adductPopFrac[i] / pdtDOS[i]);
    }
    (*CollOptr)[jj][jj] -= qd_real(rMeanOmega * DissRateCoeff * Keq);       // Loss of the source from detailed balance.
  }

  void AssociationReaction::AddReactionTerms(qdMatrix      *CollOptr,
    molMapType    &isomermap,
    const double rMeanOmega)
  {
    // Get densities of states of the adduct for detailed balance.
    vector<double> pdtDOS;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;

    // Locate isomers in system matrix.
    const int pdtLoc =      isomermap[m_pdt1] ;
    const int jj     = (*m_sourceMap)[get_pseudoIsomer()] ;

    // Get equilibrium constant.
    const double Keq = calcEquilibriumConstant() ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn ;
    vector<double> adductPopFrac ; // Population fraction of the adduct

    m_pdt1->getColl().normalizedGrnBoltzmannDistribution(adductPopFrac, MaximumGrain) ;
    qd_real DissRateCoeff(0.0) ;

    const int colloptrsize = m_pdt1->getColl().get_colloptrsize();
    //const int forwardThreshE = get_EffGrnFwdThreshold();
    const int reverseThreshE = get_EffGrnRvsThreshold();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    for ( int i = reverseThreshE, j = fluxStartIdx; i < colloptrsize; ++i, ++j) {
      int ii(pdtLoc + i) ;

      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * m_GrainFlux[j] / pdtDOS[i]);                                // Loss of the adduct to the source
      (*CollOptr)[jj][ii]  = qd_real(rMeanOmega * m_GrainFlux[j] * sqrt(adductPopFrac[i] * Keq) / pdtDOS[i]);// Reactive gain of the source
      (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                                                      // Reactive gain (symmetrization)
      DissRateCoeff       += qd_real(m_GrainFlux[j] * adductPopFrac[i] / pdtDOS[i]);
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
      CanPrtnFn = double(m_rct1->getDOS().getSpinMultiplicity() * m_rct2->getDOS().getSpinMultiplicity()) ;
    }

    return CanPrtnFn ;
  }
  double AssociationReaction::pdtsRovibronicGrnCanPrtnFn() { return m_pdt1->getDOS().rovibronicGrnCanPrtnFn();}


  // Is reaction equilibrating and therefore contributes
  // to the calculation of equilibrium fractions.
  bool AssociationReaction::isEquilibratingReaction(double &Keq, Molecule **rct, Molecule **pdt) {

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
    Qrcts *= translationalContribution(m_rct1->getStruc().getMass(), m_rct2->getStruc().getMass(), beta);

    // rovibronic partition function for product
    const double Qpdt1 = m_pdt1->getDOS().rovibronicGrnCanPrtnFn() ;

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

  //
  // Calculate grained forward and reverse k(E)s from transition state flux
  //
  void AssociationReaction::calcGrainRateCoeffs(){

    vector<double> rctGrainDOS;
    vector<double> rctGrainEne;
    vector<double> pdtGrainDOS;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtGrainDOS) ;
    calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

    calcEffGrnThresholds();
    const int fwdThreshold = get_EffGrnFwdThreshold();
    const int rvsThreshold = get_EffGrnRvsThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxFirstNonZeroIdx = get_fluxFirstNonZeroIdx();

    const int MaximumGrain = (getEnv().MaxGrn-fluxFirstNonZeroIdx);
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain , 0.0);
    m_GrainKbmc.clear();
    m_GrainKbmc.resize(MaximumGrain , 0.0);

    for (int i = rvsThreshold, j = fluxFirstNonZeroIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKbmc[i] = m_GrainFlux[j] / pdtGrainDOS[i];
    }
    for (int i = fwdThreshold, j = fluxFirstNonZeroIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i];
    }

    // the code that follows is for printing of the f & r k(E)s
    if (getFlags().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().kbEGrainsEnabled){
      ctest << "\nk_b(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKbmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().testRateConstantEnabled)
      testRateConstant();
  }

  void AssociationReaction::testRateConstant() {
    double k_f_grained(0.0), k_b_grained(0.0), k_f_cell(0.0), k_b_cell(0.0);
    vector<double> pdtGrainDOS, pdtGrainEne, pdtCellDOS, pdtCellEne;

    const int MaximumCell = (getEnv().MaxCell);
    const int MaximumGrain = (getEnv().MaxGrn-get_fluxFirstNonZeroIdx());

    m_pdt1->getDOS().getGrainDensityOfStates(pdtGrainDOS) ;
    m_pdt1->getDOS().getCellDensityOfStates(pdtCellDOS);
    m_pdt1->getDOS().getGrainEnergies(pdtGrainEne);
    getCellEnergies(MaximumCell, pdtCellEne);


    // dissociation k(E) calculated in grains.
    const double beta = getEnv().beta;
    for (int i = 0; i < MaximumGrain; ++i){
      k_b_grained += m_GrainKbmc[i] * exp( log(pdtGrainDOS[i]) - beta * pdtGrainEne[i] ) ;
    }

    // dissociation k(E) calculated in cells.
    const vector<double>& cellFlux = get_CellFlux();
    const int fluxCellZPE = int(get_fluxZPE());
    const int pdtZPE = int(get_relative_pdtZPE());
    const int rev_threshold = fluxCellZPE - pdtZPE;
    for (int i = 0; i < MaximumCell - rev_threshold; ++i){
      k_b_cell += cellFlux[i] * exp( -beta * pdtCellEne[i + rev_threshold] ) ;
    }

    // if the partition function calculated in grain level does not differ too much to the cell level, the cell level
    // calculation can be removed.
    const double prtfn_grained = canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, beta);
    const double prtfn_cell = canonicalPartitionFunction(pdtCellDOS, pdtCellEne, beta);
    k_b_grained /= prtfn_grained;
    k_b_cell /= prtfn_cell;

    set_rvsGrnCanonicalRate(k_b_grained);
    set_rvsCellCanonicalRate(k_b_cell);

    double Keq = calcEquilibriumConstant();
    k_f_grained = get_rvsGrnCanonicalRate() * Keq;
    k_f_cell = get_rvsCellCanonicalRate() * Keq;

    set_fwdGrnCanonicalRate(k_f_grained);
    set_fwdCellCanonicalRate(k_f_cell);

    const double temperature = 1. / (boltzmann_RCpK * beta);
    ctest << endl << "Canonical pseudo first order rate constant of association reaction "
      << getName() << " = " << get_fwdCellCanonicalRate() << " s-1 (" << temperature << " K)" << endl;
    ctest << "Canonical bimolecular rate constant of association reaction "
      << getName() << " = " << get_fwdCellCanonicalRate()/m_ERConc << " cm^3/mol/s (" << temperature << " K)" << endl;
    ctest << "Canonical first order rate constant for the reverse of reaction "
      << getName() << " = " << get_rvsCellCanonicalRate() << " s-1 (" << temperature << " K)" << endl;
  }


  //
  // Calculate the rovibrational density of states of reactants.
  //
  bool AssociationReaction::calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)
  {
    std::vector<double> rctsCellDOS;
    getRctsCellDensityOfStates(rctsCellDOS);

    std::vector<double> shiftedCellDOS;
    std::vector<double> shiftedCellEne;
    const int MaximumCell = getEnv().MaxCell;
    const int cellOffset = get_pseudoIsomer()->getDOS().get_cellOffset();
    std::vector<double> rctsCellEne;
    getCellEnergies(MaximumCell, rctsCellEne);
    shiftCells(MaximumCell, cellOffset, rctsCellDOS, rctsCellEne, shiftedCellDOS, shiftedCellEne);

    const string catName = m_rct1->getName() + " + " + m_rct2->getName();

    if (getFlags().cyclePrintCellDOS){
      ctest << endl << "Cell rovibronic density of states of " << catName << endl << "{" << endl;
      for (int i = 0; i < MaximumCell; ++i){
        formatFloat(ctest, rctsCellEne[i],  6,  15) ;
        formatFloat(ctest, rctsCellDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
      getFlags().cyclePrintCellDOS = false;
    }

    calcGrainAverages(getEnv().MaxGrn, getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, grainDOS, grainEne);

    if (getFlags().cyclePrintGrainDOS){
      ctest << endl << "Grain rovibronic density of states of " << catName << endl << "{" << endl;
      for (int i = 0; i < getEnv().MaxGrn; ++i){
        formatFloat(ctest, grainEne[i],  6,  15) ;
        formatFloat(ctest, grainDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
      getFlags().cyclePrintGrainDOS = false;
    }

    return true;
  }

  //
  // Get reactants cell density of states.
  //
  void AssociationReaction::getRctsCellDensityOfStates(vector<double> &cellDOS) {
    countDimerCellDOS(m_rct1->getDOS(), m_rct2->getDOS(), cellDOS);
  }

  const int AssociationReaction::get_rctsGrnZPE(){
    double grnZpe = (m_rct1->getDOS().get_zpe()+m_rct2->getDOS().get_zpe()-getEnv().EMin) / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void AssociationReaction::calcEffGrnThresholds(void){  // calculate the effective forward and reverse
    double threshold = get_ThresholdEnergy();  // threshold energy for an association reaction
    double RxnHeat   = getHeatOfReaction();
    if (threshold < RxnHeat && m_pMicroRateCalculator->getName() == "MesmerILT"){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }

    int fluxGrnBottom   = get_fluxGrnZPE();
    int pdtGrnZPE       = m_pdt1->getColl().get_grnZPE();
    int rctsGrnZPE      = get_rctsGrnZPE();
    int GrainedRxnHeat  = pdtGrnZPE - rctsGrnZPE;

    if(threshold<0.0){                          // if the forward threshold energy is negative
      set_EffGrnFwdThreshold(0);                // forward grained flux threshold energy = 0
      set_EffGrnRvsThreshold(-GrainedRxnHeat);  // reverse grained flux threshold energy = heat of reaction
    }
    else if(threshold>0.0 && threshold<RxnHeat){// if the reverse threshold energy is negative
      set_EffGrnFwdThreshold( GrainedRxnHeat);  // forward grained flux threshold energy = heat of reaction
      set_EffGrnRvsThreshold(0);                // reverse grained flux threshold energy = 0
    }
    else{
      // forward grained flux threshold energy = TS energy - rct energy
      set_EffGrnFwdThreshold(fluxGrnBottom - rctsGrnZPE);
      // reverse grained flux threshold energy = TS energy - pdt energy
      set_EffGrnRvsThreshold(fluxGrnBottom - pdtGrnZPE);
    }
  }

}//namespace
