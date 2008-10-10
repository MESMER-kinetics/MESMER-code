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

    // Read heat of reaction and rate parameters.
    return ReadRateCoeffParameters(ppReac) ;

  }

  //
  // Add (REVERSIBLE) association reaction terms to collision matrix.
  //
  //void AssociationReaction::AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega)
  //{
  //  // Locate isomers in system matrix.
  //  const int pdtLoc =      isomermap[m_pdt1] ;
  //  const int jj     = (*m_sourceMap)[get_pseudoIsomer()] ;

  //  // Get equilibrium constant.
  //  const double Keq = calcEquilibriumConstant() ;
  //  int gsz = getEnv().GrainSize;
  //  int MaximumCell = getEnv().MaxCell;

  //  const double pdtZPE = get_relative_pdtZPE();
  //  const double fluxZPE = get_fluxZPE();
  //  const int rvsThreshold = int(fluxZPE) - int(pdtZPE);

  //  // Get densities of states of the adduct in the region where flux is non-zero (pdtDOS) for detailed balance.
  //  vector<double> pdtDOS;
  //  m_pdt1->getCellDensityOfStates(pdtDOS) ;

  //  // Get Boltzmann distribution for detailed balance.
  //  vector<double> adductPopFrac;
  //  m_pdt1->normalizedCellBoltzmannDistribution(adductPopFrac);

  //  // ZPE of m_pdt1
  //  const int fluxCellOffset = int(fmod(fluxZPE, gsz));
  //  const int grnRvsThreshold  = get_EffGrnRvsThreshold();
  //  const int colloptrsize    = m_pdt1->get_colloptrsize();

  //  qd_real DissRateCoeff(0.0) ;
  //  double disRate(0.0), rg(0.0), adTs(0.0);
  //  for ( int i = grnRvsThreshold , k = 0; i < colloptrsize; ++i, ++k) {
  //    int ii(pdtLoc + i) ;
  //    int pq(k*gsz-fluxCellOffset);
  //    double sumFlux(0.0), sumDOS(0.0), grainDissFlux(0.0), sumPopFrac(0.0);
  //    for (int cx = 0; cx < gsz; ++cx){
  //      int dx(cx+pq);
  //      if (dx >= 0){
  //        sumFlux       += m_CellFlux[dx];
  //        sumPopFrac    += adductPopFrac[dx + rvsThreshold];
  //        DissRateCoeff += m_CellFlux[dx] * adductPopFrac[dx + rvsThreshold] / pdtDOS[dx + rvsThreshold];
  //      }
  //      sumDOS += pdtDOS[dx + rvsThreshold];
  //    }
  //    (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * sumFlux / sumDOS);                  // Loss of the adduct to the source

  //    (*CollOptr)[jj][ii] = qd_real(rMeanOmega * sumFlux * sqrt(sumPopFrac * Keq) / sumDOS);  // Reactive gain of the source

  //    (*CollOptr)[ii][jj] = (*CollOptr)[jj][ii];                                      // Reactive gain (symmetrization)

  //    adTs = to_double((*CollOptr)[ii][ii]);
  //    disRate = to_double(DissRateCoeff);
  //    rg = to_double((*CollOptr)[jj][ii]);
  //  }
  //  const double assR = to_double(DissRateCoeff * Keq);
  //  (*CollOptr)[jj][jj]   -= qd_real(rMeanOmega * DissRateCoeff * Keq);       // Loss of the source from detailed balance.
  //}
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

    m_pdt1->normalizedGrnBoltzmannDistribution(adductPopFrac, MaximumGrain) ;
    qd_real DissRateCoeff(0.0) ;

    const int colloptrsize = m_pdt1->get_colloptrsize();
    const int forwardThreshE = get_EffGrnFwdThreshold();
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
      CanPrtnFn = double(m_rct1->getSpinMultiplicity() * m_rct2->getSpinMultiplicity()) ;
    }

    return CanPrtnFn ;
  }
  double AssociationReaction::pdtsRovibronicGrnCanPrtnFn() { return m_pdt1->rovibronicGrnCanPrtnFn();}


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

  //
  // Calculate grained forward and reverse k(E)s from transition state flux
  //
  void AssociationReaction::calcGrainRateCoeffs(){

    vector<double> rctGrainDOS;
    vector<double> rctGrainEne;
    vector<double> pdtGrainDOS;
    m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;
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

    m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;
    m_pdt1->getCellDensityOfStates(pdtCellDOS);
    m_pdt1->getGrainEnergies(pdtGrainEne);
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
    const int cellOffset = get_pseudoIsomer()->get_cellOffset();
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

    calcGrainAverages(getEnv().MaxGrn, getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, grainDOS, grainEne, catName);

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
    get_rctsDensityOfStatesCalculator()->countDimerCellDOS(m_rct1, m_rct2, cellDOS);
  }

  const int AssociationReaction::get_rctsGrnZPE(){
    double grnZpe = (m_rct1->get_zpe()+m_rct2->get_zpe()-getEnv().EMin) / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void AssociationReaction::calcEffGrnThresholds(void){  // calculate the effective forward and reverse
    double threshold = get_ThresholdEnergy();  // threshold energy for an association reaction
    double RxnHeat   = getHeatOfReaction();
    if (threshold < RxnHeat && m_pMicroRateCalculator->getName() == "Mesmer ILT"){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }

    int fluxGrnBottom   = get_fluxGrnZPE();
    int pdtGrnZPE       = m_pdt1->get_grnZPE();
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
