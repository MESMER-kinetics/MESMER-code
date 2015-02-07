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
#include "gWellProperties.h"
#include "gStructure.h"
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
	m_rct1 = GetMolRef(ppReactant1);
	PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
	m_rct2 = GetMolRef(ppReactant2);

	// If deficientReactantLocation=true, then swap the reactant roles.

	if (m_deficientReactantLocation) {
	  Molecule *tmp = m_rct1 ;
	  m_rct1 = m_rct2;
	  m_rct2 = tmp ;
	}

	if(!m_rct1){
	  cerr << "The deficient reactant in the association reaction " << getName() << " is undefined." << endl;
	  return false;
	}
	if(!m_rct2){
	  cerr << "The excess reactant in the association reaction " << getName() << " is undefined." << endl;
	  return false;
	}

  //Save reactant ZPEs before psudoisomer use. Restore in Finish().
  m_SavedZPE1 = m_rct1->getDOS().get_zpe();
  m_SavedZPE2 = m_rct2->getDOS().get_zpe();

	//Read product details.
	PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
	if(!ppProductList)
	  ppProductList=ppReac; //Be forgiving; we can get by without a productList element

	PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
	m_pdt1 = GetMolRef(ppProduct1);
	if (!m_pdt1) {
	  cerr << "Cannot find product molecule definition for association reaction " << getName() << ".";
	  return false;
	}

	// Read heat of reaction and rate parameters.
	return ReadRateCoeffParameters(ppReac) ;

  }

  // Reset zero point energy locations of the reactants such that
  // location of the pair is entirely on the pseudoisomer.
  double AssociationReaction::resetZPEofReactants() {

    const double pseudoIsomerZPE = get_pseudoIsomer()->getDOS().get_zpe();
    const double excessReactantZPE = get_excessReactant()->getDOS().get_zpe();
    double sourceTermZPE = pseudoIsomerZPE + excessReactantZPE;
    get_pseudoIsomer()->getDOS().set_zpe(sourceTermZPE);
    get_excessReactant()->getDOS().set_zpe(0.0);

    return sourceTermZPE;
  }

  void AssociationReaction::AddReactionTerms(qdMatrix      *CollOptr,
	molMapType    &isomermap,
	const double rMeanOmega)
  {
	// Get densities of states of the adduct for detailed balance.
	vector<double> pdtDOS;
	m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;

	// Locate isomers in system matrix.
	const int pdtLoc = isomermap[m_pdt1] ;
	const int jj     = (*m_sourceMap)[get_pseudoIsomer()] ;

	// Get equilibrium constant.
	const qd_real Keq = qd_real(calcEquilibriumConstant()) ;

	// Get Boltzmann distribution for detailed balance.
	vector<double> adductPopFrac ; // Population fraction of the adduct
	const int pShiftedGrains(m_pdt1->getColl().reservoirShift());
	m_pdt1->getColl().normalizedGrnBoltzmannDistribution(adductPopFrac) ;

	qd_real DissRateCoeff(0.0), qdMeanOmega(rMeanOmega) ;

	const int pdtRxnOptPos(pdtLoc - pShiftedGrains);
	const int colloptrsize = m_pdt1->getColl().get_colloptrsize() + pShiftedGrains ;
	const int reverseThreshE = get_EffGrnRvsThreshold();
	const int fluxStartIdx = get_fluxFirstNonZeroIdx();

	// Note: reverseThreshE will always be greater than pShiftedGrains here.

	for ( int i = reverseThreshE, j = fluxStartIdx; i < colloptrsize; ++i, ++j) {
	  int ii(pdtRxnOptPos + i) ;
	  int kk (i - pShiftedGrains);
	  qd_real Flux(m_GrainFlux[j]), dos(pdtDOS[i]), addPop(adductPopFrac[kk]) ;
	  (*CollOptr)[ii][ii] -= qdMeanOmega * Flux / dos ;                  // Loss of the adduct to the source
	  (*CollOptr)[jj][ii]  = qdMeanOmega * Flux * sqrt(Keq*addPop)/dos ; // Reactive gain of the source
	  (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                       // Reactive gain (symmetrization)
	  DissRateCoeff       += Flux * addPop / dos;
	}
	(*CollOptr)[jj][jj] -= qdMeanOmega * DissRateCoeff * Keq ;           // Loss of the source from detailed balance.
  }

  //
  // Add isomer reaction terms to contracted basis reaction matrix.
  //
  void AssociationReaction::AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap)
  {
	// Get densities of states for detailed balance.
	vector<double> pdtDOS;
	m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;

	// Get equilibrium constant.
	const double Keq = calcEquilibriumConstant() ;

	// Get Boltzmann distribution for detailed balance.
	vector<double> adductPopFrac ; // Population fraction of the adduct
	m_pdt1->getColl().normalizedGrnBoltzmannDistribution(adductPopFrac) ;

	const int pdtColloptrsize = m_pdt1->getColl().get_colloptrsize();
	const int reverseThreshE  = get_EffGrnRvsThreshold();
	const int fluxStartIdx    = get_fluxFirstNonZeroIdx();

	double DissRateCoeff(0.0) ;

	vector<double> RvsMicroRateCoef(pdtColloptrsize, 0.0) ;
	vector<double> CrsMicroRateCoef(pdtColloptrsize, 0.0) ;
	for ( int i=fluxStartIdx, j = reverseThreshE, k=0; j < pdtColloptrsize; ++i, ++j, ++k) {
	  int mm = k + reverseThreshE;
	  RvsMicroRateCoef[mm] = m_GrainFlux[i] / pdtDOS[mm] ;                          // Backward loss reaction. 
	  CrsMicroRateCoef[mm] = RvsMicroRateCoef[mm] * sqrt(adductPopFrac[mm] * Keq) ; // Reactive gain from detailed balance.
	  DissRateCoeff       += RvsMicroRateCoef[mm] * adductPopFrac[mm] ;
	}

	// Calculate the elements of the product block.

	const int pdtLocation  = isomermap[m_pdt1] ;
	const int pdtBasisSize = static_cast<int>(m_pdt1->getColl().get_nbasis());
	for (int i=0, ii(pdtLocation), egvI(pdtColloptrsize-1) ; i < pdtBasisSize ; i++, ii++, --egvI) {
	  (*CollOptr)[ii][ii] -= m_pdt1->getColl().matrixElement(egvI, egvI, RvsMicroRateCoef) ;
	  for (int j=i+1, jj(pdtLocation + j), egvJ(pdtColloptrsize-j-1)  ; j < pdtBasisSize ; j++, jj++, --egvJ) {
		qd_real tmp = m_pdt1->getColl().matrixElement(egvI, egvJ, RvsMicroRateCoef) ;
		(*CollOptr)[ii][jj] -= tmp ;
		(*CollOptr)[jj][ii] -= tmp ;
	  }
	}

	// Calculate the elements of the reactant block.

	const int jj         = (*m_sourceMap)[get_pseudoIsomer()] ;
	(*CollOptr)[jj][jj] -= qd_real(DissRateCoeff * Keq);       // Loss of the source from detailed balance.

	// Calculate the elements of the cross blocks.

	vector<double> pdtBasisVector(pdtColloptrsize, 0.0) ;
	for (int i=0, pdtEgv(pdtColloptrsize-1)  ; i < pdtBasisSize ; i++, --pdtEgv) {
	  int ii(pdtLocation + i) ;
	  qd_real tmp(0.0) ;

	  if (i==0) {

		// Special case for equilibrium eigenvectors which obey a detailed balance relation.
		// SHR, 8/Mar/2009: are there other relations like this I wonder.

		qd_real elmti = (*CollOptr)[ii][ii] ;
		qd_real elmtj = (*CollOptr)[jj][jj] ;
		tmp = sqrt(elmti*elmtj) ;

	  } else {

		// General case.

		m_pdt1->getColl().eigenVector(pdtEgv, pdtBasisVector) ;
		double sum = 0.0 ;
		for (int k(pdtColloptrsize-reverseThreshE), n(pdtColloptrsize-1); k >=0 ; --n, --k) {
		  sum += pdtBasisVector[n]*CrsMicroRateCoef[n];
		}

		tmp = qd_real(sum) ;
	  }
	  (*CollOptr)[ii][jj] += tmp ;
	  (*CollOptr)[jj][ii] += tmp ;
	}

  }

  //
  // Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
  //
  double AssociationReaction::rctsRovibronicGrnCanPrtnFn() {
	vector<double> rctGrainDOS;
	vector<double> rctGrainEne;
	calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

	// Calculate the rovibronic partition function based on the grain DOS
	return canonicalPartitionFunction(rctGrainDOS, rctGrainEne, getEnv().beta) ;
  }

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
	const double HeatOfReaction = getHeatOfReaction();
	Keq *= exp(-beta * HeatOfReaction) ;

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
	if (getFlags().grainTSsosEnabled){
	  ctest << "\nN(e) for TS of " << getName() << " (referenced to " << (this->get_pseudoIsomer())->getName() << " energy):\n{\n";
	  for (int i = 0; i < MaximumGrain; ++i){
		ctest << m_GrainKfmc[i]*rctGrainDOS[i]/SpeedOfLight_in_cm << endl;
	  }
	  ctest << "}\n";
	}
	if (getFlags().testRateConstantEnabled)
	  HighPresRateCoeffs(NULL);
  }

  // Calculate high pressure rate coefficients at current T.
  void AssociationReaction::HighPresRateCoeffs(vector<double> *pCoeffs) {

	const int MaximumGrain = (getEnv().MaxGrn-get_fluxFirstNonZeroIdx());
	vector<double> pdtGrainDOS, pdtGrainEne;
	m_pdt1->getDOS().getGrainDensityOfStates(pdtGrainDOS) ;
	m_pdt1->getDOS().getGrainEnergies(pdtGrainEne);

	// dissociation k(E) calculated in grains.
	double k_f_grained(0.0), k_b_grained(0.0);
	const double beta = getEnv().beta;
	for (int i = 0; i < MaximumGrain; ++i){
	  k_b_grained += m_GrainKbmc[i] * exp( log(pdtGrainDOS[i]) - beta * pdtGrainEne[i] ) ;
	}

	const double prtfn_grained = canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, beta);
	k_b_grained /= prtfn_grained;

	double Keq = calcEquilibriumConstant();
	k_f_grained = k_b_grained * Keq;

	if (pCoeffs) {
	  pCoeffs->push_back(k_f_grained) ;
	  pCoeffs->push_back(k_b_grained) ;
	  pCoeffs->push_back(Keq) ;
	} else {
	  const double temperature = 1. / (boltzmann_RCpK * beta);
	  ctest << endl << "Canonical pseudo first order rate constant of association reaction "
		<< getName() << " = " << k_f_grained << " s-1 (" << temperature << " K)" << endl;
	  ctest << "Canonical bimolecular rate constant of association reaction "
		<< getName() << " = " << k_f_grained/m_ERConc << " cm^3/mol/s (" << temperature << " K)" << endl;
	  ctest << "Canonical first order rate constant for the reverse of reaction "
		<< getName() << " = " << k_b_grained << " s-1 (" << temperature << " K)" << endl;
	}
  }


  //
  // Calculate the rovibrational density of states of reactants.
  //
  bool AssociationReaction::calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)
  {
	std::vector<double> rctsCellDOS;
	getRctsCellDensityOfStates(rctsCellDOS);

	const int MaximumCell = getEnv().MaxCell;

	// Get the cell offset for the source term.
	const size_t cellOffset = get_cellOffset() ;

	std::vector<double> rctsCellEne;
	getCellEnergies(MaximumCell, getEnv().CellSize, rctsCellEne);

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

	calcGrainAverages(getEnv().MaxGrn, getEnv().cellPerGrain(), cellOffset, rctsCellDOS, rctsCellEne, grainDOS, grainEne) ;

	if (getFlags().cyclePrintGrainDOS){
	  ctest << endl << "Grain rovibronic density of states of " << catName << endl << "{" << endl;
	  for (size_t i(0); i < getEnv().MaxGrn; ++i){
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
	double threshold = get_ThresholdEnergy();            // threshold energy for an association reaction
	double RxnHeat   = getHeatOfReaction(); 

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
