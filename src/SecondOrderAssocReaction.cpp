//-------------------------------------------------------------------------------------------
//
// SecondOrderAssocReaction.cpp
//
// Author: Struan Robertson
// Date:   24/Aug/2014
//
// This file contains the implementation of the SecondOrderAssocReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "SecondOrderAssocReaction.h"
#include <math.h>

using namespace Constants ;
using namespace std;

namespace mesmer
{
  //
  // Read molecular data for second order association reaction from input stream.
  // Note: the convention adopted here is that there are two (one repeated) reactants
  // and one product (adduct).
  // 
  bool SecondOrderAssocReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if(!ppReactantList)
      ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
    m_rct1 = GetMolRef(ppReactant1);
    m_rct2 = m_rct1;

    if(!m_rct1){
      cerr << "The deficient reactant in the association reaction " << getName() << " is undefined." << endl;
      return false;
    }

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

  void SecondOrderAssocReaction::AddReactionTerms(qdMatrix      *CollOptr,
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
  void SecondOrderAssocReaction::AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap)
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

  // Is reaction equilibrating and therefore contributes
  // to the calculation of equilibrium fractions.
  bool SecondOrderAssocReaction::isEquilibratingReaction(double &Keq, Molecule **rct, Molecule **pdt) {

    Keq = calcEquilibriumConstant() ;

    *rct = m_rct1 ;
    *pdt = m_pdt1 ;

    return true ;
  }

  //
  // Calculate the rovibrational density of states of reactants.
  //
  bool SecondOrderAssocReaction::calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)
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
  void SecondOrderAssocReaction::getRctsCellDensityOfStates(vector<double> &cellDOS) {
    countDimerCellDOS(m_rct1->getDOS(), m_rct2->getDOS(), cellDOS);
  }

  const int SecondOrderAssocReaction::get_rctsGrnZPE(){
    double grnZpe = (m_rct1->getDOS().get_zpe()+m_rct2->getDOS().get_zpe()-getEnv().EMin) / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void SecondOrderAssocReaction::calcEffGrnThresholds(void){  // calculate the effective forward and reverse
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
