//-------------------------------------------------------------------------------------------
//
// PseudoIsomerizationReaction.cpp
//
// Author: Struan Robertson
// Date:   26/May/2013
//
// This file contains the implementation of the PseudoIsomerizationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "PseudoIsomerizationReaction.h"
#include "gWellProperties.h"
#include <math.h>

using namespace Constants ;
using namespace std;

namespace mesmer
{

  void PseudoIsomerizationReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega)
  {
    // Get densities of states for detailed balance.
    vector<double> pdtDOS;
    vector<double> pdtEne;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;
    m_pdt1->getDOS().getGrainEnergies(pdtEne) ;

    // Locate isomers in system matrix.
    const size_t rctLocation = isomermap[m_rct1] ;
    const size_t pdtLocation = isomermap[m_pdt1] ;

    // Need to know the number of grouped grains in both wells.
    const size_t rShiftedGrains(m_rct1->getColl().reservoirShift());
    const size_t pShiftedGrains(m_pdt1->getColl().reservoirShift());

    // Get equilibrium constant.
    const qd_real Keq = qd_real(calcEquilibriumConstant()) ;

    const size_t pColloptrsize  = m_pdt1->getColl().get_colloptrsize() + pShiftedGrains ;
    const size_t rColloptrsize  = m_rct1->getColl().get_colloptrsize() + rShiftedGrains ;
    const size_t reverseThreshE = get_EffGrnRvsThreshold();
    const size_t fluxStartIdx   = get_fluxFirstNonZeroIdx();

    // Get Boltzmann distribution for detailed balance.
    vector<double> pdtEq ; // Population fraction of the adduct
    m_pdt1->getColl().normalizedGrnBoltzmannDistribution(pdtEq) ;
    vector<double> rctEq ; // Population fraction of the psudoisomer
    m_rct1->getColl().normalizedGrnBoltzmannDistribution(rctEq) ;

    // Note: reverseThreshE will always be greater than pNGG here

    m_fragDist->initialize(this) ;
    for ( size_t i = reverseThreshE, j = fluxStartIdx; i < pColloptrsize; ++i, ++j) {
      size_t ii(pdtLocation + i - pShiftedGrains) ;
      qd_real disMicroRateCoeff = qd_real(rMeanOmega) * qd_real(m_GrainFlux[j]) / qd_real(pdtDOS[i]) ;
      (*CollOptr)[ii][ii] -= disMicroRateCoeff ;   // Loss from adduct to pseudoisomer.

      vector<double> fragDist ;
      m_fragDist->calculate(pdtEne[i], fragDist, rColloptrsize) ;

      // Distribute adduct loss rate according to fragmentation loss rate.
      // Use detailed balance to determine gain from association reaction.

      for (size_t k = 0; k < rColloptrsize; ++k) {
        size_t jj(rctLocation + k - rShiftedGrains) ;
        qd_real eqmRatio     = Keq*qd_real(pdtEq[i])/qd_real(rctEq[k]) ;
        (*CollOptr)[jj][ii]  = disMicroRateCoeff*qd_real(fragDist[k]) ; // Gain of pseudoisomer from adduct.
        (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii]*eqmRatio ;           // Gain of adduct from pseudoisomer.
        (*CollOptr)[jj][jj] -= (*CollOptr)[ii][jj] ;                    // Loss from pseudoisomer to adduct.

        // Symmetrize system.

        (*CollOptr)[jj][ii] *= sqrt(eqmRatio) ;  
        (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ; 
      }
    }

    // Return any resources used.
    m_fragDist->clear() ;

  }

  // Initialize the fragment distribution.
  void priorDist::initialize(PseudoIsomerizationReaction *pReaction) {

    m_pReaction = pReaction ; 

    Molecule* pXsRct = m_pReaction->get_excessReactant() ;
    Molecule* pRct   = m_pReaction->get_reactant() ;

    vector<double> xsDOS;
    pXsRct->getDOS().getCellDensityOfStates(xsDOS);

    pRct->getDOS().getCellDensityOfStates(m_rctDOS);

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    size_t Size = xsDOS.size();
    vector<double> Trans_DOS ;
    getCellEnergies(Size, pRct->getEnv().CellSize, Trans_DOS) ;
    for (size_t i(0) ; i < Trans_DOS.size() ; i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]) ;
    }

    FastLaplaceConvolution(xsDOS, Trans_DOS, m_upperConv);

    FastLaplaceConvolution(m_upperConv, m_rctDOS, m_lowerConv);

  } ;

  // Calculate dissociation distribution

  void priorDist::calculate(double Energy, std::vector<double>& dist, size_t size) {

    // Get the difference in zero point energies between the well and the adduct.
    const double DeltaH  = m_pReaction->getHeatOfReaction();

    // Calcualte threshold for reverse reaction.
    const double rvsThreshold = m_pReaction->get_ThresholdEnergy() - DeltaH ;

    // Get the excess energy available for redistribution among bimolecular species. 
    double XsE = max((Energy - rvsThreshold), 0.0) ;

    const size_t excessEnergy = static_cast<size_t>(XsE);

	const size_t cellOffSet = m_pReaction->getFluxCellOffset() ;

    dist.clear() ;
    dist.resize(size,0.0) ;
    if (excessEnergy > 0) {

      // Calculate cell distribution vector. The distribution is shifted by the cell-to-grain
	  // off set so as to match the shift applied to the reaction flux above.

      vector<double> mDist(m_rctDOS.size(), 0.0) ;
      for (size_t i(0), j(cellOffSet) ; i < excessEnergy && j < mDist.size() ; i++, j++ ) {
        mDist[j] = m_rctDOS[i]*m_upperConv[excessEnergy - i]/m_lowerConv[excessEnergy] ;
      }

      // Average cell distribution over grains.

      // Get grain size for grain averaging.
      const size_t GrainSize = m_pReaction->get_reactant()->getEnv().GrainSize ;
      double sum (0.0) ;
      for (size_t i(0), index(0) ; i < dist.size() ; i++) {       
        for (size_t j=0; j<(GrainSize) && index < mDist.size(); j++, index++ ) {
          dist[i] += mDist[index];
        }
        sum += dist[i] ;
      }

	  // Normalize fragment distribution.
	  double rSum = 1.0/sum ;
      for (size_t i(0) ; i < dist.size() ; i++) {       
        dist[i] *= rSum;
	  }

    }
    return ;
  }

  void modPriorDist::initialize(PseudoIsomerizationReaction *pReaction) {

    m_pReaction = pReaction ; 

    Molecule* pXsRct = m_pReaction->get_excessReactant() ;
    Molecule* pRct   = m_pReaction->get_reactant() ;

    vector<double> xsDOS;
    pXsRct->getDOS().getCellDensityOfStates(xsDOS);

    pRct->getDOS().getCellDensityOfStates(m_rctDOS);

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    size_t Size = xsDOS.size();
    vector<double> Trans_DOS ;
    getCellEnergies(Size, pRct->getEnv().CellSize, Trans_DOS) ;
    for (size_t i(0) ; i < Trans_DOS.size() ; i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]) ;
    }

    for (size_t i(0) ; i < m_rctDOS.size() ; i++) {
      m_rctDOS[i] = pow(m_rctDOS[i], 0.25) ;
    }

    FastLaplaceConvolution(xsDOS, Trans_DOS, m_upperConv);

    FastLaplaceConvolution(m_upperConv, m_rctDOS, m_lowerConv);

  } ;

}//namespace
