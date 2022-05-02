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
#include "Fragmentation.h"
#include "ParseForPlugin.h"
#include <math.h>

using namespace Constants ;
using namespace std;

namespace mesmer
{

  // Initialize reaction.
  bool PseudoIsomerizationReaction::InitializeReaction(PersistPtr ppReac) {

    // Determine fragmentation model.
    PersistPtr ppDstbn = ppReac->XmlMoveTo("me:FragmentDist");
    m_fragDist = ParseForPlugin<FragDist>(this, "me:FragmentDist", ppReac);
    if (!m_fragDist)
      return false;

    m_fragDist->ReadParameters(ppDstbn, this->getName());

    return AssociationReaction::InitializeReaction(ppReac);
  };

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
    const size_t forwardThreshE = get_EffGrnFwdThreshold();
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

      vector<double> fragDist(rColloptrsize, 0.0) ;
      m_fragDist->calculate(pdtEne[i], fragDist) ;

      // Distribute adduct loss rate according to fragmentation loss rate.
      // Use detailed balance to determine gain from association reaction.

      for (size_t k = forwardThreshE, l(0); k < rColloptrsize; ++k, ++l) {
        size_t jj(rctLocation + k - rShiftedGrains) ;
        qd_real eqmRatio     = Keq*qd_real(pdtEq[i])/qd_real(rctEq[k]) ;
        qd_real Cji          = disMicroRateCoeff*qd_real(fragDist[l]) ; // Gain of pseudoisomer from adduct.
        (*CollOptr)[jj][jj] -= Cji * eqmRatio;                          // Loss from pseudoisomer to adduct.

        // Symmetrize system.

        (*CollOptr)[jj][ii] += Cji*sqrt(eqmRatio) ;
        (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ; // Gain of adduct from pseudoisomer.
      }
    }

    // Return any resources used.
    m_fragDist->clear() ;

  }

}//namespace
