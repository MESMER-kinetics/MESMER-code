//-------------------------------------------------------------------------------------------
//
// PriorDistFragmentation.cpp
//
// Author: Struan Robertson
// Date:   8th October 2019
//
// This file contains the implemetation of the fragmentation abstract base class and related 
// implementation classes. 
//
//-------------------------------------------------------------------------------------------

#include "PriorDistFragmentation.h"
#include "../gWellProperties.h"

using namespace Constants;
using namespace std;

namespace mesmer
{

  // Initialize the fragment distribution.
  void priorDist::initialize(Reaction* pReaction) {

    m_pReaction = pReaction;

    Molecule* pXsRct = m_pReaction->getExcessReactant();
    Molecule* pRct = m_pReaction->get_reactant();

    vector<double> xsDOS;
    pXsRct->getDOS().getCellDensityOfStates(xsDOS);

    pRct->getDOS().getCellDensityOfStates(m_rctDOS);

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    size_t Size = xsDOS.size();
    vector<double> Trans_DOS;
    getCellEnergies(Size, pRct->getEnv().CellSize, Trans_DOS);
    for (size_t i(0); i < Trans_DOS.size(); i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]);
    }

    FastLaplaceConvolution(xsDOS, Trans_DOS, m_upperConv);

    FastLaplaceConvolution(m_upperConv, m_rctDOS, m_lowerConv);

  };

  // Calculate dissociation distribution

  void priorDist::calculate(double Energy, std::vector<double>& dist) {

    // Get the difference in zero point energies between the well and the adduct.
    const double DeltaH = m_pReaction->getHeatOfReaction();

    // Calcualte threshold for reverse reaction.
    const double rvsThreshold = m_pReaction->get_ThresholdEnergy() - DeltaH;

    // Get the excess energy available for redistribution among bimolecular species. 
    double XsE = max((Energy - rvsThreshold), 0.0);

    const size_t excessEnergy = static_cast<size_t>(XsE);

    const size_t cellOffSet = m_pReaction->getFluxCellOffset();

    if (excessEnergy > 0) {

      // Calculate cell distribution vector. The distribution is shifted by the cell-to-grain
      // off set so as to match the shift applied to the reaction flux above.

      vector<double> mDist(m_rctDOS.size(), 0.0);
      for (size_t i(0), j(cellOffSet); i < excessEnergy && j < mDist.size(); i++, j++) {
        mDist[j] = m_rctDOS[i] * m_upperConv[excessEnergy - i] / m_lowerConv[excessEnergy];
      }

      // Average cell distribution over grains.

      // Get grain size for grain averaging.
      const size_t GrainSize = m_pReaction->get_reactant()->getEnv().GrainSize;
      double sum(0.0);
      for (size_t i(0), index(0); i < dist.size(); i++) {
        for (size_t j = 0; j < (GrainSize) && index < mDist.size(); j++, index++) {
          dist[i] += mDist[index];
        }
        sum += dist[i];
      }

      // Normalize fragment distribution.
      double rSum = 1.0 / sum;
      for (size_t i(0); i < dist.size(); i++) {
        dist[i] *= rSum;
      }

    }
    return;
  }

  void modPriorDist::initialize(Reaction* pReaction) {

    m_pReaction = pReaction;

    ReactionType reactionType = m_pReaction->getReactionType();

    Molecule *pXsSpcs(NULL), *pSpcs(NULL);
    if (reactionType == PSEUDOISOMERIZATION) {
      pXsSpcs = m_pReaction->getExcessReactant();
      pSpcs = m_pReaction->get_reactant();
    }
    else if (reactionType == BIMOLECULAR_EXCHANGE) {
      vector<Molecule*> products;
      m_pReaction->get_products(products);
      pSpcs = products[0];
      pXsSpcs = products[1];
    }

    vector<double> xsDOS;
    pXsSpcs->getDOS().getCellDensityOfStates(xsDOS);

    pSpcs->getDOS().getCellDensityOfStates(m_rctDOS);

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    size_t Size = xsDOS.size();
    vector<double> Trans_DOS;
    getCellEnergies(Size, pSpcs->getEnv().CellSize, Trans_DOS);
    for (size_t i(0); i < Trans_DOS.size(); i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]);
    }

    const double beta = m_pReaction->getEnv().beta;
    double order = m_order * pow((1.0 / (beta * m_Tref * boltzmann_RCpK)), m_nexp);
    for (size_t i(0); i < m_rctDOS.size(); i++) {
      m_rctDOS[i] = pow(m_rctDOS[i], order);
    }

    FastLaplaceConvolution(xsDOS, Trans_DOS, m_upperConv);

    FastLaplaceConvolution(m_upperConv, m_rctDOS, m_lowerConv);

  };

}