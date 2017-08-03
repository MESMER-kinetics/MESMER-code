//-------------------------------------------------------------------------------------------
//
// DiffusiveLossReaction.cpp
//
// Author: Struan Robertson
// Date:   2/May/2016
//
// This file contains the implementation of the DiffusiveLossReaction class. Note reactant 
// and product are the same. 
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "DiffusiveLossReaction.h"
#include "gWellProperties.h"
#include "gStructure.h"

using namespace Constants;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool DiffusiveLossReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if (!ppReactantList)
      ppReactantList = ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1 = ppReactantList->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if (!pMol1) {
      cerr << "Cannot find molecule definition for diffusive loss " << getName() << ".";
      return false;
    }

    // Save reactant as Molecule.
    m_rct1 = pMol1;

    // Read diffusion rate parameters.
		m_diffusionRate = ppReac->XmlReadDouble("me:diffusiveLoss", true);

    return true;
  }

  //
  // Add dissociation reaction terms to collision matrix.
  //
  void DiffusiveLossReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) {
    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS);

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[m_rct1];
    const int rShiftedGrains(m_rct1->getColl().reservoirShift());

    const int colloptrsize = m_rct1->getColl().get_colloptrsize();
    const int forwardThreshE = get_EffGrnFwdThreshold();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    m_MtxGrnKf.clear();
    m_MtxGrnKf.resize(colloptrsize, 0.0);

    for (int i = fluxStartIdx, j = forwardThreshE; j < colloptrsize + rShiftedGrains; ++i, ++j) {
      int ii(rctLocation + j - rShiftedGrains);
      qd_real rtcf = qd_real(m_GrainFlux[i]) / qd_real(rctDOS[j]);
      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega) * rtcf;  // Forward loss reaction.
      m_MtxGrnKf[j - rShiftedGrains] = to_double(rtcf);
    }
  }

  //
  // Calculate grained forward k(E)s from transition state flux
  //
  void DiffusiveLossReaction::calcGrainRateCoeffs() {
    vector<double> rctGrainDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS);

    calcEffGrnThresholds();
    const int forwardTE = get_EffGrnFwdThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    const int MaximumGrain = (getEnv().MaxGrn - fluxStartIdx);
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain, 0.0);

    // calculate forward k(E)s from flux
    for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j) {
      m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i];
    }

    // the code that follows is for printing the forward k(E)s
    if (getFlags().kfEGrainsEnabled) {
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i) {
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().grainTSsosEnabled) {
      ctest << "\nN(e) for TS of " << getName() << " (referenced to " << (this->get_reactant())->getName() << " energy):\n{\n";
      for (int i = 0; i < MaximumGrain; ++i) {
        ctest << m_GrainKfmc[i] * rctGrainDOS[i] / SpeedOfLight_in_cm << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().testRateConstantEnabled)
      HighPresRateCoeffs(NULL);
  }

  // Calculate high pressure rate coefficients at current T.
  void DiffusiveLossReaction::HighPresRateCoeffs(vector<double> *pCoeffs) {

    vector<double> rctGrainDOS, rctGrainEne;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS);
    m_rct1->getDOS().getGrainEnergies(rctGrainEne);
    const size_t MaximumGrain = (getEnv().MaxGrn - get_fluxFirstNonZeroIdx());
    const double beta = getEnv().beta;

    double kf(0.0);
    for (size_t i(0); i < MaximumGrain; ++i)
      kf += m_GrainKfmc[i] * exp(log(rctGrainDOS[i]) - beta * rctGrainEne[i]);

    const double rctprtfn = canonicalPartitionFunction(rctGrainDOS, rctGrainEne, beta);
    kf /= rctprtfn;

    if (pCoeffs) {
      pCoeffs->push_back(kf);
      const double Keq = calcEquilibriumConstant();
      if (!IsNan(Keq)) {
        pCoeffs->push_back(kf / Keq);
        pCoeffs->push_back(Keq);
      }
    }
    else {
      const double temperature = 1. / (boltzmann_RCpK * beta);
      ctest << endl << "Canonical pseudo first order forward rate constant of irreversible reaction "
        << getName() << " = " << kf << " s-1 (" << temperature << " K)" << endl;
    }
  }

}//namespace
