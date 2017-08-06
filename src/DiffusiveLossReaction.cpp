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
		if (!ppReactant1) 
			return false;
		PersistPtr ppmol = ppReactant1->XmlMoveTo("molecule");
		if (!ppmol) 
			return false;

		const char* reftxt = ppmol->XmlReadValue("ref"); //using const char* in case NULL returned
		if (reftxt) {
			m_strRct1 = string(reftxt);
		}
		else
			return false;

		// Read diffusion rate parameters.
    m_diffusionRate = ppReac->XmlReadDouble("me:diffusiveLoss", true);

    return true;
  }

  //
  // Add diffusive loss terms to collision matrix.
  //
  void DiffusiveLossReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) {

		// Initialize diffusing species here after all the participating have been specified by
		// regular reaction terms.
		if (!m_rct1) {
			m_rct1 = m_pMoleculeManager->find(m_strRct1);
			if (!m_rct1) {
				stringstream msg;
				msg << "Failed to located diffusing species referred to in diffusive loss reaction" << getName() << ".";
				throw(std::runtime_error(msg.str()));
			}
		}
		
		// Search for diffusing species in isomers.
    Reaction::molMapType::iterator rctitr = isomermap.find(m_rct1);
    if (rctitr != isomermap.end()) {
      const int rctLocation = rctitr->second;
      const int colloptrsize = m_rct1->getColl().get_colloptrsize();

      m_MtxGrnKf.clear();
      m_MtxGrnKf.resize(colloptrsize, 0.0);

      for (int j = 0; j < colloptrsize; ++j) {
        int ii(rctLocation + j);
        (*CollOptr)[ii][ii] -= qd_real(m_diffusionRate*rMeanOmega);  // Forward loss reaction.
        m_MtxGrnKf[j] = to_double(m_diffusionRate);
      }
      return;
    }

    // Search for diffusing species in pseudoisomers.
    if (m_sourceMap) {
      Reaction::molMapType::iterator rctitr = m_sourceMap->find(m_rct1);
      const int rctLocation = rctitr->second;
      (*CollOptr)[rctLocation][rctLocation] -= qd_real(m_diffusionRate*rMeanOmega);  // Diffusive loss reaction.
      m_MtxGrnKf.clear();
      m_MtxGrnKf.push_back(m_diffusionRate);

      return;
    }

    // If MESMER gets this far the loss is in a sink molecule which is not currently catered for.
    stringstream msg;
    msg << "Diffusive loss reaction" << getName() << " refers to a species that is not represented in the collision matrix.";
    throw(std::runtime_error(msg.str()));
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

    if (pCoeffs) {
      pCoeffs->push_back(m_diffusionRate);
      pCoeffs->push_back(0.0);
      pCoeffs->push_back(0.0);
    }
    else {
      const double beta = getEnv().beta;
      const double temperature = 1. / (boltzmann_RCpK * beta);
      ctest << endl << "Diffusive loss constant for reaction "
        << getName() << " = " << m_diffusionRate << " m2s-1 (" << temperature << " K)" << endl;
    }
  }

}//namespace
