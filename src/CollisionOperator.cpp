//-------------------------------------------------------------------------------------------
//
// CollisionOperator.cpp
//
// Author: Struan Robertson
// Date:   26/Feb/2011
//
// This file contains implementation of the master equation collision operator class.
//
//-------------------------------------------------------------------------------------------
#include <numeric>
#include <set>
#include "CollisionOperator.h"
#include "ConditionsManager.h"

#include "AnalysisData.h"
#include "AssociationReaction.h"
#include "IrreversibleUnimolecularReaction.h"
#include "IsomerizationReaction.h"
#include "IrreversibleExchangeReaction.h"
#include "BimolecularSinkReaction.h"
#include "gWellProperties.h"
#include "gPopulation.h"

namespace mesmer
{

	CollisionOperator::CollisionOperator() : m_pMoleculeManager(0),
		m_pReactionManager(0),
		m_isomers(),
		m_sources(),
		m_sinkRxns(),
		m_meanOmega(0.0),
		m_precision(DOUBLE),
		m_reactionOperator(0),
		m_eigenvectors(0),
		m_eigenvalues(),
		m_SpeciesSequence(),
		m_eqVector(),
		m_eqVecSize(0),
		m_punchSymbolGathered(false),
		m_GrainProfileAtTimeData(),
		m_phenomenlogicalRates() {}

	CollisionOperator::~CollisionOperator() {
		if (m_reactionOperator) delete m_reactionOperator;
	}

	// Initialize the collision operator object.
	bool CollisionOperator::initialize(MoleculeManager *pMoleculeManager, ReactionManager *pReactionManager, ConditionsManager *pConditionManager) {

		if ((m_pMoleculeManager = pMoleculeManager) &&
			(m_pReactionManager = pReactionManager) &&
			(m_pConditionManager = pConditionManager)) return true;

		return false;
	};

	//
	// Main methods for constructing the Collision operator.
	//
	bool CollisionOperator::BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags, bool writeReport)
	{
		const double SUPREMUM = 9e23;
		const double INFIMUM = -SUPREMUM;
		//
		// Find all the unique wells and lowest zero point energy.
		//
		m_isomers.clear();
		m_sources.clear(); // Maps the location of source in the system matrix.

		double minEnergy(SUPREMUM); // The minimum & maximum ZPE amongst all wells, set artificially large and small
		double maxEnergy(INFIMUM);  // to guarantee that each is overwritten in setting minEnergy and maxEnergy.

		Molecule *pBathGasMolecule = m_pMoleculeManager->find(mEnv.bathGasName);
		if (!pBathGasMolecule)
		{
			cerr << "The molecular data for the bath gas " << mEnv.bathGasName << " has not been found" << endl;
			return false;
		}

		// populate molMapType with unimolecular species and determine minimum/maximum energy on the PES
		for (size_t i(0); i < m_pReactionManager->size(); ++i) {

			Reaction *pReaction = (*m_pReactionManager)[i];

			// Reset the the microcanonical re-calculation flags if required.
			if (!mFlags.useTheSameCellNumber) pReaction->resetCalcFlag();

			// Transition State
			// third check for the transition state in this reaction
			Molecule *pTransitionState = pReaction->get_TransitionState();
			double TS_ZPE(INFIMUM);
			if (pTransitionState) {
				TS_ZPE = pTransitionState->getDOS().get_zpe();
				maxEnergy = max(maxEnergy, TS_ZPE);
			}

			// unimolecular species
			vector<Molecule *> unimolecules;
			pReaction->get_unimolecularspecies(unimolecules);
			// populate molMapType with unimolecular species
			for (size_t j(0); j < unimolecules.size(); ++j) {
				// wells
				Molecule *pCollidingMolecule = unimolecules[j];
				const double collidingMolZPE(pCollidingMolecule->getDOS().get_zpe());
				if (pCollidingMolecule && m_isomers.find(pCollidingMolecule) == m_isomers.end()) { // New isomer
					m_isomers[pCollidingMolecule] = 0; //initialize to a trivial location

					minEnergy = min(minEnergy, collidingMolZPE);
					maxEnergy = max(maxEnergy, collidingMolZPE);
				}

				//calculate the lowest barrier associated with this well(species)
				if (TS_ZPE != INFIMUM) {
					const double barrierHeight = TS_ZPE - collidingMolZPE;
					if (barrierHeight < pCollidingMolecule->getColl().getLowestBarrier()) {
						pCollidingMolecule->getColl().setLowestBarrier(barrierHeight);
					}
				}
			}

			//
			// For Association reactions determine zero point energy location of the
			// associating pair.
			//
			AssociationReaction *pAReaction = dynamic_cast<AssociationReaction*>(pReaction);
			if (pAReaction) {
				double sourceTermZPE = pAReaction->resetZPEofReactants();
				minEnergy = min(minEnergy, sourceTermZPE);
				maxEnergy = max(maxEnergy, sourceTermZPE);

				// Calculate the lowest barrier associated with this well(species)
				// For association reaction, it is assumed that the barrier height is close to the source term energy
				// and in a sense, it is preferable to set this variable to the source term energy even there is an explicit
				// transition state.
				double adductZPE = unimolecules[0]->getDOS().get_zpe();
				double barrierHeight = sourceTermZPE - adductZPE;
				if (barrierHeight < unimolecules[0]->getColl().getLowestBarrier()) {
					unimolecules[0]->getColl().setLowestBarrier(barrierHeight);
				}
			}

			//
			// For irreversible exchange reactions determine zero point energy location of the
			// associating pair.
			//
			IrreversibleExchangeReaction *pIEReaction = dynamic_cast<IrreversibleExchangeReaction*>(pReaction);
			if (pIEReaction) {
				double pseudoIsomerZPE = pIEReaction->get_pseudoIsomer()->getDOS().get_zpe();
				double excessReactantZPE = pIEReaction->get_excessReactant()->getDOS().get_zpe();
				double sourceTermZPE = pseudoIsomerZPE + excessReactantZPE;
				minEnergy = min(minEnergy, sourceTermZPE);
				maxEnergy = max(maxEnergy, sourceTermZPE);

				// There is no well for this reaction
			}

			//
			// For irreversible unimolecular reactions determine zero point energy location of the barrier
			//
			IrreversibleUnimolecularReaction *pDissnRtn = dynamic_cast<IrreversibleUnimolecularReaction*>(pReaction);
			if (pDissnRtn) {
				const double rctZPE = pDissnRtn->get_reactant()->getDOS().get_zpe();
				double barrierZPE = rctZPE + pDissnRtn->get_ThresholdEnergy();
				minEnergy = min(minEnergy, barrierZPE);
				maxEnergy = max(maxEnergy, barrierZPE);

				// Calculate the lowest barrier associated with this well(species).
				if (barrierZPE < unimolecules[0]->getColl().getLowestBarrier()) {
					unimolecules[0]->getColl().setLowestBarrier(barrierZPE);
				}
			}

			// drg 15 Dec 2011
			// For bimolecular sink reactions, determine the zero point energy of the bimolecular barrier (presently zero)
			//
			BimolecularSinkReaction *pBimSinkRxn = dynamic_cast<BimolecularSinkReaction*>(pReaction);
			if (pBimSinkRxn) {
				const double rctZPE = pBimSinkRxn->get_reactant()->getDOS().get_zpe();
				double barrierZPE = rctZPE + pBimSinkRxn->get_ThresholdEnergy();
				minEnergy = min(minEnergy, barrierZPE);
				maxEnergy = max(maxEnergy, barrierZPE);
			}

			//
			// Find all source terms. Note: a source term contains the deficient reactant.
			// It is possible for there to be more than one source term.
			//
			pReaction->updateSourceMap(m_sources);

		}

		// Set grain parameters for the current Temperature/pressure condition.
		if (!SetGrainParams(mEnv, mFlags, minEnergy, maxEnergy, writeReport))
			return false;

		// Calculate flux and k(E)s
		for (size_t i(0); i < m_pReactionManager->size(); ++i) {
			if (!(*m_pReactionManager)[i]->calcGrnAvrgMicroRateCoeffs())
				return false;
		}

		if (!mFlags.rateCoefficientsOnly) {
			//
			// Shift all wells to the same origin, calculate the size of the reaction operator,
			// calculate the mean collision frequency and initialize all collision operators.
			// Set grain ZPE (with respect to the minimum of all wells)
			//
			m_meanOmega = 0.0;
			Reaction::molMapType::iterator isomeritr = m_isomers.begin();
			for (; isomeritr != m_isomers.end(); ++isomeritr) {

				Molecule *isomer = isomeritr->first;
				int colloptrsize = mEnv.MaxGrn - isomer->getColl().get_grnZPE();
				isomer->getColl().set_colloptrsize(colloptrsize);

				if (!isomer->getColl().initCollisionOperator(mEnv, pBathGasMolecule)) {
					cerr << "Failed initializing collision operator for " << isomer->getName() << endl;
					return false;
				}

				m_meanOmega += isomer->getColl().get_collisionFrequency();
			}
			m_meanOmega /= double(m_isomers.size());

			// Build reaction operator.
			//
			// One of two methods for building the reaction operator are available:
			// the conventional energy grained master equation method which is based
			// on energy grains and a contracted basis set method in which a basis
			// set is generated from the individual collision operators and a
			// representation of the reaction operator build upon this basis.

			if (!mEnv.useBasisSetMethod) {

				// Full energy grained reaction operator.
				switch (m_precision) {
				case DOUBLE_DOUBLE:
					constructGrainMatrix<dd_real>(mEnv, mFlags);
					break;
				case QUAD_DOUBLE:
					constructGrainMatrix<qd_real>(mEnv, mFlags);
					break;
				default:
					constructGrainMatrix<double>(mEnv, mFlags);
				}
			}
			else {

				// Contracted basis set reaction operator.
				constructBasisMatrix(mEnv);

			}

			// Locate all sink terms.
			locateSinks();

		}

		return true;
	}

	// Sets grain parameters and determine system environment.
	bool CollisionOperator::SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne, bool writeReport)
	{
		//  Grain size and number of grain:
		//
		//  - Either grain size or number of grains can be specified, but not both.
		//
		//  - Uses the value of grain size in the datafile, if specified.
		//
		//  - If grain size is not specified but number of grains is, use a grain size to fit the energy range.
		//  If neither is specified, the grain size is set to 100cm-1 and the number of grains set so that
		//  the energy range is sufficient.
		//
		//  Energy Range:
		//
		//  - The required total energy domain extends from the lowest zero point energy of the lowest molecule
		//  to 10 k_B T above the highest.

		static bool bcalGrainSize(false);
		static bool bcalGrainNum(false);

		//
		// Sanity check for the energy range.
		//
		const double maxERange(1.e9);
		if (fabs(maxEne - minEne) > maxERange) {
			throw(std::runtime_error("__FUNCTION__: The energy range of all species exceeds 1.e9 cm-1."));
		}

		mEnv.EMin = minEne;
		mEnv.EMax = maxEne;

		// For testing purposes, set the maxGrn based on the highest temperature we use in all calculations.
		const double MaximumTemperature = mEnv.MaximumTemperature;

		// Calculate the maximum energy cut-off based on temperature.
		const double thermalEnergy = (mFlags.useTheSameCellNumber) ? MaximumTemperature * boltzmann_RCpK : 1.0 / mEnv.beta;
		mEnv.EMax += mEnv.EAboveHill * thermalEnergy;

		// Check cut-off against population.
		if (mFlags.autoSetMaxEne) {
			SetMaximumCellEnergy(mEnv, mFlags);
		}

		//Reset max grain and grain size unless they have been specified in the conditions section.
		if (bcalGrainNum)  mEnv.MaxGrn = 0;
		if (bcalGrainSize) mEnv.GrainSize = 0;

		if (mEnv.MaxGrn > 0 && mEnv.GrainSize <= 0) {
			mEnv.GrainSize = int((mEnv.EMax - mEnv.EMin) / double(mEnv.MaxGrn)) + 1;
			bcalGrainSize = true;
		}
		else if (mEnv.GrainSize > 0 && mEnv.MaxGrn <= 0) {
			mEnv.MaxGrn = int((mEnv.EMax - mEnv.EMin) / double(mEnv.GrainSize)) + 1;
			bcalGrainNum = true;
		}
		else if (mEnv.GrainSize <= 0 && mEnv.MaxGrn <= 0) {
			mEnv.GrainSize = 100; //default 100cm-1
			cerr << "Grain size was invalid. Reset grain size to default: 100" << once << endl;
			mEnv.MaxGrn = int((mEnv.EMax - mEnv.EMin) / double(mEnv.GrainSize)) + 1;
			bcalGrainNum = true;
		}
		else if (mEnv.GrainSize > 0 && mEnv.MaxGrn > 0) {
			cerr << "Both grain size and number of grains specified. Grain size used" << once << endl;
			mEnv.MaxGrn = int((mEnv.EMax - mEnv.EMin) / double(mEnv.GrainSize)) + 1;
			bcalGrainNum = true;
		}

		mEnv.MaxCell = max(mEnv.MaxCell, (mEnv.MaxGrn * size_t(double(mEnv.GrainSize) / mEnv.CellSize)));

		if (writeReport) cinfo << "Number of cells = " << mEnv.MaxCell << ", Number of grains = " << mEnv.MaxGrn << once << endl;

		return true;
	}

	// The following method checks to see if any of the principal species has an equilibrium
	// population that is greater than that of a specified threshold. If it is, it attempts to
	// find a new upper limit of the energy cut-off based on population.
	void  CollisionOperator::SetMaximumCellEnergy(MesmerEnv &mEnv, const MesmerFlags& mFlags)
	{
		// Locate the principal species: Iterate through all isomers.
		vector<Molecule *> species;
		Reaction::molMapType::iterator isomeritr = m_isomers.begin();
		for (; isomeritr != m_isomers.end(); ++isomeritr) {
			species.push_back(isomeritr->first);
		}

		// Then through all pseudoIsomers.
		Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin();
		for (; pseudoIsomeritr != m_sources.end(); ++pseudoIsomeritr) {
			species.push_back(pseudoIsomeritr->first);
		}

		mEnv.MaxCell = int(mEnv.EMax - mEnv.EMin);
		const double populationThreshold(mFlags.popThreshold);
		double HighCell(0.0);
		for (size_t i(0); i < species.size(); ++i) {
			Molecule *pmol = species[i];
			vector<double> cellFrac;
			pmol->getColl().normalizedCellBoltzmannDistribution(cellFrac);

			// Offset cell size by relative energy of species.
			double Rel_ZPE(pmol->getDOS().get_zpe() - mEnv.EMin);
			size_t cutoffCell = mEnv.MaxCell - 1 - size_t(Rel_ZPE);
			if (cellFrac[cutoffCell] > populationThreshold) {

				// Find cell at which population threshold is reached.
				bool flag(true);
				size_t maxcell(0);
				for (size_t i(1); i < cellFrac.size() && flag; i++) {
					if (cellFrac[i] < populationThreshold && cellFrac[i] < cellFrac[i - 1]) {
						maxcell = i;
						flag = false;
					}
				}
				if (flag) {
					maxcell = cellFrac.size();
					cerr << "Warning: The equilbrum population of species " << pmol->getName()
						<< " is greater than the specefied cutt-off " << populationThreshold << endl;
				}

				HighCell = max(HighCell, double(maxcell));
			}
		}
		mEnv.EMax += HighCell;

	}

	// This method constructs a transition matrix based on energy grains.
	//
	template<class T>
	void CollisionOperator::constructGrainMatrix(MesmerEnv &mEnv, MesmerFlags& mFlags) {

		// Determine the size and location of various blocks.

		// 1. Isomers.

		size_t msize(0);
		Reaction::molMapType::iterator isomeritr = m_isomers.begin();
		for (; isomeritr != m_isomers.end(); ++isomeritr) {
			Molecule *isomer = isomeritr->first;
			isomeritr->second = static_cast<int>(msize); //set location
			msize += isomer->getColl().get_colloptrsize() - isomer->getColl().reservoirShift();
		}

		// 2. Pseudoisomers.

		Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin();
		for (; pseudoIsomeritr != m_sources.end(); ++pseudoIsomeritr) {
			pseudoIsomeritr->second = static_cast<int>(msize); //set location
			msize++;
		}

		m_eqVecSize = msize;

		// Allocate space for the full system collision operator.
		if (m_reactionOperator) delete m_reactionOperator;
		m_reactionOperator = new qdMatrix(msize, 0.0);

		// Insert collision operators to reaction operator from individual wells.
		isomeritr = m_isomers.begin();
		for (isomeritr = m_isomers.begin(); isomeritr != m_isomers.end(); ++isomeritr) {

			Molecule *isomer = isomeritr->first;
			double omega = isomer->getColl().get_collisionFrequency();
			size_t idx = isomeritr->second;

			TMatrix<T> *egme(NULL);
			if (!isomer->getColl().collisionOperator(mEnv, &egme)) {
				string errorMsg = "Failed building collision operator for " + isomer->getName() + ".";
				throw(std::runtime_error(errorMsg));
			}
			isomer->getColl().copyCollisionOperator(m_reactionOperator, egme, idx, omega / m_meanOmega);
			delete egme;

		}

		// Add connecting rate coefficients.
		for (size_t i(0); i < m_pReactionManager->size(); ++i) {
			(*m_pReactionManager)[i]->AddReactionTerms(m_reactionOperator, m_isomers, 1.0 / m_meanOmega);
		}

		// Add diffusive Loss terms if required.
		AddDiffusiveLossTerms(1.0 / m_meanOmega);

	}

	//
	// Add diffusive loss terms to collision matrix.
	//
	void CollisionOperator::AddDiffusiveLossTerms(const double rMeanOmega) {

		// Locate data for diffusion loss coefficient.

		int currentSet = m_pConditionManager->get_currentData();
		if (currentSet < 0)
			return;

		std::vector<RawDataSet>& exptDataSets = m_pConditionManager->get_experimentalrawDataSets(size_t(currentSet));

		if (exptDataSets.size() == 0) {
			// Nothing to do so return.
			return;
		}

		RawDataSet& expData = exptDataSets[0];
		string ref1 = expData.m_ref1;

		// Initialize diffusing species here, after all the participating have been specified by
		// regular reaction terms.
		Molecule* pMol = m_pMoleculeManager->find(ref1);
		if (!pMol) {
			stringstream msg;
			msg << "Failed to located diffusing species referred to in diffusive loss reaction" << pMol->getName() << ".";
			throw(std::runtime_error(msg.str()));
		}

		double diffusionRate = expData.m_pPersistPtr->XmlReadDouble("diffusiveLoss");

		// Search for diffusing species in isomers.
		Reaction::molMapType::iterator rctitr = m_isomers.find(pMol);
		if (rctitr != m_isomers.end()) {
			const int rctLocation = rctitr->second;
			const int colloptrsize = pMol->getColl().get_colloptrsize();

			for (int j = 0; j < colloptrsize; ++j) {
				int ii(rctLocation + j);
				(*m_reactionOperator)[ii][ii] -= qd_real(diffusionRate*rMeanOmega);  // Diffusive loss reaction.
			}
			return;
		}

		// Search for diffusing species in pseudoisomers.
		if (m_sources.size()) {
			Reaction::molMapType::iterator rctitr = m_sources.find(pMol);
			const int rctLocation = rctitr->second;
			(*m_reactionOperator)[rctLocation][rctLocation] -= qd_real(diffusionRate*rMeanOmega);  // Diffusive loss reaction.

			return;
		}
	}

	// This is a routine to construct the big basis matrix based on the alternative basis set method.
	// The full reaction operator is subject to a similarity transformation process by a set of eigenvectors.
	// If there are three wells and two sources in the system, and the eigenvectors of each well under the assumption
	// of the conservation of the wells are U_0, U_1 and U_2, respectively. The transformer matrix should look like
	//
	//        [  U_0   0    0   0   0 ]
	//        [   0   U_1   0   0   0 ]
	//    U = [   0    0   U_2  0   0 ]
	//        [   0    0    0   1   0 ]
	//        [   0    0    0   0   1 ]
	//
	// This transformer matrix operates on the reaction operator to give the basis matrix by doing
	//
	//     M'' = U^-1 M U
	//
	// One then needs to decide how many members of this basis matrix to include in the reduced basis matrix for
	// diagonalization.
	//
	void CollisionOperator::constructBasisMatrix(MesmerEnv &mEnv) {

		// Determine the size and location of various blocks.

		// 1. Isomers.

		size_t msize(0);
		Reaction::molMapType::iterator isomeritr = m_isomers.begin();
		for (; isomeritr != m_isomers.end(); ++isomeritr) {
			Molecule *isomer = isomeritr->first;
			isomeritr->second = static_cast<int>(msize); //set location
			msize += isomer->getColl().get_nbasis();
		}

		// 2. Pseudoisomers.

		Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin();
		for (; pseudoIsomeritr != m_sources.end(); ++pseudoIsomeritr) {
			pseudoIsomeritr->second = static_cast<int>(msize); //set location
			msize++;
			m_eqVecSize++;
		}

		// Allocate space for the reaction operator.

		if (m_reactionOperator) delete m_reactionOperator;
		m_reactionOperator = new qdMatrix(msize, 0.0);

		// Insert collision operators: in the contracted basis these are the eignvalues
		// of the isomer collision operators.
		for (isomeritr = m_isomers.begin(); isomeritr != m_isomers.end(); ++isomeritr) {

			Molecule *isomer = isomeritr->first;
			double omega = isomer->getColl().get_collisionFrequency();
			size_t idx = isomeritr->second;

			qdMatrix *egme(NULL);
			isomer->getColl().collisionOperator(mEnv, &egme);
			isomer->getColl().diagonalizeCollisionOperator(egme);
			isomer->getColl().copyCollisionOperatorEigenValues(m_reactionOperator, idx, omega);
			delete egme;
		}

		// Add rate coefficients.
		for (size_t i(0); i < m_pReactionManager->size(); ++i) {
			(*m_pReactionManager)[i]->AddContractedBasisReactionTerms(m_reactionOperator, m_isomers);
		}

		// Print out system matrix.

		//stest << endl << "System matrix:" << endl << endl ;
		//for (size_t i(0) ; i < msize ; ++i) {
		//  for (size_t j(0) ; j < msize ; ++j) {
		//    formatFloat(stest, (*m_reactionOperator)[i][j],  6,  15) ;
		//  }
		//  stest << endl ;
		//}

	}

	bool CollisionOperator::calculateEquilibriumFractions()
	{ /* Consider a three well system: e.g., A <-> B <-> C where A <-> B has Keq = K1 & B <-> C has Keq = K2.
		This routine uses the fact that the normalized equilibrated system may be described
		by a 3x3 matrix and a vector which satisfy the following:
		|-K1  1   0| |A|   |0|
		| 0  -K2  1| |B| = |0|
		| 1   1   1| |C|   |1|
		The equilibrium fraction of each isomer (or pseudo isomer, in the case of a source term) may be
		obtained by inverting the matrix shown above, and taking the elements in the final column of the inverse.
		Any system, with an arbitrary number of wells and connections, may be described by such a Matrix */

		// Determine the total number of isomers + sources from the m_isomers and m_sources
		// maps, the number of independent reactions needed to specify equilibrium and
		// intialize the matrix which holds the system of equations that describe the
		// equilibrium distribution.
		size_t nTotalNumSpecies(m_isomers.size() + m_sources.size());
		size_t nIndependRxns(nTotalNumSpecies - 1);
		qdMatrix eqMatrix(nTotalNumSpecies, qd_real(0.0));
		for (size_t i(0); i < nTotalNumSpecies; ++i) {
			eqMatrix[nTotalNumSpecies - 1][i] = qd_real(1.0);
		}

		// Initialize a map of equilibrium fractions.
		m_SpeciesSequence.clear();

		// The flag checks for a second order source and if present will invoke an
		// iterative solution. It will also cause an error to be thrown if there is 
		// more than one second order source. SHR 19/Oct/2014: I am not sure that 
		// this latter restricition is necessary but has been adopted for temporary
		// conveniance.
		bool bSecondOrderFlag(false);
		size_t idxr(0), idxs(0);

		// Counters to keep track of how many elements are in m_SpeciesSequence and reactions examined.
		size_t reactionCount(0), counter(0);

		// The followong vector keeps track of how isomers/pseudoisomers are connected. The
		// ordering of the vector is the same as m_SpeciesSequence and each set contains the
		// isomers a give isomer is connected to, either direactly or indireactly. This object
		// is used to determine if a reaction is redundant or not.
		vector<set<Molecule*> > Connections(nTotalNumSpecies);

		// Loop over the number of reactions in order to assign elements to the m_SpeciesSequence
		// map and then update the corresponding matrix elements in eqMatrix. The iteration through
		// reactions continues until enough information has been gathered to establish the unique
		// equilibrium values for all species. Note: the system may be over determined so not all
		// reactions need to be inspected to obtain the required information.
		for (size_t i(0); i < m_pReactionManager->size() && reactionCount < nIndependRxns; ++i) {

			Molecule* rct;
			Molecule* pdt;
			double Keq(0.0);

			// Only need eq fracs for species in isom & assoc rxns.
			if ((*m_pReactionManager)[i]->isEquilibratingReaction(Keq, &rct, &pdt)) {

				size_t ploc(0), rloc(0);
				bool bNewRec(false), bNewPdt(false);

				Reaction::molMapType::iterator rctitr = m_SpeciesSequence.find(rct);
				if (rctitr != m_SpeciesSequence.end()) {
					rloc = rctitr->second;
				}
				else {
					m_SpeciesSequence[rct] = counter;
					rloc = counter;
					counter++;
					bNewRec = true;
				}
				set<Molecule*> &rSet = Connections[rloc];
				rSet.insert(pdt);

				Reaction::molMapType::iterator pdtitr = m_SpeciesSequence.find(pdt);
				if (pdtitr != m_SpeciesSequence.end()) {
					ploc = pdtitr->second;
				}
				else {
					m_SpeciesSequence[pdt] = counter;
					ploc = counter;
					counter++;
					bNewPdt = true;
				}
				set<Molecule*> &pSet = Connections[ploc];
				pSet.insert(rct);

				// Determine if reaction is redundant.
				if (bNewRec || bNewPdt) {
					// Any reaction with a species not in the m_SpeciesSequence list is added.
					// Also any conections either reactant or product have should be shared.
					for (set<Molecule*>::const_iterator itr = rSet.begin(); itr != rSet.end(); itr++) {
						pSet.insert(*itr);
					}
					rSet = pSet;
				}
				else {
					// Both species are already in the list, but this does not mean reaction is redundant.
					if (pSet.find(rct) != pSet.end() || rSet.find(pdt) != rSet.end()) {
						continue;
					}
				}

				if ((*m_pReactionManager)[i]->getReactionType() == SECONDORDERASSOCIATION) {
					if (bSecondOrderFlag) {
						throw(std::runtime_error("__FUNCTION__:  Multiple second order sources defined."));
					}
					else {
						bSecondOrderFlag = true;
						idxr = reactionCount;
						idxs = rloc;
						eqMatrix[reactionCount][ploc] += qd_real(1.0);
						eqMatrix[reactionCount][rloc] -= qd_real(2.0*Keq);
						eqMatrix[nTotalNumSpecies - 1][reactionCount] = qd_real(2.0);
					}
				}
				else {
					eqMatrix[reactionCount][ploc] += qd_real(1.0);
					eqMatrix[reactionCount][rloc] -= qd_real(Keq);
				}

				reactionCount++;
			}
		}

		if (reactionCount < nIndependRxns)
			throw(std::runtime_error(string(__FUNCTION__) 
				+ string(":\n The total number of independent species is greater than the number of unique reactions - the system is under specified.")
				+ string("\n This may be because two independent sets of reactions are being defined in one system.\n")));

		// If counter==0 after the for loop above, then there are no equilibrating reactions
		// (i.e., all the reactions are irreversible).  In that case, the lone isomer has an
		// equilibrium fraction of 1. Thus, we increment counter so that the 1 is added to
		// the eqMatrix in the for loop immediately following.
		if (counter == 0) {
			if (m_isomers.size()) {
				Molecule* rct = (m_isomers.begin())->first;
				m_SpeciesSequence[rct] = counter;
			}
			else if (m_sources.size()) {
				Molecule* rct = (m_sources.begin())->first;
				m_SpeciesSequence[rct] = counter;
			}
			else {
				return false;
			}
			++counter;
		}

		qdMatrix backup(eqMatrix);  //backup EqMatrix for error reporting

		stest << endl << "Eq fraction matrix:" << endl;
		backup.showFinalBits(counter);

		vector<qd_real> eqFraction(eqMatrix.size(), 0.0);
		if (bSecondOrderFlag) {
			iterativeEquiSln(eqMatrix, eqFraction, idxr, idxs);
		}
		else {
			eqFraction.back() = qd_real(1.0);
			eqMatrix.invertLUdecomposition();
			eqFraction *= eqMatrix;
		}

		stest << "inverse of Eq fraction matrix:" << endl;
		eqMatrix.showFinalBits(counter);

		Reaction::molMapType::iterator itr1;

		// Assign Eq fraction to appropriate Molecule in the Eq. frac map.
		for (itr1 = m_SpeciesSequence.begin(); itr1 != m_SpeciesSequence.end(); ++itr1) {
			int seqMatrixLoc = itr1->second;
			Molecule* key = itr1->first;
			key->getPop().setEqFraction(to_double(eqFraction[seqMatrixLoc]));
			string speciesName = key->getName();
			stest << "Equilibrium Fraction for " << speciesName << " = " << key->getPop().getEqFraction() << endl;
		}

		// Calculate equilibrium vector.
		if (!produceEquilibriumVector()) {
			cerr << "Calculation of equilibrium vector failed.";
			return false;
		}

		return true;
	}

	// Calculate the equilibrium fraction for systems with second order terms
	// using an iterative approach - a Newton-Raphson root finding algorithm 
	// is used. 
	bool CollisionOperator::iterativeEquiSln(qdMatrix &eqMatrix, vector<qd_real> &eqFraction, size_t idxr, size_t idxs) {

		size_t niterations(100);
		vector<qd_real> delta;
		eqFraction[idxs] = qd_real(1.0);
		const qd_real tol(1.0e-05);
		qdMatrix jacobian(eqMatrix.size());
		bool converged(false);
		for (size_t i(0); i < niterations && !converged; ++i) {

			// Construct the latest estimate of the Jacobian.
			jacobian = eqMatrix;
			jacobian[idxr][idxs] *= eqFraction[idxs];

			// Construct the latest estimate of the difference vector.
			// For this the value of the vector function at the current
			// iteration point is required. The profuct of the Jacobian
			// and the current best estimate is already very close to 
			// this and manipulations below correct for factors of
			// two that happen during differenitation.
			qdMatrix fn(jacobian);
			fn[idxr][idxs] *= qd_real(0.5);
			delta = eqFraction;
			delta *= fn;
			delta.back() -= qd_real(2.0);
			jacobian.solveLinearEquationSet(&delta[0]);

			// Update estimate of equilibrium distribution and test for convergence.
			bool testConvergence(true);
			for (size_t j(0); j < eqFraction.size(); ++j) {
				testConvergence = ((delta[j] / eqFraction[j] < tol) && testConvergence);
				eqFraction[j] -= delta[j];
			}
			converged = testConvergence;
		}

		// Normalize equilbrium vector.
		for (size_t j(0); j < eqFraction.size(); ++j) {
			eqFraction[j] *= qd_real(0.5);
		}

		// Calculate the effective inverse, this is for reporting purposes only.
		jacobian = eqMatrix;
		jacobian[idxr][idxs] *= eqFraction[idxs];
		jacobian.invertLUdecomposition();
		eqMatrix = jacobian;

		return true;
	}

	void CollisionOperator::diagReactionOperator(const MesmerFlags &mFlags, const MesmerEnv &mEnv,
		AnalysisData* analysisData)
	{
		// Allocate space for eigenvalues.
		const size_t smsize = m_reactionOperator->size();
		m_eigenvalues.clear();
		m_eigenvalues.resize(smsize, 0.0);
		if (m_eigenvectors) delete m_eigenvectors;
		m_eigenvectors = new qdMatrix(smsize, 0.0);

		// This block prints Reaction Operator before diagonalization
		if (mFlags.printReactionOperatorNum) {
			stest << "Reaction operator --- ";
			printReactionOperator(mFlags);
		}

		//-------------------------------------------------------------
		// diagonalize the whole matrix
		switch (m_precision) {
		case DOUBLE:
		{
			dMatrix dDiagM(smsize);
			for (size_t i = 0; i < smsize; ++i)
				for (size_t j = 0; j < smsize; ++j)
					dDiagM[i][j] = to_double((*m_reactionOperator)[i][j]);
			vector<double>  dEigenValue(smsize, 0.0);
			dDiagM.diagonalize(&dEigenValue[0]);
			for (size_t i = 0; i < smsize; ++i)
				m_eigenvalues[i] = dEigenValue[i];
			for (size_t i = 0; i < smsize; ++i)
				for (size_t j = 0; j < smsize; ++j)
					(*m_eigenvectors)[i][j] = dDiagM[i][j];
			break;
		}
		case DOUBLE_DOUBLE:
		{
			ddMatrix ddDiagM(smsize);
			for (size_t i = 0; i < smsize; ++i)
				for (size_t j = 0; j < smsize; ++j)
					ddDiagM[i][j] = to_dd_real((*m_reactionOperator)[i][j]);
			vector<dd_real> ddEigenValue(smsize, 0.0);
			ddDiagM.diagonalize(&ddEigenValue[0]);
			for (size_t i = 0; i < smsize; ++i)
				m_eigenvalues[i] = ddEigenValue[i];
			for (size_t i = 0; i < smsize; ++i)
				for (size_t j = 0; j < smsize; ++j)
					(*m_eigenvectors)[i][j] = ddDiagM[i][j];
			break;
		}
		default: // diagonalize in quad double
		{
			(*m_eigenvectors) = (*m_reactionOperator);
			m_eigenvectors->diagonalize(&m_eigenvalues[0]);
		}

		}
		// diagonalize the whole matrix
		//-------------------------------------------------------------

		if (mFlags.printEigenValuesNum != 0)
		{
			size_t numberStarted = 0; //will apply when mFlags.printEigenValuesNum<0: print all
			if (mFlags.printEigenValuesNum > 0 && mFlags.printEigenValuesNum <= int(smsize))
				numberStarted = smsize - mFlags.printEigenValuesNum;

			stest << "\nTotal number of eigenvalues = " << smsize << endl;
			stest << "Eigenvalues\n{\n";
			for (size_t i = numberStarted; i < smsize; ++i) {
				qd_real tmp = (mEnv.useBasisSetMethod) ? m_eigenvalues[i] : m_eigenvalues[i] * m_meanOmega;
				formatFloat(stest, tmp, 6, 15);
				stest << endl;
			}
			stest << "}\n";

			if (analysisData) {
				analysisData->m_number = smsize;
				analysisData->m_selection = (mFlags.printEigenValuesNum != -1) ? toString(mFlags.printEigenValuesNum) : "all";
				for (size_t i = numberStarted; i < smsize; ++i) {
					qd_real tmp = (mEnv.useBasisSetMethod) ? m_eigenvalues[i] : m_eigenvalues[i] * m_meanOmega;
					analysisData->m_eigenvalues.push_back(to_double(tmp));
				}
			}

			if (mFlags.printEigenVectors) {
				string title("Eigenvectors:");
				m_eigenvectors->print(title, stest, -1, -1, -1, numberStarted);
			}
		}

	}

	void CollisionOperator::printReactionOperator(const MesmerFlags& mFlags)
	{
		const int smsize = int(m_reactionOperator->size());

		switch (mFlags.printReactionOperatorNum)
		{
		case -1:
			stest << "Printing all (" << smsize << ") columns/rows of the Reaction Operator:\n";
			(*m_reactionOperator).showFinalBits(smsize, mFlags.print_TabbedMatrices);
			break;
		case -2:
			stest << "Printing final 1/2 (" << smsize / 2 << ") columns/rows of the Reaction Operator:\n";
			(*m_reactionOperator).showFinalBits(smsize / 2, mFlags.print_TabbedMatrices);
			break;
		case -3:
			stest << "Printing final 1/3 (" << smsize / 3 << ") columns/rows of the Reaction Operator:\n";
			(*m_reactionOperator).showFinalBits(smsize / 3, mFlags.print_TabbedMatrices);
			break;
		default: // the number is either smaller than -3 or positive
			if (abs(mFlags.printReactionOperatorNum) > smsize) {
				stest << "Printing all (" << smsize << ") columns/rows of the Reaction Operator:\n";
				(*m_reactionOperator).showFinalBits(smsize, mFlags.print_TabbedMatrices);
			}
			else {
				stest << "Printing final " << abs(mFlags.printReactionOperatorNum) << " columns/rows of the Reaction Operator:\n";
				(*m_reactionOperator).showFinalBits(abs(mFlags.printReactionOperatorNum), mFlags.print_TabbedMatrices);
			}
		}
	}

	bool CollisionOperator::timeEvolution(MesmerFlags& mFlags, AnalysisData* analysisData)
	{
		ErrorContext c(__FUNCTION__);

		// Cut short if species profiles not needed.
		if (!mFlags.speciesProfileEnabled)
			return true;

		size_t smsize = m_eigenvectors->size();
		vector<double> r_0(smsize, 0.);
		if (!projectedInitialDistrbtn(r_0)) {
			cerr << "Projection of initial disttribution failed.";
			return false;
		}

		// Calculate time points of interest.
		vector<double> timePoints;
		timeAxisPoints(mFlags, timePoints);
		if (analysisData)
			analysisData->m_timePoints = timePoints;

		dMatrix totalEigenVecs(smsize); // copy full eigenvectors of the system
		for (size_t i = 0; i < smsize; ++i) {
			double tmp = to_double(m_eqVector[i]);
			for (size_t j = 0; j < smsize; ++j) {
				totalEigenVecs[i][j] = tmp*to_double((*m_eigenvectors)[i][j]);
			}
		}

		const size_t maxTimeStep = timePoints.size();
		vector<vector<double> > grnProfile(smsize, vector<double>(maxTimeStep));
		vector<double> p_t(smsize, 0.);

		for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
			double numColl = m_meanOmega * timePoints[timestep];
			for (size_t j(0); j < smsize; ++j) {
				p_t[j] = r_0[j] * exp(to_double(m_eigenvalues[j]) * numColl);
			} // now |p_t> = exp(Lambda*t)*V^(T)*|init> = exp(Lambda*t)*U^(-1)*|n_0>

			p_t *= totalEigenVecs;

			for (size_t j(0); j < smsize; ++j) {
				grnProfile[j][timestep] = p_t[j];
			} // now |grnProfile(t)> = |grnProfile(i)> = F*V*exp(Lambda*t)*V^(T)*|init> = U*exp(Lambda*t)*U^(-1)*|n_0>

		}

		//------------------------------
		// print grained species profile
		if (mFlags.grainedProfileEnabled) {
			stest << "\nGrained species profile:(first row is the time point in units of second & first column is the grain index)\n{\n";
			Reaction::molMapType::iterator ipos;
			// Iterate through isomer map to print out which grains are spanned by which isomers.
			for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
				Molecule* isomer = ipos->first;
				stest << " isomer " << isomer->getName() << " spans grains " << ipos->second << " to "
					<< ipos->second + isomer->getColl().get_colloptrsize() - 1 << endl;  // offset of 1 is b/c grain idx starts at 0
			}

			Reaction::molMapType::iterator spos;
			for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {  // iterate through source map
				Molecule* source = spos->first;                        // to print out which grains are spanned by which sources
				stest << " source " << source->getName() << " is in grain " << spos->second << endl;
			}

			stest << "\n\t";
			for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
				formatFloat(stest, timePoints[timestep], 6, 15);
				stest << ",";
			}
			stest << endl;
			for (size_t j(0); j < smsize; ++j) {
				stest << j << ","; // << "\t";
				for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
					formatFloat(stest, grnProfile[j][timestep], 6, 15);
					stest << ",";
				}
				stest << endl;
			}

			// Now print out the average of the grain energy in each isomer.
			stest << endl << "average energy in each isomer (kJ/mol)" << endl;
			for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
				Molecule* isomer = ipos->first;
				stest << isomer->getName() << ": " << endl;
				stest << "\t" << endl;
				if (analysisData)
					analysisData->m_aveEnergyRef.push_back(isomer->getName());
				vector<double> grnEne;
				isomer->getDOS().getGrainEnergies(grnEne);
				for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
					double averageEnergy(0.0);
					double totalIsomerPopulation(0.0);
					for (size_t grain(ipos->second), iene(0); iene < (size_t)isomer->getColl().get_colloptrsize(); ++iene, ++grain) {
						double pop = grnProfile[grain][timestep];
						totalIsomerPopulation += pop;       // Determine how much total population in each isomer at time t.
						averageEnergy += grnEne[iene] * pop;	// Calculate the average energy in each isomer at time t.
					}
					double normAvE = (totalIsomerPopulation > 0.0) ? ConvertFromWavenumbers("kJ/mol", averageEnergy / totalIsomerPopulation) : 0.0;
					formatFloat(stest, timePoints[timestep], 6, 15);
					formatFloat(stest, normAvE, 6, 15);
					stest << endl;
					if (analysisData)
						analysisData->m_aveEnergy.push_back(normAvE);
				}
				stest << endl;
			}

			stest << "}\n";
		}

		stest << "mean collision frequency = " << m_meanOmega << "/s" << endl;

		vector<double> totalIsomerPop(maxTimeStep, 0.);
		vector<double> totalPdtPop(maxTimeStep, 0.);

		for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
			for (size_t j(0); j < smsize; ++j) {
				totalIsomerPop[timestep] += grnProfile[j][timestep];
			}
			double popTime = totalIsomerPop[timestep];
			if (popTime > 1.0) {
				popTime = 1.0; // correct some numerical error
				//totalIsomerPop[timestep] = 1.0; // Not very sure if we should cover up this numerical error entirely!!?
			}
			else if (popTime < 0.0) {
				popTime = 0.0;
				//totalIsomerPop[timestep] = 0.0; // Not very sure if we should cover up this numerical error entirely!!?
			}
			totalPdtPop[timestep] = 1.0 - popTime;
		}

		if (mFlags.speciesProfileEnabled) {
			stest << endl << "Print time dependent species and product profiles" << endl << "{" << endl;
			int numberOfSpecies = static_cast<int>(m_isomers.size() + m_sources.size() + m_sinkRxns.size());

			vector<vector<double> > speciesProfile(numberOfSpecies, vector<double>(maxTimeStep));
			int speciesProfileidx(0);

			stest << setw(16) << "Timestep/s";

			// Iterate through the source map, to determine the total source 
			// density as a function of time.
			vector<string> speciesNames;
			Reaction::molMapType::iterator spos;
			for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {
				Molecule* source = spos->first;
				stest << setw(16) << source->getName();
				speciesNames.push_back(source->getName());
				int rxnMatrixLoc = spos->second;
				for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
					speciesProfile[speciesProfileidx][timestep] = grnProfile[rxnMatrixLoc][timestep];
				}
				++speciesProfileidx;
			}

			// Iterate through the isomer map, to calculate the total isomer 
			// density as a function of time.
			Reaction::molMapType::iterator ipos;
			for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
				Molecule* isomer = ipos->first;
				string isomerName = isomer->getName();
				stest << setw(16) << isomerName;
				speciesNames.push_back(isomerName);
				int rxnMatrixLoc = ipos->second;
				size_t colloptrsize = isomer->getColl().get_colloptrsize();
				for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
					for (size_t i(0); i < colloptrsize; ++i) {
						speciesProfile[speciesProfileidx][timestep] += grnProfile[i + rxnMatrixLoc][timestep];
					}
				}
				++speciesProfileidx;
			}

			// SHR 2/Jan/2015: The following vector holds the flux into each sink.
			// I am not sure how useful this quantity is and may wish to remove it 
			// at a later stage.
			vector<vector<double> > sinkFluxProfile(m_sinkRxns.size(), vector<double>(maxTimeStep));

			int pdtProfileStartIdx = speciesProfileidx;
			size_t fluxIdx(0);

			// Iterate through sink map to get product profile vs t.
			// SHR 26/Mar/2017: This section involves a very crude integration 
			// over time and is in need of review.
			sinkMap::iterator pos;
			for (pos = m_sinkRxns.begin(); pos != m_sinkRxns.end(); ++pos) {
				Reaction* sinkReaction = pos->first;
				const vector<double> KofEs = sinkReaction->get_MtxGrnKf();        // Vector to hold sink k(E)s.
				vector<Molecule*> pdts;                                           // in the sink reaction
				sinkReaction->get_products(pdts);
				string pdtName = pdts[0]->getName();
				if (KofEs.size() == 1) {
					pdtName += "(bim)";
				}
				stest << setw(16) << pdtName;
				speciesNames.push_back(pdtName);
				int rxnMatrixLoc = pos->second;  // Get sink location.
				double TimeIntegratedProductPop(0.0);

				double lastTime(0.0);
				for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
					double deltat = timePoints[timestep] - lastTime;
					lastTime = timePoints[timestep];
					for (size_t i(0); i < KofEs.size(); ++i) {
						speciesProfile[speciesProfileidx][timestep] += KofEs[i] * grnProfile[i + rxnMatrixLoc][timestep] * deltat;
						sinkFluxProfile[fluxIdx][timestep] += KofEs[i] * grnProfile[i + rxnMatrixLoc][timestep];
					}
					TimeIntegratedProductPop += speciesProfile[speciesProfileidx][timestep];
					speciesProfile[speciesProfileidx][timestep] = TimeIntegratedProductPop;
				}
				++speciesProfileidx;
				++fluxIdx;
			}

			if (pdtProfileStartIdx < speciesProfileidx) {
				for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {    // normalize product profile to account for small
					double normConst(0.0);                          // numerical errors in TimeIntegratedProductPop
					double pdtYield(0.0);
					for (int i(pdtProfileStartIdx); i < speciesProfileidx; ++i) {   // calculate normalization constant
						pdtYield += speciesProfile[i][timestep];
					}
					normConst = totalPdtPop[timestep] / pdtYield;
					for (int i(pdtProfileStartIdx); i < speciesProfileidx; ++i) {   // apply normalization constant
						speciesProfile[i][timestep] *= normConst;
					}
				}
			}

			// Write to stest and (ultimately) XML.
			stest << setw(16) << "totalIsomerPop" << setw(16) << "totalPdtPop" << endl;
			if (analysisData) {
				for (int i(0); i < speciesProfileidx; ++i)
					analysisData->m_PopRef.push_back(speciesNames[i]);
			}
			for (size_t timestep(0); timestep < maxTimeStep; ++timestep) {
				stest << setw(16) << timePoints[timestep];
				for (int i(0); i < speciesProfileidx; ++i) {
					stest << setw(16) << speciesProfile[i][timestep];
					if (analysisData) {
						analysisData->m_Pop.push_back(speciesProfile[i][timestep]);
					}
				}
				stest << setw(16) << totalIsomerPop[timestep] << setw(16) << totalPdtPop[timestep];
				if (mFlags.printSinkFluxes) {
					for (size_t i(0); i < m_sinkRxns.size(); ++i) {
						stest << setw(16) << sinkFluxProfile[i][timestep];
					}
				}
				stest << endl;
			}
			stest << "}" << endl;
		}
		return true;
	}

	// Method to calculate points of interest on the time axis.
	bool CollisionOperator::timeAxisPoints(MesmerFlags& mFlags, vector<double>& timePoints) {
		double shortestTime = 0.;
		// set the default maximum evolution time
		if (mFlags.shortestTimeOfInterest < 1.0e-20 || mFlags.shortestTimeOfInterest > 1.0)
			shortestTime = 1.0e-11;
		else
			shortestTime = mFlags.shortestTimeOfInterest;

		double maxEvoTime = 0.;
		// set the default maximum evolution time
		if (mFlags.maxEvolutionTime <= 0.001 || mFlags.maxEvolutionTime > 1.0e8)
			maxEvoTime = 1.2e5;
		else
			maxEvoTime = mFlags.maxEvolutionTime;

		// Calculates the time points
		for (int i = 0; i <= 300; ++i) {
			double thetime = pow(10., static_cast<double>(i) / 10. - 20.);
			if (thetime < shortestTime) continue;
			if (thetime > maxEvoTime) break;
			timePoints.push_back(thetime);
		}

		return true;
	}

	bool CollisionOperator::produceEquilibriumVector()
	{

		m_eqVector.clear();
		m_eqVector.resize(m_eqVecSize);

		Reaction::molMapType::iterator spos;
		for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {  // Iterate through the source map to get
			Molecule* source = spos->first;                                 // the equilibrum fractions.
			int rxnMatrixLoc = spos->second;
			qd_real eqFrac = source->getPop().getEqFraction();
			m_eqVector[rxnMatrixLoc] = sqrt(eqFrac);
		}

		Reaction::molMapType::iterator ipos;
		for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {  // Iterate through the isomer map
			Molecule* isomer = ipos->first;                                 // to get the equilibrium fractions.
			int rxnMatrixLoc = ipos->second;
			qd_real eqFrac = isomer->getPop().getEqFraction();
			const size_t colloptrsize = isomer->getColl().get_colloptrsize();
			vector<double> boltzFrac;
			isomer->getColl().normalizedGrnBoltzmannDistribution(boltzFrac);
			for (size_t i(0); i < colloptrsize; ++i) {
				m_eqVector[rxnMatrixLoc + i] = sqrt(eqFrac * qd_real(boltzFrac[i]));
			}
		}
		return true;
	}

	bool CollisionOperator::produceInitialPopulationVector(vector<double>& n_0) const {

		double populationSum = 0.0;

		Reaction::molMapType::const_iterator ipos;
		for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {  // iterate through isomer map
			Molecule* isomer = ipos->first;                        // to get isomer initial populations
			populationSum += isomer->getPop().getInitPopulation();
		}

		Reaction::molMapType::const_iterator spos;
		for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {  // iterate through source map to get
			Molecule* source = spos->first;                         // source initial populations
			populationSum += source->getPop().getInitPopulation();
		}

		for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
			Molecule* isomer = ipos->first;                        // get initial population of each isomer
			double initFrac = isomer->getPop().getInitPopulation();
			if (initFrac != 0.0) {                                           // if isomer initial populations are nonzero
				initFrac /= populationSum;                                    // normalize initial pop fraction
				int rxnMatrixLoc = ipos->second;
				const size_t colloptrsize = isomer->getColl().get_colloptrsize();

				map<int, double> grainMap;                            // get the grain pop map and check to see if any grain populations are specified
				isomer->getPop().getInitGrainPopulation(grainMap);

				if (grainMap.size() != 0) {   // if grain populations have been specified, then use them
					for (size_t i(0); i < colloptrsize; ++i) {
						n_0[i + rxnMatrixLoc] = 0.0;  // set elements of initial distribution vector to zero
					}
					map<int, double>::iterator grainIt;
					for (grainIt = grainMap.begin(); grainIt != grainMap.end(); ++grainIt) {
						if ((grainIt->first) < int(n_0.size())) {
							n_0[grainIt->first + rxnMatrixLoc - 1] = initFrac * grainIt->second;    // put populations in grain n where n=GrainMap->first	
						}
						else {
							cerr << "you requested population in grain " << grainIt->first << " of isomer " << isomer->getName() << ", which exceeds the number of grains in the isomer" << endl;
							cerr << "you must respecify the system to accomodate your request... exiting execution " << endl;
							exit(1);
						}
					}
				}
				else {	// otherwise if no grain population has been specified, use a boltzmann population
					vector<double> boltzFrac;
					isomer->getColl().normalizedInitialDistribution(boltzFrac);
					for (size_t i(0); i < colloptrsize; ++i) {
						n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
					}
				}
			}
		}

		// if there is no source term and the populationSum is still zero, set population = 1.0 for the first isomer
		int sizeSource = static_cast<int>(m_sources.size());
		if (populationSum == 0. && sizeSource == 0) {
			ipos = m_isomers.begin();
			Molecule* isomer = ipos->first;
			isomer->getPop().setInitPopulation(1.0); // set initial population for the first isomer
			double initFrac = isomer->getPop().getInitPopulation();
			cinfo << "No population was assigned, and there is no source term." << endl
				<< "Initial poupulation set to 1.0  in the first isomer." << endl;
			int rxnMatrixLoc = ipos->second;
			const size_t colloptrsize = isomer->getColl().get_colloptrsize();
			vector<double> boltzFrac;
			isomer->getColl().normalizedInitialDistribution(boltzFrac);
			for (size_t i(0); i < colloptrsize; ++i) {
				n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
			}
		}

		for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {
			Molecule* source = spos->first;
			int rxnMatrixLoc = spos->second;
			if (populationSum == 0. && spos == m_sources.begin()) {
				cinfo << "No population was assigned. Initialize the first source term to 1.0." << once << endl;
				n_0[rxnMatrixLoc] = 1.0;
			}
			else {
				double initFrac = source->getPop().getInitPopulation() / populationSum;
				n_0[rxnMatrixLoc] = initFrac;
			}
		}

		return true;
	}

	bool CollisionOperator::BartisWidomPhenomenologicalRates(qdMatrix& mesmerRates, qdMatrix& lossRates, MesmerFlags& mFlags, AnalysisData* anlysisData)
	{
		// Constants.
		const size_t smsize = m_eigenvectors->size();
		const size_t nchem = m_isomers.size() + m_sources.size();  // Number of isomers+pseudoisomers.
		const size_t nchemIdx = smsize - nchem;       // Location of chemically significant eigenvalues & vectors.
		const size_t nsinks = m_sinkRxns.size();    // Number of Sinks.

		stest << "\nBartis Widom eigenvalue/eigenvector analysis\n" << endl;
		stest << "Number of sinks in this system: " << nsinks << endl;

		if (nsinks > 0) {
			stest << "\nThere should be " << nchem << " chemically significant eigenvalues (CSEs)" << endl;
		}
		else {
			stest << "\nThere should be 1 zero eigenvalue (zero within numerical precision) and " << nchem - 1
				<< " chemically significant eigenvalues (CSEs)" << endl;
		}

		//
		// If there are no sinks, replace the equilibrium vector with the eigenvector whose
		// associated eigenvalue is zero, as this is a consistent estimate of the equilibrium 
		// with respect to the other eigenvalues. Also, as the system is conservative, set the 
		// smallest eigenvalue explicitly to zero.
		//
		if (nsinks < 1) {
			m_eigenvalues[smsize - 1] = 0.0;
			for (size_t i(0); i < smsize; ++i) {
				m_eqVector[i] = (*m_eigenvectors)[i][smsize - 1];
			}
		}

		// Construct assymmetric eigenvectors required for the z matrix.

		vector<vector<qd_real> > assymEigenVec(nchem, vector<qd_real>(smsize, 0.0)); // U
		for (size_t i(0), ii(nchemIdx); i < nchem; ++i, ++ii) {
			for (size_t j(0); j < smsize; ++j) {
				assymEigenVec[i][j] = m_eqVector[j] * (*m_eigenvectors)[j][ii]; // Calculation of U = FV
			}
		}

		// Check that the inverse matrix is correctly calculated by multiplying U*U^(-1) 
		// for CSE vectors.

		for (size_t i(0); i < nchem; ++i) {
			qd_real test = 0.0;
			for (size_t j(nchemIdx); j < smsize; ++j) {
				qd_real sm = 0.0;
				for (size_t k(0); k < smsize; ++k) {
					sm += assymEigenVec[i][k] * (*m_eigenvectors)[k][j] / m_eqVector[k]; // Calculation of U^(-1) = (FV)^-1 = V^T * F^-1
				}
				test += sm;
			}
			if (fabs(test - 1.0) > 0.001)      // test that U*U^(-1) = 1
				stest << "row " << i << " of the U*U^(-1) matrix does not equal unity. It sums to " << test << endl;
		}

		if (!mFlags.rateCoefficientsOnly) {
			qdMatrix Z_matrix(nchem);  // definitions of Y_matrix and Z_matrix taken from PCCP 2007(9), p.4085
			qdMatrix Y_matrix(max(nchem, nsinks));
			Reaction::molMapType::iterator ipos;  // Set up an iterator through the isomer map.
			Reaction::molMapType::iterator spos;  // Set up an iterator through the source map.
			sinkMap::iterator sinkpos;            // Set up an iterator through the irreversible rxn map.

			// Check the separation between chemically significant eigenvalues (CSEs) and
			// internal energy relaxation eigenvalues (IEREs); if it's not good, print a warning.

			int adsorbedCSE(0);
			const double adsorbedCSETol = 0.1;
			const double last_CSE = (to_double(m_eigenvalues[nchemIdx]));
			const double first_IERE = (to_double(m_eigenvalues[nchemIdx - 1]));
			const double CSE_IERE_ratio = last_CSE / first_IERE;
			static bool bCSE_IERE_ratio_WARN = true;
			if (CSE_IERE_ratio > adsorbedCSETol && bCSE_IERE_ratio_WARN) {
				bCSE_IERE_ratio_WARN = false; // Only issue this warning once.
				stringstream ss1;
				ss1 << "\nWARNING: Chemically significant eigenvalues (CSE) not well separated from internal energy relaxation eigenvals (IEREs)." << endl;
				ss1 << "\nThe last CSE = " << last_CSE*m_meanOmega << " and the first IERE = " << first_IERE*m_meanOmega << endl;
				ss1 << "(last CSE)/(first IERE) ratio = " << CSE_IERE_ratio << ", which is less than an order of magnitude" << endl;
				ss1 << "\nResults obtained from Bartis Widom eigenvalue-vector analysis may be unreliable" << endl;
				string s(ss1.str());
				stest << s; clog << s;

				// Replace tabs and line feeds (bad for XML) by spaces.
				replace(s.begin(), s.end(), '\t', ' ');
				replace(s.begin(), s.end(), '\n', ' ');
				if (anlysisData) anlysisData->m_warning = s;

				// Determine the number of adsorbed eigenvalues.
				adsorbedCSE++;
				for (size_t i(1); (i < nchem) && (to_double(m_eigenvalues[nchemIdx + i]) / first_IERE > adsorbedCSETol); i++) {
					adsorbedCSE++;
				}
			}

			for (size_t i(0); i < nchem; ++i) {

				// Calculate Z matrix elements for all the isomers in the system.

				for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
					qd_real sm(0.0);
					Molecule* isomer = ipos->first;
					size_t colloptrsize = isomer->getColl().get_colloptrsize(); // get colloptrsize for isomer
					int rxnMatrixLoc = ipos->second + colloptrsize - 1;         // get location for isomer in the rxn matrix
					int seqMatrixLoc = m_SpeciesSequence[isomer];               // get sequence position for isomer
					for (size_t j(0); j < colloptrsize; ++j) {
						sm += assymEigenVec[i][rxnMatrixLoc - j];
					}
					Z_matrix[seqMatrixLoc][i] = sm;
				}

				// Calculate Z_matrix matrix elements for all sources in the system.

				for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {
					Molecule* pPseudoIsomer = spos->first;
					const int rxnMatrixLoc = spos->second;
					const int seqMatrixLoc = m_SpeciesSequence[pPseudoIsomer];
					Z_matrix[seqMatrixLoc][i] = assymEigenVec[i][rxnMatrixLoc];
				}

				// Calculate Y_matrix elements for sinks.

				if (nsinks) {
					int seqMatrixLoc(0);
					for (sinkpos = m_sinkRxns.begin(); sinkpos != m_sinkRxns.end(); ++sinkpos, ++seqMatrixLoc) {
						Reaction* sinkReaction = sinkpos->first;
						const vector<double> KofEs = sinkReaction->get_MtxGrnKf();  // Vector to hold sink k(E)s.
						int rxnMatrixLoc = sinkpos->second;                         // Get sink location.
						qd_real sm(0.0);
						for (size_t j(0); j < KofEs.size(); ++j) {
							sm += assymEigenVec[i][rxnMatrixLoc + j] * KofEs[j];
						}
						Y_matrix[seqMatrixLoc][i] = sm;
					}
				}

			}

			// Print out Y_matrix for testing.
			if (nsinks) {
				string MatrixTitle("Y_matrix:");
				Y_matrix.print(MatrixTitle, stest, int(nsinks), int(m_SpeciesSequence.size()));
			}

			qdMatrix Zinv(Z_matrix);
			if (nsinks && !mFlags.bForceMacroDetailedBalance) {

				// Apply standard inversion method.

				if (Zinv.invertLUdecomposition()) {
					cerr << "Inversion of Z_matrix failed.  Matrix before inversion is: ";
					Z_matrix.showFinalBits(nchem);
				}

			}
			else {

				// Apply Gram-Schmit orthogonalization in order to invert the matrix.
				// This imposes detailed balance at the macroscopic level.
				//
				// SHR 25/Apr/2010 : It remains unclear that this is correct at the time
				// of writting, however for some systems it is difficult to realize mass
				// conservation without it.
				//
				// SHR 25/Aug/2013 : The above comment refers to conservative systems.
				// The method has been extended to non-conservative systems by using the
				// equilibrium distribution calculated previously. It is even less clear
				// that this is appropriate, as microscopic reversibility is explicitly
				// broken. However, in situations where -ve rate coefficients are observed,
				// this method can rectify the problem.

				// Decompose the reduced eigenvector matrix.

				qdMatrix Fr(nchem), Fr_inv(nchem);

				if (nsinks && mFlags.bForceMacroDetailedBalance) {

					// Non-conservative case.

					Reaction::molMapType::iterator spcitr = m_SpeciesSequence.begin();
					for (; spcitr != m_SpeciesSequence.end(); ++spcitr) {
						size_t i = spcitr->second;
						Fr[i][i] = sqrt(qd_real((spcitr->first)->getPop().getEqFraction()));
						Fr_inv[i][i] = qd_real(1.0) / Fr[i][i];
					}
				}
				else {

					// Conservative case.

					for (size_t i(0); i < nchem; ++i) {
						Fr[i][i] = sqrt(Z_matrix[i][nchem - 1]);
						Fr_inv[i][i] = qd_real(1.0) / Fr[i][i];
					}
				}


				qdMatrix Er = Fr_inv * Z_matrix;

				// Orthogonalize the reduced symmetric eigenvectro matrix.

				Er.GramSchimdt(nchem - 1);

				Z_matrix = Fr * Er;

				// Transpose the orthonormal matrix and form inverse.

				Er.Transpose();

				Zinv = Er * Fr_inv;

			}

			stest << "\nZ_matrix: ";
			Z_matrix.showFinalBits(nchem, true);

			stest << endl << "Z_matrix^(-1):" << endl;
			Zinv.showFinalBits(nchem, true);

			qdMatrix Zidentity = Z_matrix * Zinv;

			stest << "\nZ_matrix * Z_matrix^(-1) [Identity matrix]:" << endl;
			Zidentity.showFinalBits(nchem, true);

			// Construct phenomenological rate coefficient matrix.

			qdMatrix Egv(nchem);
			for (size_t i(0); i < nchem; ++i) {
				Egv[i][i] = m_eigenvalues[nchemIdx + i] * m_meanOmega;
			}
			qdMatrix Kr = Z_matrix * Egv * Zinv;

			stest << "\nKr matrix:" << endl;
			Kr.showFinalBits(nchem, true);       // Print out Kr_matrix

			// Construct loss matrix.

			qdMatrix Kp(max(nsinks, nchem), 0.0);
			if (nsinks > 0) {
				for (size_t i(0); i != nsinks; ++i) {    // calculate Kp (definition taken from PCCP 2007(9), p.4085)
					for (size_t j(0); j < nchem; ++j) {
						qd_real sm = 0.0;
						for (size_t k(0); k < nchem; ++k) {
							sm += Y_matrix[i][k] * Zinv[k][j];
						}
						Kp[i][j] = sm;
					}
				}
				string MatrixTitle("Kp matrix:");
				Kp.print(MatrixTitle, stest, nsinks, m_SpeciesSequence.size());
			}

			// If requested, write out phenomenological evolution.
			if (mFlags.printPhenomenologicalEvolution) {
				if (mFlags.bIsSystemSecondOrder) {
					cinfo << "At present it is not possible to phenomenological profiles for systems with a second order term." << once << endl;
				}
				else {
					PhenomenologicalIntegration(Z_matrix, Zinv, Egv, mFlags);
				}
			}

			// If eigenvalues adsorbed calculate an effective set of rate equations.
			if (adsorbedCSE > 0) {
			}

			mesmerRates = Kr;
			lossRates = Kp;
		}
		return true;

	}

	// Method to integrate the phenomenological rate equations using BW coefficients.
	bool CollisionOperator::PhenomenologicalIntegration(qdMatrix& Z_matrix, qdMatrix& Zinv, qdMatrix& Egv, MesmerFlags& mFlags) {

		// Write header section and setup initial concentration vector. 

		stest << endl << "Phenomenological species profiles" << endl << "{" << endl;
		stest << setw(16) << "Timestep/s";

		size_t nchem = Z_matrix.size();
		vector<qd_real> c0(nchem, 0.0);
		vector<size_t> speciesOrder;

		// Loop over isomers then sources.

		Reaction::molMapType::iterator ipos;
		for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
			Molecule* isomer = ipos->first;
			stest << setw(16) << isomer->getName();
			int seqMatrixLoc = m_SpeciesSequence[isomer];
			speciesOrder.push_back(size_t(seqMatrixLoc));
			c0[seqMatrixLoc] = isomer->getPop().getInitPopulation();
		}

		for (ipos = m_sources.begin(); ipos != m_sources.end(); ++ipos) {
			Molecule* pseudoisomer = ipos->first;
			stest << setw(16) << pseudoisomer->getName();
			int seqMatrixLoc = m_SpeciesSequence[pseudoisomer];
			speciesOrder.push_back(size_t(seqMatrixLoc));
			c0[seqMatrixLoc] = pseudoisomer->getPop().getInitPopulation();
		}
		stest << endl;

		// Calculate time points, calculate and write populations. 

		c0 *= Zinv;
		vector<double> timePoints;
		timeAxisPoints(mFlags, timePoints);
		for (size_t i(0); i < timePoints.size(); ++i) {
			qd_real time = timePoints[i];
			vector<qd_real> p(Z_matrix.size(), 0.0);
			for (size_t j(0); j < p.size(); ++j) {
				p[j] = exp(Egv[j][j] * time)*c0[j];
			}
			p *= Z_matrix;
			stest << setw(16) << time;
			for (size_t j(0); j < p.size(); ++j) {
				stest << setw(16) << p[speciesOrder[j]];
			}
			stest << endl;
		}
		stest << "}" << endl;

		return true;
	}

	// Write out phenomenological rate coefficients.
	bool CollisionOperator::PrintPhenomenologicalRates(qdMatrix& Kr, qdMatrix& Kp, MesmerFlags& mFlags, AnalysisData* analysisData) {

		Reaction::molMapType::iterator ipos;  // set up an iterator through the isomer map

		stest << "\nFirst order & pseudo first order rate coefficients for loss rxns:\n{\n";
		Reaction::molMapType::iterator lossitr, rctitr, pdtitr;

		stringstream puSymbols;
		stringstream puNumbers;
		// print pseudo 1st order k loss for isomers
		for (lossitr = m_SpeciesSequence.begin(); lossitr != m_SpeciesSequence.end(); ++lossitr) {
			Molecule* iso = lossitr->first;
			int losspos = lossitr->second;
			string isomerName = iso->getName();
			stest << isomerName << " loss = " << Kr[losspos][losspos] << endl;
			if (analysisData) {
				analysisData->m_lossRef.push_back(isomerName);
				analysisData->m_lossRateCoeff.push_back(to_double(Kr[losspos][losspos]));
			}
			puNumbers << Kr[losspos][losspos] << "\t";
			if (m_punchSymbolGathered == false) {
				puSymbols << isomerName << " loss\t";
			}
		}
		stest << "}\n";

		if (m_SpeciesSequence.size() > 1) {
			stest << "\nFirst order & pseudo first order rate coefficients for isomerization rxns:\n{\n";

			// print pseudo first order connecting ks
			for (rctitr = m_SpeciesSequence.begin(); rctitr != m_SpeciesSequence.end(); ++rctitr) {
				string rctName = rctitr->first->getName();
				int rctpos = rctitr->second;
				for (pdtitr = m_SpeciesSequence.begin(); pdtitr != m_SpeciesSequence.end(); ++pdtitr) {
					string pdtName = pdtitr->first->getName();
					int pdtpos = pdtitr->second;
					if (rctpos != pdtpos) {
						stest << rctName << " -> " << pdtName << " = " << Kr[pdtpos][rctpos] << endl;

						ostringstream reaction;
						reaction << rctName << " => " << pdtName;
						m_phenomenlogicalRates[reaction.str()] = to_double(Kr[pdtpos][rctpos]);

						if (analysisData) {
							analysisData->m_firstOrderRateCoeff.push_back(to_double(Kr[pdtpos][rctpos]));
							analysisData->m_firstOrderFromRef.push_back(rctName);
							analysisData->m_firstOrderToRef.push_back(pdtName);
							analysisData->m_firstOrderReactionType.push_back("isomerization");
						}
					}

					puNumbers << Kr[pdtpos][rctpos] << "\t";
					if (m_punchSymbolGathered == false) {
						puSymbols << rctName << " -> " << pdtName << "\t";
					}
				}
			}
			stest << "}\n";
		}

		if (m_sinkRxns.size() != 0) {
			stest << "\nFirst order & pseudo first order rate coefficients for irreversible rxns:\n{\n";
			sinkMap::iterator sinkitr = m_sinkRxns.begin();

			for (int sinkpos(0); sinkitr != m_sinkRxns.end(); ++sinkitr, ++sinkpos) {
				Reaction* sinkReaction = sinkitr->first;          // get Irreversible Rxn
				vector<Molecule*> pdts;
				sinkReaction->get_products(pdts);
				string pdtsName = pdts[0]->getName();
				if (pdts.size() >= 2) { pdtsName += "+"; pdtsName += pdts[1]->getName(); }
				if (pdts.size() >= 3) { pdtsName += "+"; pdtsName += pdts[2]->getName(); }
				for (rctitr = m_SpeciesSequence.begin(); rctitr != m_SpeciesSequence.end(); ++rctitr) {
					Molecule* rcts = rctitr->first;     // get reactants & their position
					int rctpos = rctitr->second;
					string rctName = rcts->getName();
					if (sinkReaction->getReactionType() == IRREVERSIBLE_EXCHANGE) {
						stest << rctName << " -> " << pdtsName << "(bim) = " << Kp[sinkpos][rctpos] << endl;

						ostringstream reaction;
						reaction << rctName << " => " << pdtsName;
						m_phenomenlogicalRates[reaction.str()] = to_double(Kp[sinkpos][rctpos]);

						puNumbers << Kp[sinkpos][rctpos] << "\t";
						if (!m_punchSymbolGathered) {
							puSymbols << rcts->getName() << " -> " << pdtsName << "(bim)\t";
						}
					}
					else {
						stest << rctName << " -> " << pdtsName << " = " << Kp[sinkpos][rctpos] << endl;

						ostringstream reaction;
						reaction << rctName << " => " << pdtsName;
						m_phenomenlogicalRates[reaction.str()] = to_double(Kp[sinkpos][rctpos]);

						if (analysisData) {
							analysisData->m_firstOrderRateCoeff.push_back(to_double(Kp[sinkpos][rctpos]));
							analysisData->m_firstOrderFromRef.push_back(rctName);
							analysisData->m_firstOrderToRef.push_back(pdtsName);
							analysisData->m_firstOrderReactionType.push_back("irreversible");
							puNumbers << Kp[sinkpos][rctpos] << "\t";
						}
						if (m_punchSymbolGathered == false) {
							puSymbols << rctName << " -> " << pdtsName << "\t";
						}
					}
				}
			}
			stest << "}\n\n";
		}

		if (puSymbols.str().size()) {
			puSymbols << "\n";
			mFlags.punchSymbols = puSymbols.str();
			m_punchSymbolGathered = true;
		}

		if (puNumbers.str().size()) {
			puNumbers << "\n";
			mFlags.punchNumbers = puNumbers.str();
		}

		return true;
	}

	//
	// Calculates the Bartis-Widom macroscopic rate coefficients, using the contracted basis set eigenvectors.
	//
	bool CollisionOperator::BartisWidomBasisSetRates(qdMatrix& mesmerRates, MesmerFlags& mFlags) {

		// Constants.
		const size_t smsize = m_eigenvectors->size();
		const size_t nchem = m_isomers.size() + m_sources.size();  // number of isomers+pseudoisomers
		// const size_t nchemIdx = smsize - nchem ;                       // idx for chemically significant eigenvalues & vectors

		// Print out eigenvector matrix.

		//stest << endl << "Eigenvector matrix:" << endl << endl ;
		//for (size_t i(0) ; i < smsize ; ++i) {
		//  for (size_t j(0) ; j < smsize ; ++j) {
		//    formatFloat(stest, (*m_eigenvectors)[i][j],  6,  15) ;
		//  }
		//  stest << endl ;
		//}

		qdMatrix Z(nchem), Zinv(nchem), Kr(nchem);

		if (m_sinkRxns.size() == 0) {

			//
			// Conservative system.
			//

			// 1. Isomers.

			size_t location(0);
			Reaction::molMapType::iterator isomeritr = m_isomers.begin();
			for (size_t i(0); isomeritr != m_isomers.end(); ++i, ++isomeritr) {
				location = isomeritr->second;
				for (size_t j(1); j <= nchem; ++j) {
					Z[i][nchem - j] = (*m_eigenvectors)[location][smsize - j];
				}
			}

			// Invert Z matrix. 

			stest << endl << "BW coefficient matrix:" << endl << endl;
			for (size_t i(0); i < nchem; ++i) {
				for (size_t j(0); j < nchem; ++j) {
					formatFloat(stest, Z[i][j], 6, 15);
					Zinv[j][i] = Z[i][j];
				}
				stest << endl;
			}

			// Calculate symmetric rate matrix.

			m_eigenvalues[smsize - 1] = 0.0;

			for (size_t i(0); i < nchem; ++i) {
				for (size_t j(0); j < nchem; ++j) {
					qd_real sm = 0.0;
					for (size_t k(0); k < nchem; ++k) {
						// sm += Z[i][k] * to_double(m_eigenvalues[nchemIdx+k]) * Zinv[k][j];
						sm += Zinv[i][k] * Z[k][j];
					}
					Kr[i][j] = sm; // * m_meanOmega;
				}
			}

			// Apply similarity transform. 

			//for (size_t i(0) ; i < nchem ; ++i) {
			//  for (size_t j(0) ; j < nchem ; ++j) {
			//    Kr[i][j] *= Z[i][nchem]/Z[j][nchem];
			//  }
			//}

			string rcm(string("Rate coefficient matrix:"));
			Kr.print(rcm, stest);

		}
		else {

			//
			// Non-conservative system.
			//

		}

		mesmerRates = Kr;

		return true;

	}

	int CollisionOperator::getSpeciesSequenceIndex(const std::string ref)
	{
		Reaction::molMapType::iterator spcitr;
		for (spcitr = m_SpeciesSequence.begin(); spcitr != m_SpeciesSequence.end(); ++spcitr)
		{
			if (ref == (spcitr->first)->getName())
				return spcitr->second;
		}
		cerr << "No molecule named " << ref << " is available in the reaction species.";
		return -1;
	}

	void CollisionOperator::locateSinks()
	{
		m_sinkRxns.clear();
		for (size_t i(0); i < m_pReactionManager->size(); ++i) {

			Reaction* pReaction = (*m_pReactionManager)[i];
			ReactionType reactionType = pReaction->getReactionType();

			bool Irreversible = (reactionType == IRREVERSIBLE_ISOMERIZATION || reactionType == IRREVERSIBLE_EXCHANGE ||
				reactionType == DISSOCIATION || reactionType == BIMOLECULAR_SINK);
			if (Irreversible && m_sinkRxns.find(pReaction) == m_sinkRxns.end()) {
				// Add an irreversible rxn to the map.
				Molecule* rctnt = pReaction->get_reactant();
				if (reactionType == IRREVERSIBLE_EXCHANGE) {
					m_sinkRxns[pReaction] = m_sources[rctnt];
				}
				else { // Irreversible exchange reaction.
					m_sinkRxns[pReaction] = m_isomers[rctnt];
				}
			}
		}

	}

	// Accessor to get specified eigenvalue.
	double CollisionOperator::getEigenvalue(size_t idEigenvalue) const {

		// Check id is sensible.
		if (idEigenvalue > m_eigenvalues.size()) {
			throw std::runtime_error("Eigenvalue ID greater than collision operator size.");
		}

		return -m_meanOmega*to_double(m_eigenvalues[m_eigenvalues.size() - idEigenvalue]);
	}

	// Calculate Yields.
	void CollisionOperator::calculateYields(YieldMap &yieldMap, double &time) const {

		//
		// Yields are calculated by integrating the term Sum_i ki pi.
		// This effectively involves integration of pi between 0 and
		// infinity. This integral leads to the expression Sum_i ki M^(-1)pi_0.
		// The inversion is effected by inverting the eigenvalue expression, 
		// i.e. M^(-1) = FV(eigenvalues)^(-1)V^TF^(-1).
		//
		if (m_sinkRxns.size() == 0) {
			// No Sinks so throw an error.
			throw std::runtime_error("No sinks defined, therefore no yields can be calculated.");
		}

		// Get initial distribution.
		size_t smsize = m_eigenvalues.size();
		vector<double> p_0(smsize, 0.0);
		if (!produceInitialPopulationVector(p_0)) {
			throw std::runtime_error("Calculation of initial conditions vector failed.");
		}

		vector<qd_real> wrk(smsize, 0.0);
		for (size_t j(0); j < smsize; ++j) {
			wrk[j] = p_0[j] / m_eqVector[j];
		}

		(*m_eigenvectors).Transpose();
		wrk *= (*m_eigenvectors);

		if (time > 0.0) {

			// Experimental time.

			for (size_t j(0); j < smsize; ++j) {
				wrk[j] *= (exp(m_meanOmega*m_eigenvalues[j] * time) - 1.0) / (m_meanOmega*m_eigenvalues[j]);
			}
		}
		else {

			// Infinite time limit.

			for (size_t j(0); j < smsize; ++j) {
				wrk[j] /= fabs(m_meanOmega*m_eigenvalues[j]);
			}
		}

		(*m_eigenvectors).Transpose();
		wrk *= (*m_eigenvectors);

		double sum(0.0);
		for (size_t j(0); j < smsize; ++j) {
			wrk[j] *= m_eqVector[j];
			sum += to_double(wrk[j]);
		}

		sinkMap::const_iterator sinkitr = m_sinkRxns.begin();
		for (; sinkitr != m_sinkRxns.end(); ++sinkitr) {

			// Locate the sink reaction.
			Reaction* sinkReaction = sinkitr->first;
			size_t rxnMatrixLoc = sinkitr->second;

			// Calculate the total flux through this channel.
			// First, determine the mirco rate coefficients for this channel:
			const vector<double> ktemp = sinkReaction->get_MtxGrnKf();  // Vector to hold sink k(E)s.
			// Now form the yield fraction. Note more than one channel may produce the same product.
			double yield(0.0);
			for (size_t i(0); i < ktemp.size(); ++i) {
				yield += to_double(ktemp[i] * wrk[rxnMatrixLoc + i]);
			}
			yieldMap[sinkReaction] = yield;
		}

	}

	// Calculate an experimental trace.
	void CollisionOperator::calculateTrace(const string &ref, vector<double> &times, vector<double> &signal) const {

		// Signal is assumed to be proportional to the total population of the monitored species.
		// The grain population is calculated in the usual way for the whole system and then the 
		// monitored species populations is calculated through integration. 

		// First we have to determine the location of the monitored species in the population vector.
		// As the species might be a source or isomer this requires iterating through the two maps.

		bool speciesFound(false);
		size_t lower(0), upper(0);

		// Iterate through the source map. 
		Reaction::molMapType::const_iterator spos;
		for (spos = m_sources.begin(); spos != m_sources.end(); ++spos) {
			Molecule* source = spos->first;
			if (source->getName() == ref) {
				speciesFound = true;
				lower = spos->second;
				upper = lower + 1;
			}
		}

		// Iterate through the isomer map.
		Reaction::molMapType::const_iterator ipos;
		for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos) {
			Molecule* isomer = ipos->first;
			if (isomer->getName() == ref) {
				speciesFound = true;
				lower = spos->second;
				upper = lower + isomer->getColl().get_colloptrsize();
			}
		}

		if (!speciesFound) {
			throw std::runtime_error("__FUNCTION__: No monitored species found amoung sources or isomers.");
		}

		// Get initial distribution.

		size_t smsize = m_eigenvalues.size();
		vector<double> p_0(smsize, 0.0);
		if (!produceInitialPopulationVector(p_0)) {
			throw std::runtime_error("__FUNCTION__: Calculation of initial conditions vector failed.");
		}

		// Save projected initial distribution

		vector<qd_real> c_0(smsize, 0.0);
		for (size_t j(0); j < smsize; ++j) {
			c_0[j] = p_0[j] / m_eqVector[j];
		}

		(*m_eigenvectors).Transpose();
		c_0 *= (*m_eigenvectors);
		(*m_eigenvectors).Transpose();

		// Calculate population at specific times.

		vector<qd_real> wrk(smsize, 0.0);
		for (size_t i(0); i < times.size(); ++i) {

			wrk = c_0;
			double time(times[i]);
			for (size_t j(0); j < smsize; ++j) {
				wrk[j] *= (exp(m_meanOmega*m_eigenvalues[j] * time));
			}
			wrk *= (*m_eigenvectors);

			signal[i] = 0.0;
			for (size_t j(lower); j < upper; ++j) {
				signal[i] += to_double(wrk[j] * m_eqVector[j]);
			}

		}

	}

	bool CollisionOperator::parseDataForGrainProfileAtTime(PersistPtr ppData)
	{
		//This is called from System::parse()
		//Grain Populations are now calculated at the same times for each species.
		//This means that m_GrainProfileAtTimeData is overcomplicated, but has
		//not been changed.
		PersistPtr pp = ppData, pp1;
		vector<Molecule*> refs;
		while ((pp1 = pp->XmlMoveTo("ref")) || (pp1 = pp->XmlMoveTo("me:ref")))
		{
			pp = pp1;
			const char* pRef = pp->XmlRead();
			Molecule* pMol = m_pMoleculeManager->find(pRef);
			if (!pMol)
				return false; //error message is in find()
			refs.push_back(pMol);
		}
		if (refs.empty())
		{
			cerr << " me:printGrainProfileAtTime needs one or more <me:ref> element to specify the species"
				<< endl;
			return false;
		}

		pp = ppData;
		double tim;
		vector<double> times;
		while (pp = pp->XmlMoveTo("me:time"))
		{
			const char* ptimtxt = pp->XmlRead();
			stringstream ss(ptimtxt);
			ss >> tim;
			times.push_back(tim);
		}
		if (times.empty())
		{
			cerr << "Need to specify at least one time in a <me:time> element in me:printGrainProfileAtTime";
			return false;
		}

		for (unsigned i = 0; i < refs.size(); ++i)
			m_GrainProfileAtTimeData.push_back(make_pair(refs[i], times));

		return true;
	}

	bool CollisionOperator::printGrainProfileAtTime(AnalysisData* analysisData) {

		// Check there is something to do.
		if (!m_GrainProfileAtTimeData.size())
			return true;

		// Use GrainProfileAtTimeData to calculate population
		// at each grain energy of each pMol at each time (Struan)

		size_t smsize = m_eigenvectors->size();
		vector<double> r_0(smsize, 0.); // initial distribution
		if (!projectedInitialDistrbtn(r_0)) {
			cerr << "Projection of initial disttribution failed.";
			return false;
		}

		// Copy full eigenvectors of the system.
		dMatrix totalEigenVecs(smsize);
		for (size_t i(0); i < smsize; ++i) {
			double tmp = to_double(m_eqVector[i]);
			for (size_t j(0); j < smsize; ++j) {
				totalEigenVecs[i][j] = tmp*to_double((*m_eigenvectors)[i][j]);
			}
		}

		// Iterate over species requested for output
		for (size_t iMol(0); iMol < m_GrainProfileAtTimeData.size(); ++iMol) {

			// Find the location of the species in the density vector.
			Molecule*  pMol = m_GrainProfileAtTimeData[iMol].first;
			int iLoc(-1);
			size_t slsize(0);
			if (m_isomers.find(pMol) != m_isomers.end()) {
				iLoc = m_isomers[pMol];
				slsize = pMol->getColl().get_colloptrsize();;
			}
			else if (m_sources.find(pMol) != m_sources.end()) {
				iLoc = m_sources[pMol];
				slsize = 1;
			}
			else {
				cerr << "Could not calculate species profile for " << pMol->getName() << "." << endl;
				continue;
			}

			const vector<double> Times(m_GrainProfileAtTimeData[iMol].second);
			vector<vector<double> > grnDists;

			if (analysisData)
				analysisData->m_grnTimes = Times;

			for (size_t iTime(0); iTime < Times.size(); ++iTime) {
				double numColl = m_meanOmega * Times[iTime];
				vector<double> p_t(smsize, 0.0);

				// |p_t> = exp(Lambda*t)*V^(T)*|init> = exp(Lambda*t)*U^(-1)*|n_0>
				for (size_t j(0); j < smsize; ++j) {
					p_t[j] = r_0[j] * exp(to_double(m_eigenvalues[j]) * numColl);
				}

				// |p_t> =  F*V*exp(Lambda*t)*V^(T)*|init> = U*exp(Lambda*t)*U^(-1)*|n_0> 

				p_t *= totalEigenVecs;

				// Copy densities for output.

				vector<double> density(p_t.begin() + iLoc, p_t.begin() + (iLoc + slsize - 1));
				grnDists.push_back(density);
			}
			if (analysisData)
				analysisData->m_grnDists[pMol->getName()] = grnDists;
		}

		return true;
	}

	bool CollisionOperator::projectedInitialDistrbtn(vector<double>& r_0) const {

		// This method calculates the projection of the initial distribution on to the
		// eigenspace of the collision matrix.

		vector<double> n_0 = r_0;
		if (!produceInitialPopulationVector(n_0)) {
			cerr << "Calculation of initial conditions vector failed.";
			return false;
		}

		// Convert the initial population vector into Boltzmann weighted population vector.
		// All transitions in the reaction matrix are Boltzmann weighted for symmetry.
		// |n_0> = F^(-1)*|n_0>
		for (size_t j(0); j < n_0.size(); ++j) {
			n_0[j] /= to_double(m_eqVector[j]);
		}

		// Multiply the initial population with the inverse of the eigenvector
		// which converts the populations into the "decay modes" domain.
		// |r_0> = V^(T)*F^(-1)*|n_0> = U^(-1)*|n_0>
		for (size_t i(0); i < r_0.size(); ++i) {
			double sum = 0.;
			for (size_t j(0); j < r_0.size(); ++j) {
				sum += n_0[j] * to_double((*m_eigenvectors)[j][i]);
			}
			r_0[i] = sum;
		}

		return true;
	}

}  //namespace
