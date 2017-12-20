//-------------------------------------------------------------------------------------------
//
// Stochastic.cpp
//
// Author: Robin Shannon
// Date:   26/Jan/2016
//
// This file contains implementation of the stochastic simulation class.
//
//-------------------------------------------------------------------------------------------
#include "Stochastic.h"
#include "gWellProperties.h"
#include "gPopulation.h"

using namespace std;
namespace mesmer
{
	//initialise variables i.e Number of trials, Max Time:
	Stochastic::Stochastic() {}

	Stochastic::~Stochastic() {}

	// Run Stochastic simulation.
	bool Stochastic::simulate(CollisionOperator m_collisionOperator, const MesmerEnv &m_Env, const MesmerFlags &m_Flags) {

		// Get a copy of the reaction operator, Molecule Manager and Reaction Manager?
		m_reactionOperator = m_collisionOperator.get_reactMat();
		m_pMoleculeManager = m_collisionOperator.get_molManager();
		m_pReactionManager = m_collisionOperator.get_ReacManager();
		m_isomers = m_collisionOperator.getIsomers();
		m_MeanOmega = m_collisionOperator.getMeanOmega();
		m_eqVectorfull = m_collisionOperator.getEqVec();
		m_sinkCount = m_collisionOperator.get_sinkSize();
		m_sinkNames = m_collisionOperator.getSinkNames();
		int molCount = m_isomers.size() + m_sinkCount;
		bool stochOnePass = false;

		// Get equillibrium vector for optional correction to thermalised systems
		Reaction::molMapType::iterator isomeritrb = m_isomers.begin();
		for (isomeritrb; isomeritrb != m_isomers.end(); ++isomeritrb) {
			m_eqVector.push_back(isomeritrb->first->getPop().getEqFraction());
		}

		for (size_t i(0); i < m_sinkCount; ++i) {
			m_eqVector.push_back(0.0);
		}

		//Get vector of times and make empty species profile to be populated
		vector<double> timeScale;
		timeScale.push_back(0.0);
		vector<double> Boltz;
		timeAxis(timeScale, m_Env);
		vector<double> autoData2(timeScale.size(), 0.0);

		vector<bool> firstPass(molCount, false);
		Molecule* pseudoIsomer = Isomer;
		int pseudoMolIdx = 0;
		int	pseudoGrnIdx = 0;
		double pseudoTime = 0.0;
		size_t pseudoGrndIdx = 0.0;
		size_t profileSize = 0.0;

		bool Reactive = false;
		bool Pseudo = false;
		bool Uni = false;

		//initialize vector which will hold the species time profiles
		vector< vector <double>> timeProfile(timeScale.size(), vector<double>(molCount));

		// Start Simulation
		for (int runs = 1; runs < m_Env.stochTrials + 1; ++runs)
		{
			//Set up constant parameters, may want to put all this in header
			int MolIdx = 0;
			vector< vector <double>> transitionTableNew(m_eqVector.size(), vector<double>(m_eqVector.size(), 0.0));
			vector< vector <bool>> equilibratedNew(m_eqVector.size(), vector<bool>(m_eqVector.size(), false));
			transitionTable = transitionTableNew;
			equilibrated = equilibratedNew;
			vector<double> tempAutoData2(timeScale.size(), 0.0);

			// Set the time at zero seccond
			double currentTime = 0;

			//initialize vector which will hold the species time profiles
			vector<double> timeScaleTemp;
			vector<vector<double>> timeProfileTemp;
			timeScaleTemp.push_back(0.0);

			//Get species with initial population of one and its corresponding index. 
			//Currently only one starting species possible.
			size_t grndIdx(0);
			Isomer = getStartingSpecies(grndIdx);
			Molecule* startIsomer = Isomer;
			MolIdx = get_speciesIdx(Isomer, grndIdx);
			map<int, double> grainMap;
			Isomer->getPop().getInitGrainPopulation(grainMap);
			int ene = 1;
			if (grainMap.size() >= 1) {
				map<int, double>::iterator grainIt = grainMap.begin();
				ene = grainIt->first;
			}
			else {
				//Add random ammount of energy based upon initial distribution
				int ene = getStartingIndex(Isomer);
			}

			// Start timeprofile with starting species
						//vector<double> StartEqVec(molCount, 0.0);
						//StartEqVec[MolIdx] = 1;
						//timeProfileTemp.push_back(StartEqVec);

			// Get energy grain
			size_t grnIdx = ene + grndIdx;

			// Run molecular trajectory until max timestep is reached
			while (currentTime < m_Env.stochEndTime)
			{

			end: vector<double> tempEqVec(molCount, 0.0);
				//Get all possible transition rates and indices and array in vectors
				//First get collision transfer "rates"
				vector<int> connections;
				vector<Molecule*> species;
				vector<pair<double, int>> rates(connections.size());




				// Scaling Factor depending on equilibration
				double scale = 1.0;

				// Stuff defining equilibration needs tidying
				if (equilibrated[MolIdx][MolIdx] == false) {
					firstPass[MolIdx] = false;
					tempEqVec[MolIdx] = 1;
					// Get rate coefficients for reactive processes
					if (!m_Flags.stochAutoCorrelation)
						addReactiveConnections(connections, species, rates, Isomer, ene, scale, MolIdx, grnIdx, 1);
					addCollisionalConnections(connections, species, rates, Isomer, grndIdx, grnIdx, scale, false);
				}
				//If species is thermalised, use bxd appraoch to promote reaction
				else if (!m_Flags.stochAutoCorrelation) {
					//check for species to species equilibration
					double eqTotal = 0.0;
					for (int j = 0; j < molCount; ++j) {
						tempEqVec[j] = equilibrated[MolIdx][j] == true ? to_double(m_eqVector[j]) : 0.0;
						eqTotal += equilibrated[MolIdx][j] == true ? to_double(m_eqVector[j]) : 0.0;
					}

					for (int j = 0; j < molCount; ++j) {
						tempEqVec[j] /= eqTotal;
					}

					//Get boltzman populations of various grains
					vector<double> boltzFrac;
					Isomer->getColl().normalizedInitialDistribution(boltzFrac);

					// Check types of reaction from this species
					Pseudo = false;
					Uni = false;
					GetReactTypes(Isomer, Pseudo, Uni);
					//Hack here  for MechGen Stuff
					Pseudo = true;

					// If PesudoIsomerisation then
					if (Pseudo == true && firstPass[MolIdx] == true) {

						size_t range = boltzFrac.size();

						for (size_t e2(0); e2 < boltzFrac.size(); ++e2) {
							scale = boltzFrac[e2];
							//Get rates
							addReactiveConnections(connections, species, rates, Isomer, e2, scale, MolIdx, grnIdx, 2);
						}


						if (Uni == false)
							firstPass[MolIdx] = false;
						if ((rates.size() <= 1))
							firstPass[MolIdx] = false;


					}

					if ((Uni == true && firstPass[MolIdx] == false) || Pseudo == false) {
						// Top and bottom will sum grain Boltzmann Factions above AXD threshold and below respectively
						double top = 0;
						firstPass[MolIdx] = false;
						double bottom = 0;
						m_lowestBarrier = (m_lowestBarrier - m_Env.stochAXDlimit) > 0 ? (m_lowestBarrier - m_Env.stochAXDlimit) : 0;
						size_t lower = m_lowestBarrier;
						if (ene < m_lowestBarrier)
						{
							ene = m_lowestBarrier;
							grnIdx = grndIdx + ene;
						}
						for (size_t e(0); e < lower; ++e) {
							bottom += boltzFrac[e];
						}
						size_t upper = boltzFrac.size();
						for (size_t e2(lower); e2 < upper; ++e2) {
							top += boltzFrac[e2];
						}

						//Get correction factor to transition probabilites based on ratio of Boltzman Fractions in the two AXD boxes
						scale = top / (top + bottom);
						scale *= tempEqVec[MolIdx];


						//Calculate transitions probabilites stopping collisional transfer to lowest grains
						addCollisionalConnections(connections, species, rates, Isomer, grndIdx, grnIdx, scale, true);
						addReactiveConnections(connections, species, rates, Isomer, ene, scale, MolIdx, grnIdx, 1);
					}
				}
				timeProfileTemp.push_back(tempEqVec);

				// sum rates
				double sum1 = 0;
				for (size_t i = 0; i < rates.size(); ++i) {
					sum1 += rates[i].first;
				}

				// sort rate pair vector
				std::sort(rates.begin(), rates.end());

				double random = rand();

				//Determine timesteps
				double nextSimulationTime = (-log(rand()) / sum1);
				double transitionTime = rand() * sum1;
				currentTime += nextSimulationTime;

				//then pick transition.       
				double sum = 0.0;

				//Get new grain Index
				size_t nwGrndIdx;
				int newGrnIdx;
				int speciesIdx = 0;
				bool first = firstPass[MolIdx];
				for (size_t i = 0; i < rates.size(); ++i) {
					sum += rates[i].first;
					if (sum > transitionTime) {
						newGrnIdx = connections[rates[i].second];
						// Sink Encountered
						if (newGrnIdx < 0) {
							MolIdx = newGrnIdx;
							if (stochOnePass) {
								cerr << "Ene 0" << endl;
								cerr << "Time " << currentTime << endl;
								size_t index = -MolIdx - m_isomers.size();
								cerr << "Product " << m_sinkNames[index] << endl;
								throw(std::runtime_error("one stocastic reaction performed"));
							}
							if ((currentTime >= pseudoTime) && (pseudoTime != 0)) {
								goto pseudoBit;
							}
							pseudoTime = 0;
							break;
						}

						speciesIdx = get_speciesIdx(species[rates[i].second], nwGrndIdx);
						if (first != true) {
							grndIdx = nwGrndIdx;
						}
						if (speciesIdx != MolIdx) {
							Reactive = true;
							if (stochOnePass) {
								cerr << "Ene " << (newGrnIdx - nwGrndIdx) << endl;
								cerr << "Time " << currentTime << endl;
								cerr << "Product " << ((species[rates[i].second])->getName()) << endl;
								throw(std::runtime_error("one stocastic reaction performed"));
							}
							if (equilibrated[speciesIdx][speciesIdx] == true) {
								firstPass[speciesIdx] = true;
							}
						}
						//Track transitions in table
						if (transitionTable[MolIdx][speciesIdx] != -2)
							transitionTable[MolIdx][speciesIdx] += 1;
						// Update connections table
						for (int j = 0; j < molCount; ++j) {
							//Check whether species to species equilibration threshold has been reached
							if ((transitionTable[MolIdx][j] >= m_Env.stochEquilThresh) && (speciesIdx != MolIdx)) {
								transitionTable[MolIdx][j] = -2;
								equilibrated[MolIdx][j] = true;
							}
							// Then check whether thermalisation threshold has been reached
							if ((transitionTable[MolIdx][j] >= m_Env.stochThermThresh) && (speciesIdx == MolIdx)) {
								transitionTable[MolIdx][j] = -2;
								equilibrated[MolIdx][j] = true;
								firstPass[MolIdx] = true;
							}
							//If a transition occurs rest counters for all other transitions 
							if ((j != speciesIdx)) {
								//unless a transtion is already in equillibrium
								transitionTable[MolIdx][j] = (transitionTable[MolIdx][j] != -2) ? 0 : -2;
							}
						}
						if (newGrnIdx > 0) {
							if (first == true) {
								pseudoIsomer = species[rates[i].second];
								pseudoMolIdx = speciesIdx;
								pseudoGrnIdx = newGrnIdx;
								pseudoTime = currentTime;
								pseudoGrndIdx = nwGrndIdx;
								currentTime -= nextSimulationTime;
								profileSize = timeProfileTemp.size();
								firstPass[MolIdx] = false;
								first = false;
								speciesIdx = MolIdx;
								Reactive = false;
								timeProfileTemp.pop_back();
								goto end;

							}
							else {
								Isomer = species[rates[i].second];
								MolIdx = speciesIdx;
							}
						}
						//break loop
						break;
					}
				}

				if ((Reactive == true || (currentTime > m_Env.stochEndTime)) && pseudoTime > 0) {
					if (currentTime >= pseudoTime) {
					pseudoBit:
						currentTime = pseudoTime;
						Isomer = pseudoIsomer;
						MolIdx = pseudoMolIdx;
						newGrnIdx = pseudoGrnIdx;
						grndIdx = pseudoGrndIdx;
						pseudoTime = 0;

						//remove unwanted elements from time profile vectors
						while (timeProfileTemp.size() >= profileSize) {
							timeProfileTemp.pop_back();
							timeScaleTemp.pop_back();
						}
					}
					else {
						pseudoTime = 0;
					}
				}

				Reactive = false;

				if (MolIdx < 0) {
					timeScaleTemp.push_back(currentTime);
					vector<double> tempEqVec2(molCount, 0.0);
					tempEqVec2[-MolIdx] = 1.0;
					timeProfileTemp.push_back(tempEqVec2);
					timeScaleTemp.push_back(100000000);
					break;
				}

				//assign a new grain index and update time vector
				grnIdx = newGrnIdx;
				ene = grnIdx - grndIdx;


				if (ene < 0)
					ene = 0;

				if (m_Flags.stochAutoCorrelation)
					m_auto.push_back(ene);
				timeScaleTemp.push_back(currentTime);



			}
			ctest << endl << "StochasticTrialComplete" << endl;

			// Create species profile for current stochastic trajectory
			if (m_Flags.stochAutoCorrelation)
			{
				autoData2 = m_auto;
				timeScale = timeScaleTemp;
			}
			else
			{
				size_t tp = 0;
				for (size_t te(1); te < timeScale.size(); ++te) {
					for (tp; tp < timeScaleTemp.size(); ++tp) {
						if (timeScale[te] > timeScaleTemp[tp] && timeScale[te] < timeScaleTemp[tp + 1]) {
							for (size_t miter(0); miter < molCount; ++miter) {
								timeProfile[te][miter] += timeProfileTemp[tp][miter];
							}
							break;
						}
					}
				}

			}
		}

		//Normalise time profile vector
		if (!m_Flags.stochAutoCorrelation)
		{
			for (size_t te(1); te < timeScale.size(); ++te) {
				double counts = 0;
				for (size_t isom(0); isom < molCount; ++isom) {
					counts += timeProfile[te][isom];
				}
				if (counts > 0) {
					for (size_t isom2(0); isom2 < molCount; ++isom2) {
						timeProfile[te][isom2] /= counts;
					}
				}
			}
		}
		if (!m_Flags.stochAutoCorrelation) {
			//Print Species Profile
			ctest << endl << "Print time dependent species and product profiles" << endl << "{" << endl;
			ctest << setw(16) << "Timestep/s";

			Reaction::molMapType::iterator jpos;
			for (jpos = m_isomers.begin(); jpos != m_isomers.end(); ++jpos) {
				Molecule* isomer = jpos->first;
				ctest << setw(16) << isomer->getName();
			}
			for (size_t i(0); i < m_sinkCount; ++i) {
				ctest << setw(16) << m_sinkNames[i];
			}

			ctest << endl;
			for (int i = 1; i < timeScale.size(); ++i) {
				ctest << timeScale[i] << "     ";
				for (int j = 0; j < timeProfile[i].size(); ++j) {
					ctest << timeProfile[i][j] << "      ";
				}
				ctest << endl;
			}
		}

		else {
			vector<double> finalCorrelation = correlation(autoData2);
			//Print AutoCorrelation function of internal energy
			ctest << endl << "AutoCorrelation Function" << endl << "{" << endl;

			for (size_t npos(0); npos < timeProfile.size(); npos++) {
				ctest << timeScale[npos] << "     ";
				ctest << setw(16) << finalCorrelation[npos];
				ctest << endl;
			}
			return true;
		}
	}
}


//Get the starting grain based upon initial populations and energy distributions
//Locate appropriate index in Reaction matrix
Molecule* Stochastic::getStartingSpecies(size_t& grndIdx) {
	//loop over isomers to find starting species. Sources not yet included for stochastic simulations
	Reaction::molMapType::iterator isomeritr = m_isomers.begin();
	for (isomeritr; isomeritr != m_isomers.end(); ++isomeritr) {
		Molecule* mol = isomeritr->first;
		size_t ind = isomeritr->second;
		if (mol->getPop().getInitPopulation() == 1.0) {
			grndIdx = ind;
			return mol;
		}
	}

}

//Get the starting grain based upon initial populations and energy distributions
//Locate appropriate index in Reaction matrix
size_t Stochastic::getStartingIndex(Molecule* isomer) {
	//Create vector of grain fractions and use a random number to select the appropriate one
	vector<double> boltzFrac;
	double r = rand();
	double sum = 0.0;
	size_t colloptrsize = isomer->getColl().get_colloptrsize();
	isomer->getColl().normalizedInitialDistribution(boltzFrac);
	for (size_t i(0); i < colloptrsize; ++i) {
		sum += boltzFrac[i];
		if (r <= sum)
			return i;
	}
}

//Add connections for collisional energy transfer to all wells of the same species
//void Stochastic::addCollisionalConnections(vector<int>& connections, vector<Molecule*>& species, vector<pair<double, int>>& rates, Molecule *isomer, size_t grndIdx, size_t Idx, double scale, bool grainBxd) {
//	size_t size = isomer->getColl().get_colloptrsize();
//	vector<vector <double>> Transition = isomer->getColl().get_rawTransition();
//	double col = isomer->getColl().get_collisionFrequency();
//	size_t grain = grainBxd ? (Idx - grndIdx) : 0;
//	if (grainBxd && (grain > m_lowestBarrier))
//		grain = m_lowestBarrier;
//	double sumC = 0;
//	for (grain; grain < size; ++grain) {
//		if (grain != (Idx - grndIdx)) {
//			size_t index = connections.size();
//			connections.push_back(grain + grndIdx);
//			species.push_back(isomer);
//			double element = to_double((*m_reactionOperator)[grain + grndIdx][Idx]);
//			double fraction = (to_double(m_eqVectorfull[Idx]) / to_double(m_eqVectorfull[grain + grndIdx]));
//			rates.push_back({ ((element / fraction) * m_MeanOmega) * scale, static_cast<int>(index) });
//			sumC += ((element / fraction) * m_MeanOmega) * scale;
//		}
//		}
//	//size_t index = connections.size();
//	//connections.push_back(Idx);
//	//species.push_back(isomer);
//	//double r = col - sumC;
// //   rates.push_back({ r, index });
//}

////Add connections for collisional energy transfer to all wells of the same species
void Stochastic::addCollisionalConnections(vector<int>& connections, vector<Molecule*>& species, vector<pair<double, int>>& rates, Molecule *isomer, size_t grndIdx, size_t Idx, double scale, bool grainBxd) {
	size_t size = isomer->getColl().get_colloptrsize();
	vector<vector <double>> Transition = isomer->getColl().get_rawTransition();
	double col = isomer->getColl().get_collisionFrequency();
	size_t grain = grainBxd ? (Idx - grndIdx) : 0;
	double sum = 0;
	if (grainBxd && (grain > m_lowestBarrier))
		grain = m_lowestBarrier;
	for (grain; grain < size; ++grain) {
		if (grain != Idx - grndIdx) {
			size_t index = connections.size();
			connections.push_back(grain + grndIdx);
			species.push_back(isomer);
			rates.push_back({ Transition[grain][Idx - grndIdx] * col * scale, static_cast<int>(index) });
			sum += Transition[grain][Idx - grndIdx];
		}
	}
}

////Add connections for collisional energy transfer to all wells of the same species
//void Stochastic::addCollisionalConnections(vector<int>& connections, vector<Molecule*>& species, vector<pair<double, int>>& rates, Molecule *isomer, size_t grndIdx, size_t Idx, double scale, bool grainBxd) {
//	size_t size = isomer->getColl().get_colloptrsize();
//	vector<vector <double>> Transition = isomer->getColl().get_rawTransition();
//	double collision = isomer->getColl().get_collisionFrequency();
//	size_t grain = grainBxd ? (Idx - grndIdx) : 0;
//	double N_y = 0.0;
//
//	for (grain; grain < size; ++grain) {
//		Norm[size] += Transition[size][grain];
//	}
//	Norm[size] *= grain_size;
//
//	for (size_t j(size - 1); j > 0; --grain) {
//		for (size_t i(0); i < j; ++i) {
//			N1[j] += Transition[j][i];
//		}
//		for (size_t k(j+1); k < size; ++k) {
//			N2[j] += Transition[k][j];
//		}
//		Norm[j] =
//	}
//}

//Find all possible reactive transitions for given species and add rate coefficients for these
// to a vector of possible transitions
void Stochastic::addReactiveConnections(vector<int>& connections, vector<Molecule*>& species, vector<pair<double, int>>& rates, Molecule *isomer, size_t ene, double scale, int MolIdx, size_t grnIdx, int TypeSwitch) {
	int barrier;
	m_lowestBarrier = 5000;
	for (size_t i(0); i < m_pReactionManager->size(); ++i) {
		vector<Molecule *> uni;
		size_t size = isomer->getColl().get_colloptrsize();
		(*m_pReactionManager)[i]->get_unimolecularspecies(uni);

		//search for reactions where isomer is the reactant and then get the product
		Reaction::molMapType::iterator isomeritr = m_isomers.begin();
		if (uni[0] == isomer) {
			vector<double> forTrans = (*m_pReactionManager)[i]->get_forwardTransition();
			size_t counter(0);
			//Searching through species to check which are equilibrated to current
			for (isomeritr; isomeritr != m_isomers.end(); ++isomeritr) {
				Molecule *product = isomeritr->first;
				bool isEquil = false;
				if ((equilibrated[MolIdx][counter] == true))
					isEquil = true;
				++counter;
				if ((*m_pReactionManager)[i]->getReactionType() != IRREVERSIBLE_ISOMERIZATION && (*m_pReactionManager)[i]->getReactionType() != DISSOCIATION && (*m_pReactionManager)[i]->getReactionType() != BIMOLECULAR_SINK) {
					//If a species is the product and the reaction is not equilibrated add transition rates to list
					if (uni.size() == 2 && product == uni[1] && !(isEquil)) {
						size_t sizeP = product->getColl().get_colloptrsize();
						int diff = static_cast<int>(size) - static_cast<int>(sizeP);
						if (((*m_pReactionManager)[i]->getReactionType() == PSEUDOISOMERIZATION) && (TypeSwitch != 3))
							barrier = (*m_pReactionManager)[i]->get_EffGrnRvsThreshold();
						else if ((*m_pReactionManager)[i]->getReactionType() != PSEUDOISOMERIZATION)
							barrier = (*m_pReactionManager)[i]->get_EffGrnFwdThreshold();
						else
							barrier = -1;
						if (barrier < m_lowestBarrier && barrier > 0)
							m_lowestBarrier = barrier;
						// if reaction is a pseudo-isomerisation then need a vector of rates not a single value
						if (((*m_pReactionManager)[i]->getReactionType() == PSEUDOISOMERIZATION) && (TypeSwitch != 3)) {
							if (ene > diff) {
								vector<vector<double>> forTransMat = (*m_pReactionManager)[i]->get_reverseTransitionMatrix();
								for (size_t jdx(0); jdx < (ene - diff); ++jdx) {
									size_t idx = isomeritr->second + jdx;
									size_t index = connections.size();
									connections.push_back(idx);
									species.push_back(product);
									double op = to_double((*m_reactionOperator)[idx][grnIdx]);
									double frac = (to_double(m_eqVectorfull[grnIdx]) / to_double(m_eqVectorfull[idx]));
									rates.push_back({ forTransMat[ene][jdx] * scale, index });
								}
							}
						}

						//Else we already have the rates in forTrans, Species index set to -1 in case of sink
						else if (TypeSwitch != 2 && (*m_pReactionManager)[i]->getReactionType() != PSEUDOISOMERIZATION) {
							size_t idx = isomeritr->second + ene - diff;
							if ((*m_pReactionManager)[i]->getReactionType() == IRREVERSIBLE_ISOMERIZATION) {
								size_t idx = -1;
							}
							size_t index = connections.size();
							vector<double> forTrans = (*m_pReactionManager)[i]->get_forwardTransition();

							if (static_cast<int>(ene - diff) > 0) {
								connections.push_back(idx);
								species.push_back(product);
								double op = to_double((*m_reactionOperator)[grnIdx][idx]);
								double op2 = to_double((*m_reactionOperator)[idx][grnIdx]);
								double frac = (to_double(m_eqVectorfull[idx]) / to_double(m_eqVectorfull[grnIdx]));
								double frac2 = (to_double(m_eqVectorfull[grnIdx]) / to_double(m_eqVectorfull[idx]));
								double sfrac = sqrt(to_double(m_eqVectorfull[idx]) / to_double(m_eqVectorfull[grnIdx]));
								double sfrac2 = sqrt(to_double(m_eqVectorfull[grnIdx]) / to_double(m_eqVectorfull[idx]));
								rates.push_back({ forTrans[ene] * scale, index });
							}
						}

					}
				}

				else if (counter == 1 && (TypeSwitch != 2)) {
					for (int s(0); s < m_sinkNames.size(); ++s) {
						vector<Molecule *> prod;
						(*m_pReactionManager)[i]->get_products(prod);
						string name = prod[0]->getName();
						//If a species is the product and the reaction is not equilibrated add transition rates to list
						if (name == m_sinkNames[s]) {
							int sinkPos = m_isomers.size() + s;
							barrier = (*m_pReactionManager)[i]->get_EffGrnFwdThreshold();
							if (barrier < m_lowestBarrier)
								m_lowestBarrier = barrier;
							int idx = -sinkPos;
							size_t index = connections.size();
							connections.push_back(idx);
							species.push_back(uni[0]);
							rates.push_back({ forTrans[ene] * scale, index });
						}
					}
				}
			}
		}

		//As above for reverse reactions as defined by the reactionManager
		else if (uni.size() == 2 && uni[1] == isomer) {
			vector<double> revTrans = (*m_pReactionManager)[i]->get_reverseTransition();
			size_t counter(0);
			for (isomeritr; isomeritr != m_isomers.end(); ++isomeritr) {
				Molecule *product = isomeritr->first;
				bool isEquil = false;
				if ((equilibrated[MolIdx][counter] == true))
					isEquil = true;
				++counter;
				if ((product == uni[0]) && !(isEquil)) {
					size_t sizeP = product->getColl().get_colloptrsize();
					int diff = static_cast<int>(size) - static_cast<int>(sizeP);
					if (((*m_pReactionManager)[i]->getReactionType() == PSEUDOISOMERIZATION) && (TypeSwitch != 3))
						barrier = (*m_pReactionManager)[i]->get_EffGrnFwdThreshold();
					else if ((*m_pReactionManager)[i]->getReactionType() != PSEUDOISOMERIZATION)
						barrier = (*m_pReactionManager)[i]->get_EffGrnRvsThreshold();
					else
						barrier = -1;
					if ((barrier < m_lowestBarrier || m_lowestBarrier == 0) && barrier > 0)
						m_lowestBarrier = barrier;
					if (((*m_pReactionManager)[i]->getReactionType() == PSEUDOISOMERIZATION) && (TypeSwitch != 3)) {
						if (i != 2) {
							for (size_t jdx(ene - diff); jdx < (sizeP - barrier); ++jdx) {
								vector<vector<double>> revTransMat = (*m_pReactionManager)[i]->get_forwardTransitionMatrix();
								size_t index = connections.size();
								size_t idx = isomeritr->second + jdx;
								connections.push_back(idx);
								species.push_back(product);
								double op = to_double((*m_reactionOperator)[idx][grnIdx]);
								double frac = (to_double(m_eqVectorfull[grnIdx]) / to_double(m_eqVectorfull[idx]));
								rates.push_back({ revTransMat[ene][jdx + diff] * scale, index });
							}
						}
					}
					else if (TypeSwitch != 2 && ((*m_pReactionManager)[i]->getReactionType() != PSEUDOISOMERIZATION)) {
						size_t idx = isomeritr->second + ene - diff;
						size_t index = connections.size();
						if (static_cast<int>(ene - diff) > 0) {
							connections.push_back(idx);
							species.push_back(product);
							double op = to_double((*m_reactionOperator)[grnIdx][idx]);
							double op2 = to_double((*m_reactionOperator)[idx][grnIdx]);
							double frac = (to_double(m_eqVectorfull[grnIdx]) / to_double(m_eqVectorfull[idx]));
							double frac2 = (to_double(m_eqVectorfull[idx]) / to_double(m_eqVectorfull[grnIdx]));
							rates.push_back({ revTrans[ene] * scale, index });;
						}
					}
				}
			}
		}
	}
}

//Find all possible reactive transitions for given species and add rate coefficients for these
// to a vector of possible transitions
void Stochastic::GetReactTypes(Molecule *isomer, bool& Pseudo, bool& Uni) {
	for (size_t i(0); i < m_pReactionManager->size(); ++i) {
		vector<Molecule *> uni;
		(*m_pReactionManager)[i]->get_unimolecularspecies(uni);
		//search for reactions where isomer is the reactant and then get the product
		Reaction::molMapType::iterator isomeritr = m_isomers.begin();
		if ((uni.size() > 1) && ((uni[1] == isomer) || (uni[0] == isomer))) {
			if ((*m_pReactionManager)[i]->getReactionType() == PSEUDOISOMERIZATION)
				Pseudo = true;
			else
				Uni = true;
		}
		if ((uni.size() == 1) && ((uni[0] == isomer))) {
			Uni = true;
		}
	}
}


//Random number generator
double Stochastic::rand() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	return dis(gen);

}

int Stochastic::get_speciesIdx(Molecule* specie, size_t & grndIdx) {
	Reaction::molMapType::iterator isomeritr = m_isomers.begin();
	int counter = 0;
	for (isomeritr; isomeritr != m_isomers.end(); ++isomeritr) {
		if (specie == isomeritr->first) {
			grndIdx = isomeritr->second;
			return counter;
		}

		++counter;
	}
}

bool Stochastic::timeAxis(vector<double>& timePoints, const MesmerEnv &m_Env) {
	// Calculates the time points
	for (int i = 0; i <= 300; ++i) {
		double thetime = pow(10., static_cast<double>(i) / 10. - 20.);
		if (thetime < m_Env.stochStartTime) continue;
		if (thetime > m_Env.stochEndTime) break;
		timePoints.push_back(thetime);
	}

	return true;
}

vector<double> Stochastic::correlation(vector<double> autoCorr) {
	vector<double> correlation(autoCorr.size(), 0.0);
	double totalAve = 0;
	for (size_t k(0); k < autoCorr.size(); k++) {
		totalAve += autoCorr[k];
	}
	totalAve /= autoCorr.size();
	for (size_t i(0); i < autoCorr.size(); i++) {
		vector<double> firstAve((autoCorr.size() - i), 0.0);
		double sum(0);
		for (size_t j(0); j < ((autoCorr.size()) - i); j++) {
			firstAve[j] = autoCorr[j] * autoCorr[j + i];
			sum += firstAve[j];
		}
		sum /= ((autoCorr.size()) - i);
		correlation[i] = sum - (totalAve*totalAve);
	}
	return correlation;
}


//Function for sorting a vector of pairs according to the first element
bool comparator(pair<double, int> l, pair<double, int> r)
{
	return l.first < r.first;
}


