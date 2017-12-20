#ifndef GUARD_Stochastic_h
#define GUARD_Stochastic_h

//-------------------------------------------------------------------------------------------
//
// Stochastic.h
//
// Author: Robin Shannon
// Date:   25/Jan/2016
//
// This header file contains the declaration of the stochastic class  
//
//-------------------------------------------------------------------------------------------
#include <numeric>
#include <set>
#include <random>
#include "CollisionOperator.h"


namespace mesmer
{
	class Stochastic
	{
	public:

		// Constructor
		Stochastic();

		// Destructor.
		virtual ~Stochastic();

		// Run the stochastic simulation
		bool simulate(CollisionOperator m_collisionOperator, const MesmerEnv &m_Env, const MesmerFlags &m_Flags);

	private:

		Molecule* getStartingSpecies(size_t& grndIdx);
		size_t getStartingIndex(Molecule* startSpecies);
		double rand();
		void addCollisionalConnections(vector<int>& connections, vector<Molecule*>& species, vector<pair<double, int>>& rates, Molecule *isomer, size_t grndIdx, size_t Idx, double Scale, bool grainBxd);
		void addReactiveConnections(vector<int>& connections, vector<Molecule*>& species, vector<pair<double, int>>& rates, Molecule *isomer, size_t ene, double scale, int MolIdx, size_t grnIdx, int TypeSwitch);
		void GetReactTypes(Molecule *isomer, bool& Pseudo, bool& Uni);
		bool timeAxis(vector<double>& timePoints, const MesmerEnv &m_Env);
		vector<double> correlation(vector<double> autoCorr);
		int get_speciesIdx(Molecule* specie, size_t& grndIdx);
		bool comparator(pair<double, int> l, pair<double, int> r) {return l.first < r.first;};

		int m_lowestBarrier;

		Molecule *Isomer;

		// Location of the molecule manager.
		MoleculeManager *m_pMoleculeManager;
		// Location of the reaction mananger.
		ReactionManager *m_pReactionManager;
		// Maps the location of individual reactant collision operator and source terms in the reaction operator.
		Reaction::molMapType    m_isomers;
		// Reaction Matrix used to track molecular trajectory
		qdMatrix               *m_reactionOperator;
		// Mean Collision Frequency
		double m_MeanOmega;
		// Mean Collision Frequency
		vector<qd_real>    m_eqVector;
		vector<qd_real>    m_eqVectorfull;
		//Matrix tracking whether isomers are equillibrated
		vector< vector <double>> transitionTable;
		vector< vector <bool>> equilibrated;
		vector <string> m_sinkNames;
		vector<double> m_auto;
		size_t m_sinkCount;
	};
}
#endif // GUARD_Stochastic_h