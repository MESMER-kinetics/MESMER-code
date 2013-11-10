#ifndef GUARD_CollisionOperator_h
#define GUARD_CollisionOperator_h

//-------------------------------------------------------------------------------------------
//
// CollisionOperator.h
//
// Author: Struan Robertson
// Date:   26/Feb/2011
//
// This header file contains the declaration of the CollisionOperator class.
// This class will implement the master equation collision operator.
//
//-------------------------------------------------------------------------------------------

#include "dMatrix.h"
#include "ReactionManager.h"

namespace mesmer
{

  typedef std::map<Reaction* , double, Reaction::ReactionPtrLess> YieldMap ;

  class CollisionOperator
  {
  public:

    // Constructor
    CollisionOperator() ;

    // Destructor.
    virtual ~CollisionOperator() ;

    // Initialize the collision operator object.
    bool initialize(MoleculeManager *pMoleculeManager, ReactionManager *pReactionManager) ;

    // Calculate the equilibrium fraction of each species in the system.
    bool calculateEquilibriumFractions() ;

    // Build reaction operator for system.
    bool BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags, bool writeReport = false) ;

    // Diagonalize the reaction operator.
    void diagReactionOperator(const MesmerFlags &mFlags, const MesmerEnv &mEnv,
      Precision precision, PersistPtr ppAnalysis = NULL) ;

    // Calculate the time evolution of the system
    bool timeEvolution(MesmerFlags& mFlags,PersistPtr ppAnalysis, PersistPtr ppPopList, PersistPtr ppAvEList);

    // Calculates the Bartis-Widom macroscopic rate coefficients.
    bool BartisWidomPhenomenologicalRates(qdMatrix& rates, qdMatrix& lossRates, MesmerFlags& mFlags, PersistPtr ppBase = NULL);

    // Calculates the Bartis-Widom macroscopic rate coefficients, using the contracted basis set eigenvectors.
    bool BartisWidomBasisSetRates(qdMatrix& rates, MesmerFlags& mFlags);

    // Write out phenomenological rate coefficients.
    bool PrintPhenomenologicalRates(qdMatrix& Kr, qdMatrix& Kp, MesmerFlags& mFlags, PersistPtr ppList) ;

    int getSpeciesSequenceIndex(const std::string ref);

    // Accessor to get specified eigenvalue.
    double getEigenvalue(size_t idEigenvalue) const ;

    // Calculate Yields
    void calculateYields (YieldMap &yieldMap, double &time) const ;

    bool parseDataForGrainProfileAtTime(PersistPtr pp);

    bool printGrainProfileAtTime(PersistPtr ppGrainList);

    bool printAverageEnergies(PersistPtr ppAvList);

    bool hasGrainProfileData() { return !m_GrainProfileAtTimeData.empty(); }

	// Accessor for phenomenological rates.
	void get_phenRates(std::map<std::string, double> &phenRates) const {phenRates = m_phenomenlogicalRates ; } ;

  private:

    // Sets grain parameters and determines system environment.
    bool SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne, bool writeReport);

    // Set up max energy above the top hill automatically
    void SetMaximumCellEnergy(MesmerEnv &mEnv, const MesmerFlags& mFlags);

    // Construct a transition matrix based on grains.
    void constructGrainMatrix(int msize);

    // Construct a transition matrix based on collision operator eigenfunctions.
    void constructBasisMatrix(void);

    // Locate all sinks in the relevant isomer or source map. 
    void locateSinks() ;

    void printReactionOperator(const MesmerFlags &mFlags);

    bool produceEquilibriumVector();

    bool produceInitialPopulationVector(vector<double>& initDist) const ;

    bool projectedInitialDistrbtn(vector<double>& initDist) const ;

    // Location of the molecule manager.
    MoleculeManager *m_pMoleculeManager;

    // Location of the reaction mananger.
    ReactionManager *m_pReactionManager ;

    // Maps the location of individual reactant collision operator and source terms in the reaction operator.
    Reaction::molMapType    m_isomers;
    Reaction::molMapType    m_sources;

    typedef std::map<Reaction* , int, Reaction::ReactionPtrLess> sinkMap ;

    sinkMap                 m_sinkRxns;

    // Mean collision frequency.
    double                  m_meanOmega;

    // The system transition matrix and associated eigenvalues and eigenvectors.
    qdMatrix               *m_reactionOperator ;
    qdMatrix               *m_eigenvectors;
    std::vector<qd_real>    m_eigenvalues;

    // Map modelled molecules (isomers + sources) with their sequence in the transition matrix.
    Reaction::molMapType    m_SpeciesSequence ;

    // Equilibrium distribution.
    std::vector<qd_real>    m_eqVector;
    size_t                  m_eqVecSize ;

    bool                    m_punchSymbolGathered;

    // Species, times for printDataForGrainProfileAtTime
    std::vector<std::pair<Molecule*, std::vector<double> > > m_GrainProfileAtTimeData;

    // Map relating reactions with phenomenological rate coefficients.
    std::map<std::string, double> m_phenomenlogicalRates ;
  } ;

}
#endif // GUARD_CollisionOperator_h
