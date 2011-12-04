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
    bool BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags, bool writeReport) ;

    // Diagonalize the reaction operator.
    void diagReactionOperator(const MesmerFlags &mFlags, const MesmerEnv &mEnv,
      Precision precision, PersistPtr ppAnalysis) ;

    // Calculate the time evolution of the system
    bool timeEvolution(MesmerFlags& mFlags,PersistPtr ppAnalysis, PersistPtr ppPopList);

    // Calculates the Bartis-Widom macroscopic rate coefficients.
    bool BartisWidomPhenomenologicalRates(qdMatrix& rates, MesmerFlags& mFlags, PersistPtr ppBase);

    // Calculates the Bartis-Widom macroscopic rate coefficients, using the contracted basis set eigenvectors.
    bool BartisWidomBasisSetRates(qdMatrix& rates, MesmerFlags& mFlags);

    // Write out phenomenological rate coefficients.
    bool PrintPhenomenologicalRates(qdMatrix& Kr, qdb2D& Kp, int numberOfCemeteries, MesmerFlags& mFlags, PersistPtr ppList) ;

    int getSpeciesSequenceIndex(const std::string ref);

    // Accessor to get specified eigenvalue.
    double getEigenvalue(size_t idEigenvalue) const ;

    // This method locates all sinks and determines their location in the relevant
    // isomer or source map. 
    void locateSinks() ;

    // Calculate Yields
    void calculateYields (YieldMap &yieldMap, double &time) const ;

    bool parseDataForGrainProfileAtTime(PersistPtr pp);

    bool printGrainProfileAtTime();

  private:

    // Sets grain parameters and determines system environment.
    bool SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne, bool writeReport);

    // Construct a transition matrix based on grains.
    void constructGrainMatrix(int msize);

    // Construct a transition matrix based on collision operator eigenfunctions.
    void constructBasisMatrix(void);

    void printReactionOperator(const MesmerFlags &mFlags);

    void printEigenvectors(const MesmerFlags &mFlags, std::ostream& os);

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

    bool                    m_punchSymbolGathered;

    // Species, times for printDataForGrainProfileAtTime
    std::vector<std::pair<Molecule*, std::vector<double> > > m_GrainProfileAtTimeData;

  } ;

}
#endif // GUARD_CollisionOperator_h
