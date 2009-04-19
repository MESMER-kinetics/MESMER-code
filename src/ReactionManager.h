#ifndef GUARD_ReactionManager_h
#define GUARD_ReactionManager_h

//-------------------------------------------------------------------------------------------
//
// ReactionManager.h
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This header file contains the declaration of the ReactionManager class.
// This class will contain the reactions that go to make up a system.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

namespace mesmer
{
  struct locationIdx{
    Molecule* mol;
    int fml; // full matrix location
    int rml; // reduced matrix location
    int fms; // full matrix size
    int rms; // reduced matrix size
  };

  struct divisionIdx{
    Molecule* mol;
    int fml; // full matrix location
    int rml; // reduced matrix location
    int ass; // active state size 
             // The first grain location of the active state, with respect to the full matrix location. 
             // asl == 0 if the whole well is active.)
    int fms; // full matrix size
  };

  class ReactionManager
  {
  public:

    // Type defs
    typedef std::map<Reaction* , int, Reaction::ReactionPtrLess> sinkMap ;

    ReactionManager(MoleculeManager *pMoleculeManager);

    // Destructor.
    ~ReactionManager()
    {
      if (m_reactionOperator) delete m_reactionOperator;
      vector<Reaction*>::iterator iter;
      for(iter=m_reactions.begin();iter!=m_reactions.end();++iter)
        delete *iter;
      m_reactions.clear();
    }

    // Add a new reaction to the map.
    bool addreactions(PersistPtr ReacList, const MesmerEnv& mEnv, MesmerFlags& mFlags) ;

    // Remove a reaction from the map.
    void remove(){} ;

    void resetCalcFlags();

    // Total number of reaction in map.
    size_t size() const {return m_reactions.size() ; } ;

    // Find a particular reaction.
    Reaction*       operator[](const size_t i)       { return m_reactions[i] ; } ;
    const Reaction* operator[](const size_t i) const { return m_reactions[i] ; } ;

    // Find a reaction from its id
    Reaction* find(const std::string& id) const ;

    // Build reaction operator for system.
    bool BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags) ;

    // Diagonalize the reaction operator.
    void diagReactionOperator(const MesmerFlags &mFlags, const int precision) ;

    // Calculate the time evolution of the system
    bool timeEvolution(MesmerFlags& mFlags);

    // Set Initial population for individual species
    void setInitialPopulation(PersistPtr);

    bool calculateEquilibriumFractions(const double beta);

    bool BartisWidomPhenomenologicalRates(dMatrix& rates, MesmerFlags& mFlags,PersistPtr ppBase);

    double calcChiSquare(const dMatrix& mesmerRates, vector<conditionSet>& expRates);

  private:

    // Construct a transition matrix based on grains.
    void constructGrainMatrix(int msize);
    
    // Construct a transition matrix based on collision operator eigenfunctions.
    void constructBasisMatrix(void);

    std::vector<Reaction *> m_reactions ;

    MoleculeManager        *m_pMoleculeManager ;

    qdMatrix               *m_reactionOperator ;
    qdMatrix               *m_eigenvectors;
    std::vector<qd_real>    m_eigenvalues;

    // Reaction operator after similarity transformation by block diagonal U, which is U^-1 M' U.
    qdMatrix               *m_basisMatrix;

    qdMatrix               *m_reducedBasisMatrix;
    qdMatrix               *m_reducedEigenvectors;
    std::vector<qd_real>    m_reducedEigenvalues;

    std::vector<double>     m_eqVector;
    std::vector<double>     m_reducedEqVector;
    std::vector<locationIdx> m_locSizeMap;
    std::vector<divisionIdx> m_divMap;

    // Maps the location of individual reactant collision operator and source terms in the reaction operator.
    Reaction::molMapType    m_isomers;
    Reaction::molMapType    m_sources;
    sinkMap                 m_sinkRxns;

    // map modelled molecules (isomers + sources) with their sequence in the EqMatrix and Rate Coefficient matrix
    Reaction::molMapType    m_SpeciesSequence;

    sinkMap m_SinkSequence;

    double m_meanOmega;

    bool punchSymbolGathered;

    // Extract molecule information from XML stream.
    bool GetMoleculeInfo(PersistPtr pp, std::string& MolName, std::string& MolType) ;

    // sets grain parameters and determines system environment
    bool SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne);

    bool produceInitialPopulationVector(vector<double>& initDist);

    bool produceEquilibriumVector();

    bool produceReducedEquilibriumVector();

    void printReactionOperator(const MesmerFlags &mFlags);

    void printEigenvectors(const MesmerFlags &mFlags);

  } ;
}//namespace

#endif // GUARD_ReactionManager_h
