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

#include "AssociationReaction.h"
#include "DissociationReaction.h"
#include "IsomerizationReaction.h"
#include "ExchangeReaction.h"
#include "marray.h"

namespace mesmer
{
  class ReactionManager
  {
  public:

    // Type defs
    typedef  size_t  size_type ;
    typedef std::map<std::string , double> populationMap ;

    ReactionManager(MoleculeManager *pMoleculeManager);

    // Destructor.
    ~ReactionManager(){} ;

    // Add a new reaction to the map.
    bool addreactions(PersistPtr ReacList, const MesmerEnv& Env) ;

    // Remove a reaction from the map.
    void remove(){} ;

    void resetCalcFlags();

    // Total number of reaction in map.
    size_type size() const {return m_reactions.size() ; } ;

    // Find a particular reaction.
    Reaction*       operator[](const size_type i)       { return m_reactions[i] ; } ;
    const Reaction* operator[](const size_type i) const { return m_reactions[i] ; } ;

    // Find a reaction from its id
    Reaction* find(const std::string& id) const ;

    // Build collision operator for system.
    bool BuildSystemCollisionOperator(MesmerEnv &Env) ;

    // Diagonalize the collision operator.
    void diagCollisionOperator(const MesmerEnv &Env) ;

    // Calculate the time evolution of the system
    bool timeEvolution(int maxTimeStep, const MesmerEnv mEnv);

    // Set Initial population for individual species
    void setInitialPopulation(PersistPtr);

    bool calculateEquilibriumFractions(const double beta);

  private:

    std::vector<Reaction *> m_reactions ;

    MoleculeManager        *m_pMoleculeManager ;

    dMatrix                *m_pSystemCollisionOperator ;

    std::vector<double>     m_eigenvalues;

    // Maps the location of individual reactant collision operator and source terms in the system matrix.
    Reaction::isomerMap    m_isomers;
    Reaction::sourceMap    m_sources;
    populationMap          m_populations;
    
    double m_meanOmega;

    // Default Constructor.
    //ReactionManager() {} ;

    // Extract molecule information from XML stream.
    bool GetMoleculeInfo(PersistPtr pp, std::string& MolName, std::string& MolType) ;

    // sets grain parameters and determines system environment
    bool SetGrainParams(MesmerEnv &Env, const double minEne, const double maxEne);

    bool produceInitialPopulationVector(vector<double>& eqFracCoeff, vector<double>& initDist);

  } ;
}//namespace

#endif // GUARD_ReactionManager_h
