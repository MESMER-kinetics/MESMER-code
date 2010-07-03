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

    ReactionManager(MoleculeManager *pMoleculeManager);

    // Destructor.
    ~ReactionManager()
    {
      vector<Reaction*>::iterator iter;
      for(iter=m_reactions.begin();iter!=m_reactions.end();++iter)
        delete *iter;
      m_reactions.clear();
    }

    // Add a new reaction to the map.
    bool addreactions(PersistPtr ReacList, const MesmerEnv& mEnv, MesmerFlags& mFlags) ;

    // Remove a reaction from the map.
    void remove(){} ;

    // Total number of reaction in map.
    size_t size() const {return m_reactions.size() ; } ;

    // Find a particular reaction.
    Reaction*       operator[](const size_t i)       { return m_reactions[i] ; } ;
    const Reaction* operator[](const size_t i) const { return m_reactions[i] ; } ;

    // Find a reaction from its id
    Reaction* find(const std::string& id) const ;

    // Set Initial population for individual species
    void setInitialPopulation(PersistPtr);

  private:

    std::vector<Reaction *> m_reactions ;

    MoleculeManager        *m_pMoleculeManager ;

    // Extract molecule information from XML stream.
    bool GetMoleculeInfo(PersistPtr pp, std::string& MolName, std::string& MolType) ;

  } ;
}//namespace

#endif // GUARD_ReactionManager_h
