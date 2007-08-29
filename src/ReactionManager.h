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
class System;//cannot include System.h because of circular references when compiling

class ReactionManager
{
public:

   // Type defs
   typedef  size_t  size_type ;

   ReactionManager(System* pSys, MoleculeManager *pMoleculeManager) : m_reactions(),
                                                        m_pMoleculeManager(pMoleculeManager), 
                                                        m_pSystemCollisionOperator(0),
                                                        m_pSys(pSys){} ;
   // Destructor.
   ~ReactionManager(){} ;

   System* GetSys() { return m_pSys; }

   // Add a new reaction to the map.
   bool addreactions(PersistPtr ReacList) ;

   // Remove a reaction from the map.
   void remove(){} ;

   // Total number of reaction in map.
   size_type size() const {return m_reactions.size() ; } ;

   // Find a particular reaction.
   Reaction*       operator[](const size_type i)       { return m_reactions[i] ; } ;
   const Reaction* operator[](const size_type i) const { return m_reactions[i] ; } ;

   // Find a reaction from its id
   Reaction* find(const std::string& id) const ;

  // Interrogates the (virtual) connectivity matrix, returning the reaction
  // index of the reaction (one based) connecting pProduct and pReactant
  // if both are CollidingMolecules
  // 0 if reactant and product are the same and are a CollidingMolecule
  // -1 otherwise
   int Connectivity(Molecule* pReactant, Molecule* pProduct);

   // Build collision operator for system.
   void BuildSystemCollisionOperator(const double beta, const double conc) ;

   // Diagonalize the collision operator.
   void diagCollisionOperator() ;

private:

   std::vector<Reaction *> m_reactions ;

   MoleculeManager        *m_pMoleculeManager ;

   dMatrix                *m_pSystemCollisionOperator ;

   System* m_pSys;

   // Default Constructor.
   ReactionManager() {} ;

} ;
}//namespace

#endif // GUARD_ReactionManager_h
