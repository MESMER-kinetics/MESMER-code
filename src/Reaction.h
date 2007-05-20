#ifndef GUARD_Reaction_h
#define GUARD_Reaction_h

//-------------------------------------------------------------------------------------------
//
// Reaction.h 
//
// Author: Struan Robertson 
// Date:   1/Feb/2003
//
// This header file contains the declaration of the Reaction class.
//
//-------------------------------------------------------------------------------------------

#include <vector>
#include "Molecule.h"
#include "MoleculeManager.h"
#include "Persistence.h"

namespace mesmer
{
	
class Reaction
{
   //
   // Declare System as a friend class. Not sure that this is the best or 
   // most OOP way to go as it clearly defeats the point of encapsulation
   // but the System needs to know a lot about the reactions it is
   // combining. Review latter. SHR 2/Apr/2003.
   //

   friend class System ;

public:

   typedef enum ReactionType{ ASSOCIATION,
                              DISSOCIATION,
                              ISOMERIZATION,
                              ERROR_REACTION } ;

   // Constructors.
   Reaction(){} ;

   Reaction(MoleculeManager *pMoleculeManager);

   // Destructor.
   ~Reaction() ;

   // Copy constructor.
//   Reaction(const Reaction& reaction) ;

   // Assignment operator.
//   Reaction& operator=(const Reaction& reaction) ;

   // Initialize reaction.
   bool Initialize(PersistPtr ppReac) ;

   // Modifier for reaction type.
   void put_Reactiontype(ReactionType reactiontype) ;

   // Accessor for reaction type.
   ReactionType get_Reactiontype() const { } ;
  
   // Detail balance the microcanonical rate coefficients.
   void DetailedBalance() { } ;  

   // Determine the equilibrium constant.
   void CalcEquilConst() { } ;

   // Calculate the forward microcanoincal rate coeffcients. 
   void CalcMicroRateCoeffs() ;

   // Access microcanoincal rate coeffcients. 
   void get_MicroRateCoeffs(std::vector<double> &kmc, int ngrn) ;

   // Reactant information:

   int get_NumberOfReactants() const { return m_Reactant2 ? 2 : 1 ; } ;
   int get_NumberOfProducts()  const { return m_Product2 ? 2 : 1 ; } ;
   void get_unimolecularspecies(std::vector<CollidingMolecule *> &unimolecularspecies) const ;

private:
   ///Read a molecule name from the XML file and look it up
   Molecule* GetMolRef(PersistPtr pp);

   // Test the forward microcanoincal rate coeffcients. 
   void testMicroRateCoeffs() ;

   MoleculeManager   *m_pMoleculeManager ; // Pointer to molecule manager.

   //
   // Reaction composition.
   //
   CollidingMolecule *m_Reactant ;         // Reactant Molecule.
   Molecule          *m_Reactant2 ;        // Subsidiary reactant molecule.
   TransitionState   *m_TransitionState ;  // Transition State.
   CollidingMolecule *m_Product ;          // Product Molecule.
   Molecule          *m_Product2 ;         // Subsidiary product molecule.

   //
   // Reaction Rate data.
   //
   PersistPtr         m_ppPersist;         // Conduit for I/O
   ReactionType       m_reactiontype ;     // Type of reaction.
   double             m_E0 ;               // Reaction Threshold energy (measured from zero point of reactant).
   double             m_kfwd ;             // Forward canonical (high pressure) rate coefficient.
   double            *m_kfmc ;             // Forward microcanonical rate coefficients.

   //
   // Memory management.
   //
   std::allocator<double>        m_alloc ;

} ;
}//namespace
#endif // GUARD_Reaction_h
