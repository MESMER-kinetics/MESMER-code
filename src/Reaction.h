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

using namespace std;
namespace mesmer
{
	
class Reaction : public IPersistObject
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
   bool ReadFromXML(TiXmlElement* pnReac) ;

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
   void get_MicroRateCoeffs(vector<double> &kmc, int ngrn) ;

   // Reactant information:

   std::string get_ReactantName() const { return m_Reactant->getName() ; } ;
   int get_NumberOfReactants()    const { return m_Reactant2 ? 2 : 1 ; } ;
   std::string get_ProductName()  const { return m_Product->getName() ; } ;
   int get_NumberOfProducts()     const { return m_Product2 ? 2 : 1 ; } ;

private:
   //Parse XML file for a product or a reactant and assigne to member variables
   bool GetMolRef(TiXmlElement* pnReac,const char* molRole, int index);

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
   TiXmlElement      *m_pXmlEl ;           // Its XML node address
   ReactionType       m_reactiontype ;     // Type of reaction.
   double             m_E0 ;               // Reaction Threshold energy (measured from zero point of reactant).
   double             m_kfwd ;             // Forward canonical (high pressure) rate coefficient.
   double            *m_kfmc ;             // Forward microcanonical rate coefficients.

   //
   // Memory management.
   //
   allocator<double>        m_alloc ;

} ;
}//namespace
#endif // GUARD_Reaction_h
