#ifndef GUARD_MoleculeManager_h
#define GUARD_MoleculeManager_h

//-------------------------------------------------------------------------------------------
//
// MoleculeManager.h 
//
// Author: Struan Robertson 
// Date:   11/Feb/2003
//
// This header file contains the declaration of the MoleculeManager class. This class holds
// a collection of molecular specifications that may be part of a reaction. 
//
//-------------------------------------------------------------------------------------------

#include <string>
#include <map>
#include "Molecule.h"

namespace mesmer
{
class MoleculeManager {

public:

   // Default constructor.
   MoleculeManager() : m_molmap() { } ;

   // Default destructor.
   ~MoleculeManager() { } ;

   // Add a new molecule to the list.
   bool addmols(PersistPtr ppMolList) ;

   // Find a molecule in the list.
   Molecule *find(const std::string& name) const ;

   // Remove a molecule from the list.
   void remove() {} ;
  
   // Total number of molecules in the list.
   void size() const {} ;

private:

   std::map<std::string, Molecule*> m_molmap ;
 
} ;
}//namespace
#endif // GUARD_MoleculeManager_h

