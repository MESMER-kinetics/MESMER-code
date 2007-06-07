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

  /* Iterate through the molecules and return only those of a specified class.
   When id is empty, starts at beginning,
   When id contains the name of a molecule, the function recovers the
   next molecule of a class castable to that of pmol and returns true.
   Returns false when no molecules of the requested class has been found.
   Use like:
     std::string id;
..   ModelledMolecule* pmol;
     while(GetNextMolecule(id, pmol))
     { //use id and pmol }
 */
  template<typename T>
  bool GetNextMolecule(std::string& id, T*& pmol)
  {
    pmol=NULL;
    molIter iter;
    if(id.empty())
      iter=m_molmap.begin();
    else
    {
      //Set iter to the molecule after the one with name==id
      iter=m_molmap.find(id);
      if(iter!=m_molmap.end())
        ++iter;
    }
    
    for(;iter!=m_molmap.end();++iter)
    {
      pmol = dynamic_cast<T*>(iter->second);
      if(pmol) //molecule of requested class found
      {
        id = iter->first;
        return true;
      }
    }
    return false;
  }

   // Accessors and Modifers for bath gas molecule.

   Molecule *get_BathGasMolecule() {return m_molmap[m_BathGasMolecule]; } ;
   void set_BathGasMolecule(const std::string &BathGasMolecule){m_BathGasMolecule = BathGasMolecule ; } ;

private:

   std::map<std::string, Molecule*> m_molmap ;
   typedef std::map<std::string, Molecule*>::iterator molIter ;
   std::string m_BathGasMolecule ;
 
} ;
}//namespace
#endif // GUARD_MoleculeManager_h

