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

  std::map<std::string, Molecule*>                          m_molmap ;
  typedef std::map<std::string, Molecule*>::iterator        molIter ;
  typedef std::map<std::string, Molecule*>::const_iterator  constMolIter ;
  std::string                                               m_BathGasMolecule ;
  PersistPtr                                                m_ppPersist;
  int                                                       sourceNumber;
private:
   //pointers to molecules from library that have been successfully used
  std::vector<PersistPtr> m_libMols;
  void clear(void);

public:

  // Default constructor.
  MoleculeManager() : m_molmap(), m_BathGasMolecule(), m_ppPersist(NULL), sourceNumber(0) { } ;

  // Default destructor.
  ~MoleculeManager();

  // Add a new molecule to the list.
  Molecule*  addmol(string molName, string molType, PersistPtr ppMolList, const MesmerEnv& Env) ;

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
  
  /*Give the molecular name returns pointer*/
  template<typename T>
  bool GetThisMolecule(std::string& id, T*& pmol)
  {
    pmol=NULL;
    molIter iter;
    if(id.empty()){
      pmol = NULL; return false;
    }
    else
    {
      iter=m_molmap.find(id);
      if (iter != m_molmap.end()){
        pmol = dynamic_cast<T*>(iter->second);
        return true;
      }
      else return false;
    }
    return false;
  }

  ///Search the library of molecules, initialize pmolecule with the data
  ///Returns true if found
  bool LookinLibrary(const std::string molName, Molecule* pmolecule);
  vector<PersistPtr>& getLibraryMols(){return m_libMols;};

  // Accessors and Modifers for bath gas molecule.

  PersistPtr get_PersistPtr() {return m_ppPersist;}
  void set_PersistPtr(PersistPtr value) {m_ppPersist = value;}
  Molecule *get_BathGasMolecule() {return m_molmap[m_BathGasMolecule]; } ;
  void set_BathGasMolecule(const std::string &BathGasMolecule){m_BathGasMolecule = BathGasMolecule ; } ;

} ;
}//namespace
#endif // GUARD_MoleculeManager_h

