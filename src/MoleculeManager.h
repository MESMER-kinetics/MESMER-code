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
  const std::string                                         m_libfile;
private:
  void clear(void);

public:

  // Default constructor.
  MoleculeManager(const std::string& libraryfilename) 
    : m_molmap(), m_BathGasMolecule(), m_ppPersist(NULL), sourceNumber(0), m_libfile(libraryfilename) { } ;

  // Default destructor.
  ~MoleculeManager();

  // Add a new molecule to the list.
  Molecule*  addmol(string molName, string molType, PersistPtr ppMolList, const MesmerEnv& Env, MesmerFlags& Flags) ;

  // Find a molecule in the list.
  Molecule *find(const std::string& name) const ;

  // Remove a molecule from the list.
  void remove() {} ;

  // Total number of molecules in the list.
  void size() const {} ;

  //Return the Energy convention if all  molecules with _gDOS components have the same,
  //and an empty string otherwise
  std::string checkEnergyConventions();
  

/* I don't think these functions are needed within a flat Molecule
  Iterate through the molecules
   When id is empty, starts at beginning,
   When id contains the name of a molecule, the function recovers the
   next molecule of a class castable to that of pmol and returns true.
   Returns false when no molecules has been found.
   Use like:
     std::string id;
..   Molecule* pmol;
     while(GetNextMolecule(id, pmol))
     { //use id and pmol }
 */
/*
  bool GetNextMolecule(std::string& id, Molecule*& pmol)
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
      pmol = iter->second;
      if(pmol) //molecule of requested class found
      {
        id = iter->first;
        return true;
      }
    }
    return false;
  }
*/  
  /*Give the molecular name returns pointer*/
  bool GetThisMolecule(std::string& id, Molecule*& pmol)
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
        pmol = iter->second;
        return true;
      }
      else return false;
    }
    return false;
  }

  ///Search the library of molecules, copy to the main XML file,
  // replacing an existing molecule of same name. Return a pointer to the copy.
  PersistPtr GetFromLibrary(const std::string molName, PersistPtr ppMolList);

  // Accessors and Modifers for bath gas molecule.

  PersistPtr get_PersistPtr() {return m_ppPersist;}
  void set_PersistPtr(PersistPtr value) {m_ppPersist = value;}
  Molecule *get_BathGasMolecule() {return m_molmap[m_BathGasMolecule]; } ;
  void set_BathGasMolecule(const std::string &s_bgm){m_BathGasMolecule = s_bgm ; } ;

} ;
}//namespace
#endif // GUARD_MoleculeManager_h

