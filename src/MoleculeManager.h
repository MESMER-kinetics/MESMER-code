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

  typedef std::map<std::string, Molecule*>::const_iterator  constMolIter ;

  // Default constructor.
  MoleculeManager() 
    : m_molmap(), m_BathGasMolecule(), m_ppPersist(NULL), sourceNumber(0) { } ;

  // Default destructor.
  ~MoleculeManager();

  // Add a new molecule to the list.
  Molecule*  addmol(string molName, string molType, const MesmerEnv& Env, MesmerFlags& Flags) ;

  // Find a molecule in the list.
  Molecule *find(const std::string& name) const ;

  // Remove a molecule from the list.
  void remove() {} ;

  // Total number of molecules in the list.
  size_t size() const { return m_molmap.size() ; } ;

  constMolIter begin() const { return m_molmap.begin() ; } ;

  constMolIter end() const { return m_molmap.end() ; } ;

  void WriteEnergyConvention(); //to convention attribute on <moleculeList>

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

  // Accessors and Modifers for bath gas molecule.

  PersistPtr get_PersistPtr() {return m_ppPersist;}
  void set_PersistPtr(PersistPtr value) {m_ppPersist = value;}
  Molecule *get_BathGasMolecule() {return m_molmap[m_BathGasMolecule]; } ;
  void set_BathGasMolecule(const std::string &s_bgm){m_BathGasMolecule = s_bgm ; } ;
  const std::string& get_BathGasName() { return m_BathGasMolecule; } ;
private:

  std::map<std::string, Molecule*>                          m_molmap ;
  typedef std::map<std::string, Molecule*>::iterator        molIter ;
  std::string                                               m_BathGasMolecule ;
  PersistPtr                                                m_ppPersist;
  int                                                       sourceNumber;

  void clear(void);

} ;
}//namespace
#endif // GUARD_MoleculeManager_h

