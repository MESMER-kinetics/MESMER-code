//-------------------------------------------------------------------------------------------
//
// MoleculeManager.cpp
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This file contains the implementation of the MoleculeManager class.
//
//-------------------------------------------------------------------------------------------

#include <iostream>
#include "MoleculeManager.h"

using namespace std ;
namespace mesmer
{
//
// Add a new molecule to the list.
//
bool MoleculeManager::addmols(PersistPtr ppMolList, const MesmerEnv& Env) {

  PersistPtr ppmol = ppMolList->XmlMoveTo("molecule");
  while(ppmol)
  {
    //-------------
    // Create a new Molecule of the required type.
    const char* ptype = ppmol->XmlReadValue("me:type", false);//some product mols don't have type
    /*need to think about this again, if a molecule has no type, it cannot store zpe, spin, etc...
    A mol has no type may probably for user to prepare some incomplete/unused molecules.
      Maybe I am wrong, but I thought that product molecules of irreversible reations need only a name.
      If some properties of a molecule are not going to be used then we should not make them
      madatory.
    */
    string moltype;
    if (ptype) {
      moltype = ptype;
    }
    else{
      moltype = "undefined";
    }

    Molecule *pmolecule;
    if     (moltype=="source")
      pmolecule = static_cast<Molecule*>(new SuperMolecule(Env));
    else if(moltype=="colliding")
      pmolecule = static_cast<Molecule*>(new CollidingMolecule(Env));
    else if(moltype=="modelled")
      pmolecule = static_cast<Molecule*>(new ModelledMolecule(Env));
    else if(moltype=="transitionState")
      pmolecule = static_cast<Molecule*>(new TransitionState(Env));
    else if(moltype=="bathGas")
      pmolecule = static_cast<Molecule*>(new BathGasMolecule(Env));
    else
      pmolecule = static_cast<Molecule*>(new Molecule(Env));

    //-------------
    // Initialize Molecule from input stream.
    // Each molecule type has its own set of mandatory parameters
    if(!pmolecule->InitializeMolecule(ppmol) || pmolecule->getName().empty()){
      stringstream errorMsg;
      errorMsg << "Failed in initializing the molecule. moltype = " << moltype;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      delete pmolecule;
      return false;
    } /*for the case of a SuperMolecule, if the element of source exists, the name has also to be specified.*/

    //strName has to go after InitializeMolecule()
    string strName = pmolecule->getName() ;

//    {stringstream errorMsg;
//    errorMsg << "Molecule " << strName << ", size of name = " << strName.size() << ", molecular type = " << moltype;
//    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);}

    //check if the molecule exists in m_molmap
    constMolIter it = m_molmap.find(strName) ;
    if (it != m_molmap.end()) {
      stringstream errorMsg;
      errorMsg << "Molecule " << strName << " already defined. Check input file to remove duplicates.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    else{
      //pmolecule->put_verbosity(true) ;
//      {stringstream errorMsg; 
//      errorMsg << "Adding Molecule " << strName << " into m_molmap, molecular type = " << moltype;
//      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);}

      //-------------
      // Add molecule to map.
      m_molmap[strName] = pmolecule ;
    }
    ppmol = ppmol->XmlMoveTo("molecule");
  }

  //check if there is a SuperMolecule (must at least have one)
  string superName = "source";
  constMolIter it = m_molmap.find(superName) ;
  if (it == m_molmap.end()) { // if there is no source, create an instance
    PersistPtr ppSuper = ppMolList->XmlWriteElement("molecule");
    ppSuper->XmlWriteAttribute("id", "source");
    ppSuper->XmlWriteAttribute("me:type", "source");
    Molecule *pmolecule = static_cast<Molecule*>(new SuperMolecule(Env));
    pmolecule->InitializeMolecule(ppSuper);
    m_molmap[superName] = pmolecule;
  }

  return true;
}

////
//// Add a new molecule to the list.
////
//SuperMolecule* MoleculeManager::addSuperMol(PersistPtr value) {
//  string myName = "source";
//  SuperMolecule *pmolecule = new SuperMolecule();
//  pmolecule->InitializeMolecule(value);
//  m_molmap[myName] = pmolecule;
//  return pmolecule;
//}

//
// Find a molecule in the list.
//
Molecule *MoleculeManager::find(const std::string& name) const {

  constMolIter it ;

  it = m_molmap.find(name) ;

  if (it == m_molmap.end()) {
    stringstream errorMsg;
    errorMsg << name << " is not a known Molecule.";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);

    return NULL;
  }

  return it->second ;
}

}//namespace
