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
  
void MoleculeManager::clear(void){
  m_BathGasMolecule.clear();
  for (molIter i = m_molmap.begin(); i != m_molmap.end(); ++i) delete i->second;
  m_molmap.clear();
}

MoleculeManager::~MoleculeManager(){
  clear();
};


//
// Add a new molecule to the list.
//
Molecule* MoleculeManager::addmol(string& molName, string& molType, PersistPtr ppMolList, const MesmerEnv& Env) {

  Molecule *pmolecule;
  PersistPtr ppmol = ppMolList->XmlMoveTo("molecule");
  while (ppmol){
    string id = ppmol->XmlReadValue("id");
    if (id == molName) break; // found the molecule
    ppmol = ppmol->XmlMoveTo("molecule");
  }

  //check if the molecule exists in m_molmap
  constMolIter it = m_molmap.find(molName) ;
  if (it != m_molmap.end()) { // found a molecule with the same name --> should check its type as well later.
    return it->second;
  }
  else{
    //-------------
    // Create a new Molecule of the required type.
    if(molType=="modelled")
      pmolecule = static_cast<Molecule*>(new CollidingMolecule(Env));
    else if(molType=="reactant")
      pmolecule = static_cast<Molecule*>(new CollidingMolecule(Env));
    else if(molType=="excessReactant")
      pmolecule = static_cast<Molecule*>(new ModelledMolecule(Env));
    else if(molType=="transitionState")
      pmolecule = static_cast<Molecule*>(new TransitionState(Env));
    else if(molType=="bathGas")
      pmolecule = static_cast<Molecule*>(new BathGasMolecule(Env));
    else
      pmolecule = static_cast<Molecule*>(new Molecule(Env));

    //-------------
    // Initialize Molecule from input stream.
    // Each molecule type has its own set of mandatory parameters
    if(!pmolecule->InitializeMolecule(ppmol) || pmolecule->getName().empty()){
      stringstream errorMsg;
      errorMsg << "Failed in initializing the molecule. molType = " << molType;
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      delete pmolecule; return NULL;
    } /*for the case of a SuperMolecule, if the element of source exists, the name has also to be specified.*/

    //-------------
    // Add molecule to map.
    m_molmap[molName] = pmolecule ;
  }

//    {stringstream errorMsg;
//    errorMsg << "Molecule " << molName << ", size of name = " << molName.size() << ", molecular type = " << molType;
//    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);}


  if (molType == "reactant"){
    stringstream superId; superId << "source_" << sourceNumber; ++sourceNumber;
    cinfo << "source name = " << superId.str() << endl;
    PersistPtr ppSuper = NULL;
    //find if this source term is there
    PersistPtr ppMol = ppMolList->XmlMoveTo("molecule");
    while (ppMol){
      string myId = ppMol->XmlReadValue("id", false);
      if (myId == superId.str()){
        ppSuper = ppMol; break;
      }
      ppMol = ppMol->XmlMoveTo("molecule");
    }
    if (!ppSuper){
      ppSuper = ppMolList->XmlWriteElement("molecule");
      ppSuper->XmlWriteAttribute("id", superId.str());
      if (!ppSuper) {stringstream errorMsg;
        errorMsg << "Cannot get a persistent pointer.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }
    }
    SuperMolecule *pSMolecule = new SuperMolecule(Env);
    pSMolecule->InitializeMolecule(ppSuper);
    Molecule *pmolecule = static_cast<Molecule*>(pSMolecule);
    m_molmap[superId.str()] = pmolecule;
  }
  return pmolecule;
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
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);

    return NULL;
  }

  return it->second ;
}

}//namespace
