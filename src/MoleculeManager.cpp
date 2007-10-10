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
bool MoleculeManager::addmols(PersistPtr ppMolList) {

  PersistPtr ppmol = ppMolList->XmlMoveTo("molecule");
  while(ppmol)
  {
    //
    // Create a new Molecule of the required type.
    //
    const char* ptype = ppmol->XmlReadValue("me:type", false);//some product mols don't have type
    string moltype;
    if (ptype) moltype = ptype;

    Molecule *pmolecule;
    if(moltype=="modelled")
      pmolecule = static_cast<Molecule*>(new CollidingMolecule());
    else if(moltype=="transitionState")
      pmolecule = static_cast<Molecule*>(new TransitionState);
    else if(moltype=="bathGas")
      pmolecule = static_cast<Molecule*>(new BathGasMolecule);
    else
      pmolecule = static_cast<Molecule*>(new Molecule);

    //
    // Initialize Molecule from input stream.
    //Each molecule type has its own set of mandatory of parameters
    if(!pmolecule->InitializeMolecule(ppmol)){
    	delete pmolecule;
      return false;
    }

    string strName = pmolecule->getName() ;

    //pmolecule->put_verbosity(true) ;

    //
    // Add molecule to map.
    //
    m_molmap[strName] = pmolecule ;

    ppmol = ppmol->XmlMoveTo("molecule");
  }
  return true;
}

//
// Find a molecule in the list.
//
Molecule *MoleculeManager::find(const std::string& name) const { 

  map<string, Molecule*>::const_iterator it ;

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
