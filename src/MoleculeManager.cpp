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
#include <stdlib.h>
#include "System.h"
#include "MoleculeManager.h"
#include "../gDensityOfStates.h"

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
  Molecule* MoleculeManager::addmol(string molName, string molType, const MesmerEnv& mEnv, MesmerFlags& mFlags) {

    ErrorContext c(molName);

    //check if the molecule exists in m_molmap
    constMolIter it = m_molmap.find(molName) ;
    if (it != m_molmap.end()) { 
	  // found a molecule with the same name --> should check its type as well later.
      // Check if the related properties specified by molType is activated to allow
	  // the molecule to play the specific role. If they are not activated, activate them.
      if (!((it->second)->activateRole(molType))){
        cerr << "Failed to initialize some molecular properties";
        return NULL;
      }
      return it->second;
    }

    //Construct a new Molecule
    Molecule *pmolecule = new Molecule(mEnv, mFlags, molType, this);

    //Look for it by name in the datafile
    PersistPtr ppmol = m_ppPersist;
    do{
      ppmol = ppmol->XmlMoveTo("molecule");
    }
    while(ppmol && ppmol->XmlReadValue("id")!= molName);

    // Initialize Molecule from input stream.
    while(!ppmol || !pmolecule->InitializeMolecule(ppmol)){
      //If molecule with correct name not found, or if it does not initiallize properly,
      // or is just a placeholder with no content, try the Library
      PersistPtr ppmol2 = GetFromLibrary(molName, m_ppPersist);
      if(!ppmol && !ppmol2){
        //No such molecule in datafile or library
        cerr << "Failed to find missing molecule " << molName << " in data file or library" << endl;
        delete pmolecule;
        return NULL;
      }
      if(!ppmol2)
        // Accept placeholder from datafile
        break;
      //Use the library version. Go back to initialize it.
      ppmol = ppmol2;
    }

    // Activate specified properties for the molecule
    if (!(pmolecule->activateRole(molType))){
      cerr << "Failed to initialize some molecular properties in " << molName;
      return NULL;
    }

    // Add molecule to map.
    m_molmap[molName] = pmolecule ;

    return pmolecule;
  }

  //
  // Find a molecule in the list.
  //
  Molecule *MoleculeManager::find(const std::string& name) const {

    constMolIter it ;

    it = m_molmap.find(name) ;

    if (it == m_molmap.end()) {
      cerr << name << " is not a known Molecule." << endl;
      return NULL;
    }

    return it->second ;
  }

  //Return the Energy convention if all  molecules with _gDOS components have the same,
  //and an empty string otherwise
  string MoleculeManager::checkEnergyConventions()
  {
    constMolIter iter;
    string firstConvention;
    for(iter=m_molmap.begin();iter!=m_molmap.end();++iter)
    {
      if((iter->second)->hasDOSProperties())//only "modelled" molecules
      {
        string thisConvention = iter->second->getDOS().getEnergyConvention();
        if(firstConvention.empty())
          firstConvention = thisConvention;
        if(firstConvention!=thisConvention){
          cerr << "cannot find an energy convention for molecule " << iter->second->getName() << endl;
          return string();
        }
      }
    }
    return firstConvention;
  }


}//namespace
