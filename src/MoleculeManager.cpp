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
  Molecule* MoleculeManager::addmol(string molName, string molType, PersistPtr ppMolList, const MesmerEnv& mEnv, MesmerFlags& mFlags) {

    ErrorContext c(molName);

    //check if the molecule exists in m_molmap
    constMolIter it = m_molmap.find(molName) ;
    if (it != m_molmap.end()) { // found a molecule with the same name --> should check its type as well later.
      // Check if the related properties specified by molType is activated to allow the molecule to play the specific role. 
      // If they are not activated, activate them.
        if (!((it->second)->activateRole(molType))){
          cerr << "Failed to initialize some molecular properties";
          return NULL;
        }
      return it->second;
    }

    //Construct a new Molecule
    Molecule *pmolecule = new Molecule(mEnv, mFlags);


    //Look for it by name in the datafile
    PersistPtr ppmol = ppMolList;
    do{
      ppmol = ppmol->XmlMoveTo("molecule");
    }
    while(ppmol && ppmol->XmlReadValue("id")!= molName);

    // Initialize Molecule from input stream.
    while(!ppmol || !pmolecule->InitializeMolecule(ppmol)){
      //If molecule with correct name not found, or if it does not initiallize properly,
      // or is just a placeholder with no content, try the Library
      PersistPtr ppmol2 = GetFromLibrary(molName, ppMolList);
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
      cerr << name << " is not a known Molecule.";
      return NULL;
    }

    return it->second ;
  }

  PersistPtr MoleculeManager::GetFromLibrary(const std::string molName, PersistPtr ppMolList)
  {
    //Search the library of molecules, copy to the main XML file and return a pointer to the copy.
    PersistPtr ppNewMol;
    if(molName.empty())
      return ppNewMol;
    static PersistPtr ppLib; //initiallized by system to NULL
    if(!ppLib)
    {
      ppLib = XMLPersist::XmlLoad(m_libfile,"");
      if(m_libfile.empty() || !ppLib)
      {
        cwarn << "Could not find Library file to search it for missing molecule(s)."<<endl;
        return false;
      }
    }
    PersistPtr ppLibMolList   = ppLib->XmlMoveTo("moleculeList");
    if(!ppLibMolList)
      ppLibMolList = ppLib; //Can do without <moleculeList>
    PersistPtr ppMol = ppLibMolList->XmlMoveTo("molecule");
    while(ppMol)
    {
      if(molName==ppMol->XmlReadValue("id", false))
      {
        //Delete a molecule of same name in datafile, if present
        PersistPtr ppOldMol = ppMolList;
        while(ppOldMol = ppOldMol->XmlMoveTo("molecule")) {
          if(molName == ppOldMol->XmlReadValue("id", false))
            break;
        }
        //Copy a matching molecule from library to the main XML file
        //Replace old version if present
        ppMol = ppMolList->XmlCopy(ppMol, ppOldMol);
        cinfo << molName << " copied from " << m_libfile << endl;
        //Write its provenance
        ppMol->XmlWriteMetadata("source", m_libfile);
        return ppMol;
      }
      ppMol = ppMol->XmlMoveTo("molecule");
    }
    cinfo << "Could not find " << molName << " in " << m_libfile << endl;
    return ppMol; // empty, no suitable molecule found
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
        if(firstConvention!=thisConvention)
          return string();
      }
    }
    return firstConvention;
  }


}//namespace
