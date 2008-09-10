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

        //check if the molecule exists in m_molmap
        constMolIter it = m_molmap.find(molName) ;
        if (it != m_molmap.end()) { // found a molecule with the same name --> should check its type as well later.
            return it->second;
        }

        //Need to construct a new Molecule of the required type
        Molecule *pmolecule;
        if(molType=="modelled")
            pmolecule = static_cast<Molecule*>(new CollidingMolecule(mEnv, mFlags));
        else if(molType=="deficientReactant")
            pmolecule = static_cast<Molecule*>(new ModelledMolecule(mEnv, mFlags));
        else if(molType=="excessReactant")
            pmolecule = static_cast<Molecule*>(new ModelledMolecule(mEnv, mFlags));
        else if(molType=="transitionState")
            pmolecule = static_cast<Molecule*>(new TransitionState(mEnv, mFlags));
        else if(molType=="bathGas")
            pmolecule = static_cast<Molecule*>(new BathGasMolecule(mEnv, mFlags));
        else if(molType=="sink")
            pmolecule = static_cast<Molecule*>(new ModelledMolecule(mEnv, mFlags));
        else
            pmolecule = static_cast<Molecule*>(new Molecule(mEnv, mFlags));

        //Look for it by name in the datafile
        PersistPtr ppmol = ppMolList;
        do{
            ppmol = ppmol->XmlMoveTo("molecule");
        }while(ppmol && ppmol->XmlReadValue("id")!= molName);

        // Initialize Molecule from input stream.
        // Each molecule type has its own set of mandatory parameters
        if(!ppmol || !pmolecule->InitializeMolecule(ppmol)){
            //If molecule with correct name not found, or if it does not initiallize properly, try the Library
            if(!LookinLibrary(molName, pmolecule)){
                //Still failed
                cerr << "Failed to find or initialize " << molName;
                delete pmolecule;
                return NULL;
            }
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

    bool MoleculeManager::LookinLibrary(const std::string molName, Molecule* pmolecule)
    {
        if(molName.empty())
            return false;
        static PersistPtr ppLib; //initiallized by system to NULL
        string libFilename("LibraryMols.xml");
        if(!ppLib)
        {
            const char* pdir = getenv("MESMER_LIBRARY");
            string libname;
            if(pdir)
            {
                libname = pdir;
                libname += FileSeparatorChar;
            }
            libname += libFilename;
            ppLib = XMLPersist::XmlLoad(libname);
            if(!ppLib)
            {
                cwarn << "Could not find Library file to search it for missing molecule(s)."<<endl;
                return false;
            }
        }
        cinfo << "Looked for " << molName << " in " << libFilename << endl;
        PersistPtr ppMolList   = ppLib->XmlMoveTo("moleculeList");
        if(!ppMolList)
            ppMolList = ppLib; //Can do without <moleculeList>
        PersistPtr ppMol = ppMolList->XmlMoveTo("molecule");
        while(ppMol)
        {
            if(molName==ppMol->XmlReadValue("id"))
            {
                //Check that the molecule can be initialised ok
                if(pmolecule->InitializeMolecule(ppMol))
                {
                    //If ok, write provenance and save in list so that it can be added to datafile later
                    //        ppMol->XmlWriteMainElement("me:source","Copied from " + libFilename);
                    ppMol->XmlWriteMetadata("source", libFilename);
                    m_libMols.push_back(ppMol);
                    return true;
                }
            }
            ppMol = ppMol->XmlMoveTo("molecule");
        }
        return false; // no suitable molecule found
    }

}//namespace
