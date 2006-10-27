//-------------------------------------------------------------------------------------------
//
// ReactionManager.cpp 
//
// Author: Struan Robertson 
// Date:   11/Feb/2003
//
// This file contains the implementation of the ReactionManager class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"
#include "ReactionManager.h"
#include "Molecule.h"
#include "Constants.h"

using namespace Constants ;
namespace mesmer
{
// 
// Add a new reaction to the map.
//
bool ReactionManager::addreactions(TiXmlElement* pnReacList)
{
  TiXmlElement* pnReac = pnReacList->FirstChildElement("reaction");
  while(pnReac)
  {
     //
     // Create a new Reaction.
     //
     Reaction *preaction = new Reaction(m_pMoleculeManager) ;

     //
     // Initialize Reaction from input stream.
     //
     if(!preaction->ReadFromXML(pnReac))
       return false;;


//     preaction->put_verbosity(true) ;

     //
     // Add reaction to map.
     //
     m_reactions.push_back(preaction) ;

     pnReac = pnReac->NextSiblingElement("reaction");
  }
  return true;
}

int ReactionManager::Connectivity(Molecule* pReactant, Molecule* pProduct)
{
  return -1;
}

}//namespace
