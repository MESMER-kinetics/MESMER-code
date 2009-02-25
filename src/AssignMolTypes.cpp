#include "System.h"
using namespace std;

namespace mesmer
{
  // Deduce the me:type of each molecule and add it to the XML file 
  bool System::assignMolTypes(PersistPtr ppIOPtr) 
  {
    cerr << "Assigning mesmer molecule types..." << endl;

    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    
    PersistPtr ppReacList = ppIOPtr->XmlMoveTo("reactionList");
    if(!ppReacList)
    {
      cerr << "No reactantList." << endl;
      return false;
    }
    PersistPtr ppReac = ppReacList->XmlMoveTo("reaction");
    while(ppReac)//for all reactions
    {
      vector<PersistPtr> reactants, products;
      vector<PersistPtr>::iterator itr;
      
      PersistPtr pp = ppReac->XmlMoveTo("reactantList");
      if(!pp)
        pp=ppReac; //Be forgiving; we can get by without a reactantList element
      while(pp = pp->XmlMoveTo("reactant"))
        reactants.push_back(pp);

      pp = ppReac->XmlMoveTo("productList");
      if(!pp)
        pp=ppReac; //Be forgiving; we can get by without a productList element
      while(pp = pp->XmlMoveTo("product"))
        products.push_back(pp);

      pp = ppReac->XmlMoveTo("me:transitionState");
      if(pp)
      {
        pp = pp->XmlMoveTo("molecule");
        if(pp && !pp->XmlReadValue("me:type", false))//ignore if already present
          pp->XmlWriteAttribute("me:type", "transitionState");
      }

      //---------------------------------
      //With two reactants, the first radical is "deficientReactant" and the other species is "excessReactant" 
      //---------------------------------
      if(reactants.size()==2)
      {
        bool HasDeficientReactant =false;
        for(itr=reactants.begin(); itr!=reactants.end();++itr)
        {
          PersistPtr pprmol = (*itr)->XmlMoveTo("molecule");
          if(pprmol->XmlReadValue("me:type", false))
            continue;//ignore if has me:type already

          const char* pmolname;
          pmolname = pprmol->XmlReadValue("ref");
          if(pmolname)
          {
            PersistPtr ppMol = ppMolList;
            while(ppMol = ppMol->XmlMoveTo("molecule"))
            {
              const char* id = ppMol->XmlReadValue("id", false);
              if(id && !strcmp(id, pmolname))
                break;
            }
            if(!ppMol) //not found, try library. (Adds to datafile if found.)
              ppMol = m_pMoleculeManager->GetFromLibrary(pmolname, ppMolList);
            if(!ppMol) //if still not found cannot add types
              break;

            if(!HasDeficientReactant && ppMol->XmlReadPropertyDouble("me:spinMultiplicity", false)==2)
            {
              pprmol->XmlWriteAttribute("me:type", "deficientReactant"); //first radical
              HasDeficientReactant=true;
            }
            else
              pprmol->XmlWriteAttribute("me:type", "excessReactant");
          }
        }
      }

      //---------------------------------
      //With one reactant, it is "modelled"
      //---------------------------------
      else
      {
        PersistPtr ppmol = reactants[0]->XmlMoveTo("molecule");
        if(ppmol && !ppmol->XmlReadValue("me:type", false)) //ignore if has me:type already
          ppmol->XmlWriteAttribute("me:type", "modelled");
      }

      //---------------------------------
      //With two products, they are "sink" 
      //---------------------------------
      if(products.size()==2)
      {
        for(itr=products.begin(); itr!=products.end();++itr)
        {
          PersistPtr ppmol = (*itr)->XmlMoveTo("molecule");
          if(ppmol && !ppmol->XmlReadValue("me:type", false)) //ignore if has me:type already
            ppmol->XmlWriteAttribute("me:type", "sink");
        }
      }
      else
      //---------------------------------
      //With one product, it is "modelled" if the reaction is reversible and "sink" otherwise
      //---------------------------------
      {
        PersistPtr ppmol = products[0]->XmlMoveTo("molecule");
        if(ppmol && !ppmol->XmlReadValue("me:type", false)) //ignore if has me:type already
          ppmol->XmlWriteAttribute("me:type", ppReac->XmlReadBoolean("reversible") ? "modelled" : "sink");
      }

      ppReac = ppReac->XmlMoveTo("reaction");
    }
    cerr << "Check in the XML file that the assigned types are reasonable." << endl;
    return true;
  }
}//namespace mesmer

