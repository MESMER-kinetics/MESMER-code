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
#include "ReactionManager.h"

#include "AssociationReaction.h"
#include "IrreversibleUnimolecularReaction.h"
#include "IsomerizationReaction.h"
#include "IrreversibleExchangeReaction.h"
#include "BimolecularSinkReaction.h"
#include "PseudoIsomerizationReaction.h"

using namespace Constants ;
using namespace std ;

namespace mesmer
{
  ReactionManager::ReactionManager(MoleculeManager *pMoleculeManager)
    :m_reactions(),
    m_pMoleculeManager(pMoleculeManager)
  {};

  //
  // Add a new reaction to the map.
  //
  bool ReactionManager::addreactions(PersistPtr ppReacList, const MesmerEnv& mEnv, MesmerFlags& mFlags)
  {
    bool readStatus(true), allReactionsOK(true);
    PersistPtr ppReac = ppReacList;
    while(ppReac = ppReac->XmlMoveTo("reaction"))
    {
      if(!readStatus)
        allReactionsOK = false; //parsing of previous reaction failed, but keep parsing
      readStatus =true; //for current reaction

      //ignore reactions with attribute active="false"
      const char* active = ppReac->XmlReadValue("active", optional);
      if(active && !strcmp(active, "false"))
        continue;

      //Read reaction ID
      const char* id = ppReac->XmlReadValue("id", false);
      if(!id){
        cinfo << "Reaction ID not found.\n";
        return false;
      }
      ErrorContext c(id);
      cinfo << "Parsing reaction..." << endl;
      cinfo.flush();

      // Read reactant and product types.

      string rct1Name, rct1Type, rct2Name, rct2Type ;
      string pdt1Name, pdt1Type, pdt2Name, pdt2Type ;
      bool bRct2(false), bPdt1(false), bPdt2(false) ;

      PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
      if(!ppReactantList)
        ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

      PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
      if(ppReactant1) {
        readStatus = readStatus && GetMoleculeInfo(ppReactant1, rct1Name, rct1Type) ;

        PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        if(ppReactant2) {
          readStatus = readStatus && (readStatus && GetMoleculeInfo(ppReactant2, rct2Name, rct2Type)) ;
          bRct2 = true;
        }
      }

      PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
      if(!ppProductList)
        ppProductList=ppReac; //Be forgiving; we can get by without a productList element

      PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
      if (ppProduct1) {
        readStatus = (readStatus && GetMoleculeInfo(ppProduct1, pdt1Name, pdt1Type)) ;
        bPdt1 = true ;

        PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
        if (ppProduct2){
          readStatus = (readStatus && GetMoleculeInfo(ppProduct2, pdt2Name, pdt2Type)) ;
          bPdt2 = true ;
        }
      }
      //
      // Create a new Reaction.  For association & exchange reactions, if rct1Type == reactant,
      // bool is true and rct1 is the pseudoisomer.  if not, bool is false, and rct1 is the excess
      //
      Reaction *preaction ;
      if     (!bRct2 && bPdt1 && pdt1Type == "modelled" && !bPdt2){
        preaction = new IsomerizationReaction(m_pMoleculeManager, mEnv, mFlags, id) ;
      }
      else if( bRct2 && (rct1Type == "modelled" || rct2Type == "modelled" ) && bPdt1 && (pdt1Type == "sink")){
        preaction = new BimolecularSinkReaction(m_pMoleculeManager, mEnv, mFlags, id, (rct1Type == "excessReactant")) ;
      }
      else if( bRct2 && (rct1Type == "deficientReactant" || rct2Type == "deficientReactant" ) && bPdt1 && !bPdt2){
        preaction = new AssociationReaction(m_pMoleculeManager, mEnv, mFlags, id, (rct1Type != "deficientReactant")) ;
      }
      else if( bRct2 && (rct1Type == "modelled" || rct2Type == "modelled" ) && bPdt1 && (pdt1Type == "modelled") && !bPdt2){
        preaction = new PseudoIsomerizationReaction(m_pMoleculeManager, mEnv, mFlags, id, (rct1Type != "modelled")) ;
      }
      else if(!bRct2 && bPdt1 && (pdt1Type == "sink" || pdt2Type == "sink")){
        preaction = new IrreversibleUnimolecularReaction(m_pMoleculeManager, mEnv, mFlags, id) ;
      }
      else if( bRct2 && bPdt1 && (pdt1Type == "sink" || pdt2Type == "sink")){
        preaction = new IrreversibleExchangeReaction(m_pMoleculeManager, mEnv, mFlags, id, (rct1Type == "deficientReactant")) ;
      }
      else {
        cinfo << "Unknown reaction type.\n";
        return false ;
      }

      // The information of the products of a dissociation reaction is necessary, as in
      // the xml output, Mesmer needs to know the products to draw the potential energy
      // surface. In addition, for dissociation reaction with QM tunneling, Mesmer also
      // needs to know the barrier height on the products side.

      readStatus = readStatus && preaction->InitializeReaction(ppReac);
      if(!readStatus)
        cerr << endl << "UNSATISFACTORY REACTION" << endl; //but keep parsing

      //
      // Add reaction to map.
      //

      // Need to check if there is duplicate reaction name/species: CHL

      m_reactions.push_back(preaction) ;
    }

    return allReactionsOK;
  }

  //
  // Extract molecule information from XML stream.
  //
  bool ReactionManager::GetMoleculeInfo(PersistPtr pp, string& MolName, string& MolType)
  {
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) {
      cerr << "Ill formed molecule tag." << endl;
      return false;
    }
    const char *pmoltype(NULL);
    const char *pmolname = ppmol->XmlReadValue("ref", false);
    if(pmolname) {
      ErrorContext c(pmolname);
      MolName = pmolname;
      pmoltype = ppmol->XmlReadValue("me:type", optional);
      if(!pmoltype)
        pmoltype = ppmol->XmlReadValue("role");
      if(pmoltype)
        MolType = pmoltype;
    }
    return pmolname && pmoltype;
  }

  // Set Initial population for individual species and/or grains
  void ReactionManager::setInitialPopulation(PersistPtr anchor)
  {
    PersistPtr pp=anchor;
    double populationSum = 0.0;
    PersistPtr ppInitMol = pp->XmlMoveTo("molecule");
    if(!ppInitMol)
      ppInitMol = pp->XmlMoveTo("me:molecule");
    while(ppInitMol){
      string sRef = ppInitMol->XmlReadValue("ref");
      if(sRef.size()){ // if got the name of the molecule
        Molecule* pMolecule = m_pMoleculeManager->find(sRef) ;
        double population = ppInitMol->XmlReadDouble("me:population", optional);
        if(IsNan(population))
          population= ppInitMol->XmlReadDouble("population", optional);
        int grainIdx = ppInitMol->XmlReadInteger("me:grain",optional);
        if(grainIdx==0)
          grainIdx = ppInitMol->XmlReadInteger("grain",optional);
        if (!IsNan(population) && grainIdx==0){  // what to do if population, but not grain has been specified
          populationSum += population;
          pMolecule->getPop().setInitPopulation(population);
          ctest << "Initial population of " << pMolecule->getName() << " = " << population << endl;
        }
        else if(!IsNan(population) && !IsNan(grainIdx)){  // what to do if population and grain have been specified
          populationSum += population;       
          pMolecule->getPop().setInitGrainPopulation(grainIdx,population);
          ctest << "Initial population of grain " << grainIdx << " in " << pMolecule->getName() << " = " << population << endl;
        }
        else
          population = 0.0;
      }
      PersistPtr pp = ppInitMol->XmlMoveTo("molecule");  // note: the same molecule should occur more than once ONLY to specify different grain populations
      if(!pp)
        pp = ppInitMol->XmlMoveTo("me:molecule");
      ppInitMol = pp;
    }
    if (populationSum != 1.0){
      if (populationSum > 0.0){
        // Populations need to be Normalized.
      } else if (populationSum == 0.0){
        // Issue warning that there are no populations set, and calculate only rate coefficients.
        populationSum += 1.0;
        //      pPseudoIsomer->getPop().setInitPopulation(populationSum);
      } else {
        // Issue error and stop.
      }
    }
  }

  // Find a reaction from its id
  Reaction* ReactionManager::find(const std::string& id) const 
  {

    bool found(false) ;
    size_t idx(0) ;
    for (size_t i(0) ; i < m_reactions.size() && !found; i++) {
      idx = i ;
      found = (m_reactions[i]->getName() == id) ; 
    }

    return found ?  m_reactions[idx] : NULL ;

  }

pair<string,Reaction*> ReactionManager::getCompleteReactantsAndProducts(const string& modelledMols) const
{
    //assume mol ids do not contain any of =><+ or spaces
    string rname, pname, result;
    string::size_type pos = modelledMols.find_first_of(" <=>");
    if(pos!=string::npos)
      rname = modelledMols.substr(0, pos);
    pos = modelledMols.find_first_not_of(" <=>", pos);
    if(pos!=string::npos)
    {
      pname = modelledMols.substr(pos);
      //remove second and third product names (added back below) 
      pos = pname.find('+');
      if(pos!=string::npos)
        pname.erase(pos);
    }

    Reaction* r(NULL); //of reactants
    //Find the reactant isomer in some reaction's reactants and use complete reactants
    for (unsigned ir=0; ir<m_reactions.size(); ir++)
    {
      if(m_reactions[ir]->get_reactant()->getName()==rname)
      {
        result = m_reactions[ir]->getReactionString(Reaction::reactantsOnly);
        r = m_reactions[ir];
        break;
      }
    }

    //Find the product isomer in some reaction's products
    bool found(false);
    for (unsigned i=0; i<m_reactions.size() && !found; i++)
    {
      vector<Molecule*> products;
      m_reactions[i]->get_products(products);
      for(unsigned j=0; j<products.size() && !found; ++j)
      {
        if(products[j]->getName()==pname)
        {
          result += " => ";
          result += m_reactions[i]->getReactionString(Reaction::productsOnly);
          found = true;
        }
      }
    }

    if(!found)
    {
      //Find the product isomer in some reaction's reactants
      for (unsigned ir=0; ir<m_reactions.size(); ir++)
      {
        if(m_reactions[ir]->get_reactant()->getName()==pname)
        {
          result += " => ";
          result += m_reactions[ir]->getReactionString(Reaction::reactantsOnly);
          break;
        }
      }
    }
    return make_pair(result, r);
}

vector<Reaction*> ReactionManager::getReactionsWithExcessReactant()
{
  vector<Reaction*> reactions;
  for (unsigned ir=0; ir<m_reactions.size(); ir++)
  {
    if(m_reactions[ir]->get_concExcessReactant()!=0.0)
      reactions.push_back(m_reactions[ir]);
  }
  return reactions;
}


}//namespace


