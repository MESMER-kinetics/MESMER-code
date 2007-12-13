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

using namespace Constants ;
using namespace std ;

namespace mesmer
{
  //
  // Add a new reaction to the map.
  //
  bool ReactionManager::addreactions(PersistPtr ppReacList, const MesmerEnv& Env)
  {
    PersistPtr ppReac = ppReacList->XmlMoveTo("reaction");
    while(ppReac)
    {
      //
      // Create a new Reaction.
      //
      Reaction *preaction = new Reaction(m_pMoleculeManager, Env) ;

      //
      // Initialize Reaction from input stream.
      //
      if(!preaction->InitializeReaction(ppReac)){
        delete preaction;
        return false;
      }

      //preaction->put_verbosity(true) ;

      //
      // Add reaction to map.
      //

      //need to check if there is duplicate reaction name/species: CHL


      m_reactions.push_back(preaction) ;

      ppReac = ppReac->XmlMoveTo("reaction");
    }

    return true;
  }

  int ReactionManager::Connectivity(Molecule* pReactant, Molecule* pProduct)
  {
    return -1;
  }

  bool ReactionManager::BuildSystemCollisionOperator(const double beta, const MesmerEnv &m_Env)
  {
    //
    // Find all the unique wells and lowest zero point energy.
    //
    Reaction::isomerMap isomers ; // Maps the location of reactant collision operator in the system matrix.
    double minEnergy(0) ; //this is the minimum of ZPE amongst all wells
    Molecule *pBathGasMolecule = m_pMoleculeManager->get_BathGasMolecule();

    // populate isomerMap with unimolecular species
    for (size_t i(0) ; i < size() ; ++i) {

      vector<ModelledMolecule *> unimolecules ;

      int flag1 = m_reactions[i]->get_unimolecularspecies(unimolecules) ;

      if (flag1 > 0){ // association, dissociation or isomerization
        for (size_t i(0) ; i < unimolecules.size() ; ++i) {
          CollidingMolecule *pCollidingMolecule = dynamic_cast<CollidingMolecule*>(unimolecules[i]) ;
          if(isomers.find(pCollidingMolecule) == isomers.end()){ // New isomer
            isomers[pCollidingMolecule] = 0 ; //initialize to a trivial location
            minEnergy = min(minEnergy, pCollidingMolecule->get_zpe()) ;
          }
        }
      }
      else{
        stringstream errorMsg;
        // Does this kind of reaction ever exist? CHL
        errorMsg << "Reaction " << m_reactions[i]->getName() << " is invalid without any product and reactant." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      }
    }

    //
    // Shift all wells to the same origin, calculate the size of the system collision operator,
    // calculate the mean collision frequency and initialize all collision operators.
    //
    int msize(0) ; // size of the collision matrix
    double meanomega(0.0) ;
    Reaction::isomerMap::iterator isomeritr = isomers.begin() ;
    for (; isomeritr != isomers.end() ; ++isomeritr) {

      CollidingMolecule *isomer = isomeritr->first ;
      isomeritr->second = msize ; //set location

      double zpe = (isomer->get_zpe()) - minEnergy ; // cell zpe related with the minimum of all wells
      int grnZpe = int(zpe * KcalPerMolToRC / m_Env.GrainSize) ; //convert to grain

      isomer->set_grnZpe(grnZpe) ; //set grain ZPE (related with the minimum of all wells)
      int colloptrsize = m_Env.MaxGrn  - grnZpe ;
      isomer->set_colloptrsize(colloptrsize) ;
      msize += colloptrsize ;

      isomer->initCollisionOperator(beta, pBathGasMolecule) ;
      meanomega += isomer->get_collisionFrequency() ;
    }
    meanomega /= isomers.size();

    //
    // Find all source terms. 
    // Note: 1. A source term is probably the only deficient reactant that initiates the whole process of reactions 
    //          in the master equation. In this case we think there may be more than one source terms.
    //       2. In the current construction of Mesmer, the source is a SuperMolecule representing both reactants.
    Reaction::sourceMap sources ; // Maps the location of source in the system matrix.
    for (size_t i(0) ; i < size() ; ++i) {

      SuperMolecule *pSuperMolecule = NULL; // get the sourceMap returnHere
      pSuperMolecule = m_reactions[i]->get_bi_molecularspecies();
      if (pSuperMolecule && sources.find(pSuperMolecule) == sources.end()){ // New source
        sources[pSuperMolecule] = msize ;
        ++msize ;
      }
    }

    {
      stringstream errorMsg;
      errorMsg << "Size of the collision matrix: " << msize;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    // Allocate space for system collision operator.

    m_pSystemCollisionOperator = new dMatrix(msize) ;

    // Insert collision operators for individual wells.
    for (isomeritr = isomers.begin() ; isomeritr != isomers.end() ; ++isomeritr) {

      CollidingMolecule *isomer = isomeritr->first ;
      int colloptrsize = isomer->get_colloptrsize() ;
      double omega = isomer->get_collisionFrequency() ;
      int idx = isomeritr->second ;

      isomer->copyCollisionOperator(m_pSystemCollisionOperator, colloptrsize, idx, omega/meanomega) ;

    }

    // Add connecting rate coefficients.
    for (size_t i(0) ; i < size() ; ++i) {
        m_reactions[i]->AddMicroRates(m_pSystemCollisionOperator,isomers,sources,1.0/meanomega) ;
    }

    return true;
  }

  void ReactionManager::diagCollisionOperator()
  {
    // Allocate space for eigenvalues.
    const int smsize = int(m_pSystemCollisionOperator->size()) ;
    vector<double> rr(smsize, 0.0);
    {
      stringstream errorMsg;
      errorMsg << "Size of the collision operator = " << smsize;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    m_pSystemCollisionOperator->diagonalize(&rr[0]) ;

    // Print out the first ten eigen values.
    {
      stringstream errorMsg;
      errorMsg << "Printing the first ten eigen values";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    cout << endl ;

    for (int i = 0 ; i < smsize; ++i) {
      formatFloat(cout, rr[i], 6, 15) ;
      cout << endl ;
    }

    cout << endl ;
  }

}//namespace
