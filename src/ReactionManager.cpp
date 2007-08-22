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

#include <algorithm>
#include "System.h"


using namespace Constants ;
using namespace std ;

namespace mesmer
{
  // 
  // Add a new reaction to the map.
  //
  bool ReactionManager::addreactions(PersistPtr ppReacList)
  {
    PersistPtr ppReac = ppReacList->MoveTo("reaction");
    while(ppReac)
    {
      //
      // Create a new Reaction.
      //
      Reaction *preaction = new Reaction(m_pMoleculeManager) ;

      //
      // Initialize Reaction from input stream.
      //
      if(!preaction->Initialize(ppReac))
          return false;;

      //     preaction->put_verbosity(true) ;

      //
      // Add reaction to map.
      //
      m_reactions.push_back(preaction) ;

      ppReac = ppReac->MoveTo("reaction");
    }

    return true;
  }

  int ReactionManager::Connectivity(Molecule* pReactant, Molecule* pProduct)
  {
    return -1;
  }

  void ReactionManager::BuildSystemCollisionOperator(const double beta, const double conc)
  {
    // 
    // Find all the unique wells and lowest zero point energy.
    //
    Reaction::isomerMap isomers ; // Maps the location of reactant collision operator in the system matrix.
    double minEnergy(0) ;
    Molecule *pBathGasMolecule = m_pMoleculeManager->get_BathGasMolecule();
    for (size_t i(0) ; i < size() ; i++) {

      vector<CollidingMolecule *> unimolecularspecies ;
      m_reactions[i]->get_unimolecularspecies(unimolecularspecies) ;

      for (size_t i(0) ; i < unimolecularspecies.size() ; i++) {

        CollidingMolecule *pcollidingmolecule = unimolecularspecies[i] ;

        if(isomers.find(pcollidingmolecule) == isomers.end()){ // New isomer
          isomers[pcollidingmolecule] = 0 ;
          minEnergy = min(minEnergy,pcollidingmolecule->get_zpe()) ;
        }
      }
    }

		//
    // Shift all wells to the same origin, calculate the size of the system collision operator,
    // calculate the mean collision frequency and initialize all collision operators. 
    //
    int msize(0) ;
    double meanomega(0.0) ;
    Reaction::isomerMap::iterator isomeritr = isomers.begin() ;
    for (; isomeritr != isomers.end() ; ++isomeritr) {

      CollidingMolecule *isomer = isomeritr->first ;
      isomeritr->second = msize ;

      double zpe = isomer->get_zpe() - minEnergy ;
      int grnZpe = int(zpe * KCMLTOPCM / pSys->igsz()) ; 
      int colloptrsize = pSys->MAXGrn()  - grnZpe ;
      isomer->set_grnZpe(grnZpe) ;
      isomer->set_colloptrsize(colloptrsize) ;
      msize += colloptrsize ;

      isomer->initCollisionOperator(beta, conc, pBathGasMolecule) ;
      meanomega += isomer->get_collisionFrequency() ;
    }
    meanomega /= isomers.size();

		//
		// Find all source terms.
		//
    Reaction::sourceMap sources ; // Maps the location of source in the system matrix.
    for (size_t i(0) ; i < size() ; i++) {

      CollidingMolecule *pseudoIsomer = m_reactions[i]->get_pseudoIsomer() ;

      if(pseudoIsomer && sources.find(pseudoIsomer) == sources.end()){ // New source
        msize++ ;
        sources[pseudoIsomer] = msize ;
      }
    }

    cout << endl << "The size of the collision matrix is: " << msize << endl << endl ;

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
    for (size_t i(0) ; i < size() ; i++) {

        m_reactions[i]->AddMicroRates(m_pSystemCollisionOperator,isomers,sources,1.0/meanomega) ;

    }
  }

  void ReactionManager::diagCollisionOperator()
  {
    // Allocate space for eigenvalues.
    const int smsize = int(m_pSystemCollisionOperator->size()) ;
    vector<double> rr(smsize, 0.0) ;

    m_pSystemCollisionOperator->diagonalize(&rr[0]) ;

    // Print out the first ten eigen values.

    cout << endl ;

    for (int i = smsize-10 ; i < smsize; i++) {
      formatFloat(cout, rr[i], 6, 15) ; 
      cout << endl ;
    }

    cout << endl ;
  }

}//namespace
