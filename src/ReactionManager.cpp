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
  bool ReactionManager::addreactions(PersistPtr ppReacList)
  {
    PersistPtr ppReac = ppReacList->XmlMoveTo("reaction");
    while(ppReac)
    {
      //
      // Create a new Reaction.
      //
      Reaction *preaction = new Reaction(m_pMoleculeManager) ;

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
      m_reactions.push_back(preaction) ;

      ppReac = ppReac->XmlMoveTo("reaction");
    }

    return true;
  }

  int ReactionManager::Connectivity(Molecule* pReactant, Molecule* pProduct)
  {
    return -1;
  }

  void ReactionManager::BuildSystemCollisionOperator(const double beta, const MesmerEnv &mEnv)
  {
    //
    // Find all the unique wells and lowest zero point energy.
    //
    Reaction::isomerMap isomers ; // Maps the location of reactant collision operator in the system matrix.
    double minEnergy(0) ; //this is the minimum of ZPE amongst all wells
    Molecule *pBathGasMolecule = m_pMoleculeManager->get_BathGasMolecule();
    for (size_t i(0) ; i < size() ; ++i) {

      vector<CollidingMolecule *> unimolecularspecies ;
      m_reactions[i]->get_unimolecularspecies(unimolecularspecies) ;

      for (size_t i(0) ; i < unimolecularspecies.size() ; ++i) {

        CollidingMolecule *pCollidingMolecule = unimolecularspecies[i] ;

        if(isomers.find(pCollidingMolecule) == isomers.end()){ // New isomer
          isomers[pCollidingMolecule] = 0 ; //initialize to a trivial location
          minEnergy = min(minEnergy,pCollidingMolecule->get_zpe()) ;
        }
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

      double zpe = isomer->get_zpe() - minEnergy ; // cell zpe related with the minimum of all wells
      int grnZpe = int(zpe * KcalPerMolToRC / mEnv.GrainSize) ; //convert to grain
      //cout << "_2007_11_09__14_10_21_ grnZpe = " << grnZpe << ", zpe = " << zpe << ", minEnergy = " << minEnergy << endl;
      isomer->set_grnZpe(grnZpe) ; //set grain ZPE (related with the minimum of all wells)
      int colloptrsize = mEnv.MaxGrn  - grnZpe ;
      isomer->set_colloptrsize(colloptrsize) ;
      msize += colloptrsize ;

      isomer->initCollisionOperator(beta, pBathGasMolecule, mEnv) ;
      meanomega += isomer->get_collisionFrequency() ;
    }
    meanomega /= isomers.size();

    //
    // Find all source terms.
    //
    Reaction::sourceMap sources ; // Maps the location of source in the system matrix.
    for (size_t i(0) ; i < size() ; ++i) {

      CollidingMolecule *pseudoIsomer = m_reactions[i]->get_pseudoIsomer() ;

      if(pseudoIsomer && sources.find(pseudoIsomer) == sources.end()){ // New source
        ++msize ;
        sources[pseudoIsomer] = msize ;
      } // what happens if one species is both source and isomer? Does it occupy two locations in the collision matrix? CHL
    }

    cout << endl << "Size of the collision matrix: " << msize << endl << endl ;

    // Allocate space for system collision operator.

    m_pSystemCollisionOperator = new dMatrix(msize) ;

    // Insert collision operators for individual wells.
    for (isomeritr = isomers.begin() ; isomeritr != isomers.end() ; ++isomeritr) {

      CollidingMolecule *isomer = isomeritr->first ;
      int colloptrsize = isomer->get_colloptrsize() ;
      double omega = isomer->get_collisionFrequency() ;
      int idx = isomeritr->second ;

      isomer->copyCollisionOperator(m_pSystemCollisionOperator, colloptrsize, idx, omega/meanomega, mEnv) ;

    }

    // Add connecting rate coefficients.
    for (size_t i(0) ; i < size() ; ++i) {

        m_reactions[i]->AddMicroRates(m_pSystemCollisionOperator,isomers,sources,1.0/meanomega, mEnv) ;

    }
  }

  void ReactionManager::diagCollisionOperator()
  {
    // Allocate space for eigenvalues.
    const int smsize = int(m_pSystemCollisionOperator->size()) ;
    vector<double> rr(smsize, 0.0);

    m_pSystemCollisionOperator->diagonalize(&rr[0]) ;

    // Print out the first ten eigen values.

    cout << endl ;

    for (int i = smsize-10 ; i < smsize; ++i) {
      formatFloat(cout, rr[i], 6, 15) ;
      cout << endl ;
    }

    cout << endl ;
  }

}//namespace
