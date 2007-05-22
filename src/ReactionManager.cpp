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
#include <algorithm>

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
    // To begin with we need to: 
	// 1) Find all the individual wells.
    // 2) Initialize all collision operators. 
	// 3) Find well with the lowest energy.
	// 4) Calculate the mean collision frequency.

    vector<double> zpes ;
    vector<double> omegas ;
    vector<CollidingMolecule *> isomers ;
    double minEnergy(0) ;
	Molecule *pBathGasMolecule = m_pMoleculeManager->get_BathGasMolecule();
    for (size_t i(0) ; i < size() ; i++) {

		vector<CollidingMolecule *> unimolecularspecies ;
		m_reactions[i]->get_unimolecularspecies(unimolecularspecies) ;

		for (size_t i(0) ; i < unimolecularspecies.size() ; i++) {

			CollidingMolecule *pcollidingmolecule = unimolecularspecies[i] ;

			if(find(isomers.begin(),isomers.end(),pcollidingmolecule) == isomers.end()){
				isomers.push_back(pcollidingmolecule) ;

				pcollidingmolecule->initCollisionOperator(beta) ;
	        
				double zpe = pcollidingmolecule->get_zpe() ;
				zpes.push_back(zpe) ;

				minEnergy = min(minEnergy,zpe) ;

				omegas.push_back(pcollidingmolecule->collisionFrequency(beta, conc, pBathGasMolecule)) ;
			}
		}
    }

    // Shift all wells to the same origin and calculate the size of 
    // the system collision operator.

    vector<int> colloptrsizes ;
    int msize(0) ;
	double meanomega(0.0) ;
    for (size_t i(0) ; i < isomers.size() ; i++) {
        zpes[i] -= minEnergy ;
        int colloptrsize = MAXGRN  - int(zpes[i] * KCMLTOPCM / double(MAXCELL/MAXGRN)) ;
        msize += colloptrsize ;
        colloptrsizes.push_back(colloptrsize) ;
		meanomega += omegas[i] ;
    }
	meanomega /= isomers.size();

    cout << endl << "The size of the collision matrix is: " << msize << endl << endl ;

    // Allocate space for system collision operator.
    
    m_pSystemCollisionOperator = new dMatrix(msize) ;

    // Insert collision operators for individual wells.

    int idx(0) ;
    for (size_t m(0) ; m < isomers.size() ; m++) {
        int colloptrsize = colloptrsizes[m] ;

        isomers[m]->copyCollisionOperator(m_pSystemCollisionOperator, colloptrsize, idx, omegas[m]/meanomega) ;

        idx += colloptrsize ;
    }

    // Add connecting rate coefficients.
    for (size_t i(0) ; i < size() ; i++) {

		vector<CollidingMolecule *> unimolecularspecies ;
		m_reactions[i]->get_unimolecularspecies(unimolecularspecies) ;

		if (unimolecularspecies.size() == 2){ 
			// Isomerization
			CollidingMolecule *reactant = unimolecularspecies[0] ; 

		} else if (unimolecularspecies.size() == 1) { 
			//Association/Dissociation
		} else {
			// Bimolecular reaction
		}

    }
}

void ReactionManager::diagCollisionOperator()
{
    // Allocate space for eigenvalues.
    vector<double> rr(m_pSystemCollisionOperator->size(), 0.0) ;

    m_pSystemCollisionOperator->diagonalize(&rr[0]) ;
}

}//namespace
