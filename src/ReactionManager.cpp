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

     //
	 // Update connectivity table.
	 //
	 m_ReactionConnectivity.addreaction(preaction) ;

     ppReac = ppReac->MoveTo("reaction");
  }

  m_ReactionConnectivity.printconnectivitytable() ;

  return true;
}

int ReactionManager::Connectivity(Molecule* pReactant, Molecule* pProduct)
{
  return -1;
}

void ReactionManager::BuildSystemCollisionOperator(const double beta)
{
    // Find all the individual wells.
    vector<string> species ;
    m_ReactionConnectivity.get_unimolecularspecies(species)  ;

    // Initialize all collision operators and find well with the lowest energy.

    vector<double> zpes ;
    vector<CollidingMolecule *> isomers ;
    double minEnergy(0) ;
    for (size_t i(0) ; i < species.size() ; i++) {

        Molecule* pmolecule ;
        try {

            pmolecule = m_pMoleculeManager->find(species[i]) ;

        } catch (...) {
            cout << "Molecule does not exist in Molecular database." << endl ;
            exit (1) ;
        }
        //
        // SHR 2/Jan/2007: The following line has a down cast which is bad news,
        // we should review this as soon as possible.
        //
        CollidingMolecule* pcollidingmolecule = dynamic_cast<CollidingMolecule *>(pmolecule) ;
        isomers.push_back(pcollidingmolecule) ;
        pcollidingmolecule->initCollisionOperator(beta) ;
        
        double zpe = pcollidingmolecule->get_zpe() ;
        zpes.push_back(zpe) ;
        minEnergy = min(minEnergy,zpe) ;
    }

    // Shift all wells to the same origin and calculate the size of 
    // the system collision operator.

    vector<int> colloptrsizes ;
    int msize(0) ;
    for (size_t i(0) ; i < zpes.size() ; i++) {
        zpes[i] -= minEnergy ;
        int colloptrsize = MAXGRN  - int(zpes[i] * KCMLTOPCM / double(MAXCELL/MAXGRN)) ;
        msize += colloptrsize ;
        colloptrsizes.push_back(colloptrsize) ;
    }

    cout << endl << "The size of the collision matrix is: " << msize << endl << endl ;

    // Allocate space for system collision operator.
    
    m_pSystemCollisionOperator = new dMatrix(msize) ;

    // Insert collision operators for individual wells.

    int idx(0) ;
    for (size_t m(0) ; m < colloptrsizes.size() ; m++) {
        int colloptrsize = colloptrsizes[m] ;
//        const dMatrix *colloptr = isomers[m]->collisionOperator() ;

        isomers[m]->copyCollisionOperator(m_pSystemCollisionOperator, colloptrsize, idx) ;

/*      for (int i(0) ; i < colloptrsize ; i++) {
            int ii(idx + i) ;
            for (int j(0) ; j < colloptrsize ; j++) {
                int jj(idx + j) ;
                (*m_pSystemCollisionOperator)[ii][jj] = (*colloptr)[i][j] ;
            }   
        }    */

        idx += colloptrsize ;
    }

    // Add connecting rate coefficients.
}

void ReactionManager::diagCollisionOperator()
{
    // Allocate space for eigenvalues.
    vector<double> rr(m_pSystemCollisionOperator->size(), 0.0) ;

    m_pSystemCollisionOperator->diagonalize(&rr[0]) ;
}

}//namespace
