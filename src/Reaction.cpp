//-------------------------------------------------------------------------------------------
//
// Reaction.cpp 
//
// Author: Struan Robertson 
// Date:   23/Feb/2003
//
// This file contains the implementation of the Reaction class.
//
//-------------------------------------------------------------------------------------------

#include <math.h>
#include "system.h"
#include "Constants.h"
#include "Persistence.h"
#include "Reaction.h"

using namespace Constants ;
using namespace std;

namespace mesmer
{
Reaction::Reaction(MoleculeManager *pMoleculeManager): m_pMoleculeManager(pMoleculeManager),
                                                       m_Reactant(NULL), 
                                                       m_Reactant2(NULL), 
                                                       m_Product(NULL), 
                                                       m_Product2(NULL) {}

Reaction::~Reaction() 
{
    // Free any memory assigned for calculating densities of states. 
   if (m_kfmc != NULL) m_alloc.deallocate(m_kfmc, pSys->MAXCell()) ;
}
/*
Reaction::Reaction(const Reaction& reaction) {
    // Copy constructor - define later SHR 23/Feb/2003
}

Reaction& Reaction::operator=(const Reaction& reaction) {
    // Assignment operator - define later SHR 23/Feb/2003
    
    return *this ;
}
*/ 
//
// Read the Molecular data from inout stream.
//
bool Reaction::Initialize(PersistPtr ppReac) 
{
  m_ppPersist = ppReac;

  //Read reaction ID
  const char* id = ppReac->ReadValue("id");
  if(id)
    m_Name = id; //Continues if reaction id not found

  Molecule* pMol1,*pMol2=NULL;
  //Read reactants
  PersistPtr ppReactant1  = ppReac->MoveTo("reactant");
  pMol1 = GetMolRef(ppReactant1);
  if(!pMol1) return false;

  PersistPtr ppReactant2  = ppReactant1->MoveTo("reactant");
  if(ppReactant2)
  {
    pMol2 = GetMolRef(ppReactant2);
    if(!pMol2) return false;
    m_Reactant2 = pMol2;
  }
  //Put the Colliding Molecule into m_Reactant, even if it is second in datafile
  CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
  if(pColMol)
  {
    m_Reactant = pColMol;
    m_Reactant2 = pMol2;
  }
  else
  {  
    pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
    if(!pColMol)
    {
      cerr << "Either " << pMol1->getName() << " or " 
        << pMol2->getName() <<" has to be a modelled molecule" <<endl;
      return false;
    }
    m_Reactant = pColMol;
    m_Reactant2 = pMol1;
  }

  //Read products
  pMol2=NULL;
  PersistPtr ppProduct1  = ppReac->MoveTo("product");
  pMol1 = GetMolRef(ppProduct1);
  if(!pMol1) return false; //there must be at least one product

  PersistPtr ppProduct2  = ppProduct1->MoveTo("product");
  pMol2 = GetMolRef(ppProduct2);

  //Put the Colliding Molecule into m_Product, even if it is second in datafile
  pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
  if(pColMol)
  {
    m_Product = pColMol;
    m_Product2 = pMol2;
  }
  else
  {  
    pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
    if(!pColMol)
    {
      cerr << "Either " << pMol1->getName() << " or " 
        << pMol2->getName() <<" has to be a modelled molecule" <<endl;
      return false;
    }
    m_Product = pColMol;
    m_Product2 = pMol1;
  }

  //Read TransitionState
  pMol1 = GetMolRef(ppReac->MoveTo("me:transitionState"));
  TransitionState* pTSMol = dynamic_cast<TransitionState*>(pMol1);
  if(!pTSMol)
  {
    cerr << "The molecule " << pMol1->getName() << " needs to be of transition state type" << endl;
    return false;
  }
  m_TransitionState = pTSMol;

  const char* pthreshtxt = ppReac->ReadValue("me:threshold",false);
  if(pthreshtxt)
  {
    stringstream ss(pthreshtxt);
    ss >> m_E0;
  }

//Check that there is a transition state and either a modelled reactant or product
  if((m_Reactant || m_Product) && m_TransitionState)
   return true;

  cerr << "Difficulty with missing, unknown or ill-formed reactant, product or TS of a reaction" << endl;
  return false;
}


Molecule* Reaction::GetMolRef(PersistPtr pp)
{
  Molecule* pMol;
  if(!pp)
    return NULL;
  PersistPtr ppmol = pp->MoveTo("molecule");
  if(!ppmol) return false;
  const char* pRef = ppmol->ReadValue("ref");
  if(pRef)
    pMol = m_pMoleculeManager->find(pRef);
  if(!pMol)
  {
    cerr << "Unknown molecule: " << pRef <<endl;
    return NULL;
  }
  return pMol;
}

//
// Returns the unimolecular species in each reaction, i.e. for association
// (source term) or dissociation (sink term) reaction one species is returned,
// for an isomerization reaction two species are returned.
//
void Reaction::get_unimolecularspecies(vector<CollidingMolecule *> &unimolecularspecies) const
{
	if(m_Reactant2 == NULL){ // Possible dissociation or isomerization.
		unimolecularspecies.push_back(m_Reactant) ;
	}

	if(m_Product2 == NULL){	// Possible association or isomerization.
		unimolecularspecies.push_back(m_Product) ;
	}
}

//
// Calculate the forward microcanonical rate coefficients, using RRKM theory.
//
void Reaction::CalcMicroRateCoeffs() {

    double plancksConst = 1.0/2.998e+10 ; // Planck's constant.

    // Allocate space to hold Micro-canonical rate coefficients.
    // (Freed in destructor.)
    m_kfmc  = m_alloc.allocate(pSys->MAXCell()) ; 

    // Initialize microcanoincal rate coefficients.

    int i, j ;
    for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
        m_kfmc[i] = 0.0 ;
    }

    // Allocate some work space for density of states.

    double *ddosTS = m_alloc.allocate(pSys->MAXCell()) ; // Transistion state density of states.
    double *ddos   = m_alloc.allocate(pSys->MAXCell()) ; // Density of states of equilibrim molecule.

    // Extract densities of states from molecules.

    m_Reactant->densityOfStates(ddos) ;
    m_TransitionState->densityOfStates(ddosTS) ;

    double SumOfStates  = 0.0 ;
    int thresholdEnergy = int(m_E0 * KCMLTOPCM) ;
    for (i = thresholdEnergy, j = 0 ; i < pSys->MAXCell() ; ++i, ++j ) {
   
        // Integrate transition state density of states.

        SumOfStates += ddosTS[j] ;

        // Calculate microcanonical rate coefficients using RRKM expression.

        m_kfmc[i]    = SumOfStates / (plancksConst*ddos[i]) ;
    }

    // Test microcanonical rate coefficients.

//    if ( get_verbosity() ) 
        testMicroRateCoeffs() ;

    // Free work space.

    m_alloc.deallocate(ddosTS, pSys->MAXCell()) ;
    m_alloc.deallocate(ddos,   pSys->MAXCell()) ;
    
}

//
// Access microcanonical rate coefficients - cell values are averaged
// to give grain values. This code is similar to that in Molecule.cpp
// and this averaging should be done there. SHR 19/Sep/2004.
//
void Reaction::get_MicroRateCoeffs(vector<double> &kmc, int ngrn) {

    // Extract density of states of equilibrium molecule.

    double *ddos = m_alloc.allocate(pSys->MAXCell()) ; 
    m_Reactant->densityOfStates(ddos) ;

//    int igsz = pSys->MAXCell()/ngrn ;

    // Check that there are enough cells.

    if (pSys->igsz() < 1) {
       cout << "     ********* Not enought Cells to produce ************" << endl
            << "     ********* requested number of Grains.  ************" << endl ;
       exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i(0) ; i < ngrn ; i++ ) {

        int idx3 = idx1 ;

        // Calculate the number of states in a grain.

        double smt(0.0) ;
        for (int j(0) ; j < pSys->igsz() ; j++, idx1++ ) 
            smt += ddos[idx1] ;

        // Calculate average energy of the grain if it contains sum states.

        if ( smt > 0.0 ) {

            double smat(0.0) ;
            for (int j(0) ; j < pSys->igsz() ; j++, idx3++ ) 
                smat += m_kfmc[idx3] * ddos[idx3] ;
             
            kmc[idx2] = smat/smt ;
            idx2++ ;
        }
    }

    // Issue warning if number of grains produced is less that requested.

    if ( idx2 < ngrn ) {
       cout <<  endl
            <<  "     WARNING: Number of grains produced is less than requested" << endl
            <<  "     Number of grains requested: " << ngrn << endl
            <<  "     Number of grains produced : " << idx2 << endl ;
    }

    // Free work space.

    m_alloc.deallocate(ddos, pSys->MAXCell()) ;
}

//
// Test the forward microcanonical rate coefficients.
//
void Reaction::testMicroRateCoeffs() {

    cout << endl << "Test of microcanonical rate coefficients" << endl << endl ;
  string comment("Microcanonical rate coefficients");
    PersistPtr ppList = m_ppPersist->WriteMainElement("me:microRateList", comment );

    // Allocate some work space for density of states.

    double *decll = m_alloc.allocate(pSys->MAXCell()) ;
    double *ddos  = m_alloc.allocate(pSys->MAXCell()) ;

    m_Reactant->cellEnergies(decll) ;
    m_Reactant->densityOfStates(ddos) ;

    for ( int n = 0 ; n < 29 ; ++n ) {

       double temp = 100.0*static_cast<double>(n + 2) ;
       double beta = 1.0/(0.695029*temp) ;

       double sm1 = 0.0 ;
       double sm2 = 0.0 ;
       double tmp = 0.0 ;
       for ( int i = 0 ; i < pSys->MAXCell() ; ++i ) {
           tmp  = ddos[i] * exp(-beta * decll[i]) ;
           sm1 += m_kfmc[i] * tmp ;
           sm2 +=             tmp ;
       }
       sm1 /= sm2 ; 
       formatFloat(cout, temp, 6, 7) ;
       formatFloat(cout, sm1,  6, 15) ;
       cout << endl ;
 
       //Add to XML document
       PersistPtr ppItem = ppList->WriteElement("me:microRate");
       ppItem->WriteValueElement("me:T",   temp, 6);
       ppItem->WriteValueElement("me:val", sm1,  6) ;
 
    }

    // Free work space.

    m_alloc.deallocate(decll, pSys->MAXCell()) ;
    m_alloc.deallocate(ddos,  pSys->MAXCell()) ;

}
}//namespace
