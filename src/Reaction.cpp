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

#include "Persistence.h"
#include <sstream>
#include <math.h>
#include "Reaction.h"
#include "Constants.h"

using namespace Constants ;

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
   if (m_kfmc != NULL) m_alloc.deallocate(m_kfmc, MAXCELL) ;
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
bool Reaction::ReadFromXML(TiXmlElement* pnReac) 
{
  m_pXmlEl = pnReac;
  if(GetMolRef(pnReac, "reactant",0))
  {
    GetMolRef(pnReac, "reactant",1);
    GetMolRef(pnReac, "product",0);
    GetMolRef(pnReac, "product",1);
    GetMolRef(pnReac, "me:transitionState",0);

    TiXmlElement* pnThresh = pnReac->FirstChildElement("me:threshold");
    if(pnThresh)
    {
      istringstream idata(pnThresh->GetText());
      idata >> m_E0 ;
    }

  //Check that there is a transition state and either a modelled reactant or product
    if((m_Reactant || m_Product) && m_TransitionState)
        return true;
  }
  cerr << "Difficulty with missing, unknown or ill-formed reactant of a reaction" << endl;
  return false;

}

///Reads
bool Reaction::GetMolRef(TiXmlElement* pnReac,const char* molRole, int index)
{
  const char* pRef;
  Molecule* pMol;
  bool isReactant = !strcmp(molRole,"reactant");

  TiXmlHandle rHandle(pnReac);
  TiXmlElement* pnMol = rHandle.FirstChildElement(molRole).ChildElement("molecule",index).ToElement();
  if(pnMol 
    && (pRef=pnMol->Attribute("ref")))
  {
    pMol = m_pMoleculeManager->find(pRef);
    if(!pMol)
    {
      cerr << "Unknown molecule: " << pRef <<endl;
      return false;
    }

    if(!strcmp(molRole,"me:transitionState"))
    {
      m_TransitionState = dynamic_cast<TransitionState*>(pMol);
      return true;
    }

    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol);
    if(!strcmp(molRole,"reactant"))
    {
      if(pColMol && m_Reactant==NULL)
        m_Reactant = pColMol;
      else
        m_Reactant2 = pMol;
    }
    else
    {
      if(pColMol && m_Product==NULL)
        m_Product = pColMol;
      else
        m_Product2 = pMol;
    }
    return true;
  }
  return false;
}
/*
    int bracket_count = 1 ;                         // When zero end of reactant.

    while (bracket_count > 0) {                     // Main input loop begins.

       //
       // Find first delinater character.
       //
       char c ;
       while (in.get(c) && c != '@' && c != '}' ) ;

       if (c == '}') {                              // End of data.
          bracket_count-- ;
          continue ;
       }

       string keyword ;
       while (in.get(c) && c != '{' )
          keyword += c ;

       if        (keyword == "Reactant") {          // Reactant name.

          string data = collectData(in) ;

          m_Reactant = m_pMoleculeManager->find(data) ;

          continue ;

       } else if (keyword == "Transition State") {  // Transition State name.

          string data = collectData(in) ;

          m_TransitionState = m_pMoleculeManager->find(data) ;

          continue ;


       } else if (keyword == "Threshold") {         // Threshold Energy.

          istringstream idata(collectData(in)) ;

          idata >> m_E0 ;

          continue ;

       }

    }                                               // Main input loop ends.
*/

//
// Calculate the forward microcanonical rate coefficients, using RRKM theory.
//
void Reaction::CalcMicroRateCoeffs() {

    double plancksConst = 1.0/2.998e+10 ; // Planck's constant.

    // Allocate space to hold Micro-canonical rate coefficients.
    // (Freed in destructor.)
    m_kfmc  = m_alloc.allocate(MAXCELL) ; 

    // Initialize microcanoincal rate coefficients.

    int i, j ;
    for (i = 0 ; i < MAXCELL ; ++i ) {
        m_kfmc[i] = 0.0 ;
    }

    // Allocate some work space for density of states.

    double *ddosTS = m_alloc.allocate(MAXCELL) ; // Transistion state density of states.
    double *ddos   = m_alloc.allocate(MAXCELL) ; // Density of states of equilibrim molecule.

    // Extract densities of states from molecules.

    m_Reactant->densityOfStates(ddos) ;
    m_TransitionState->densityOfStates(ddosTS) ;

    double SumOfStates  = 0.0 ;
    int thresholdEnergy = int(m_E0 * KCMLTOPCM) ;
    for (i = thresholdEnergy, j = 0 ; i < MAXCELL ; ++i, ++j ) {
   
        // Integrate transition state density of states.

        SumOfStates += ddosTS[j] ;

        // Calculate microcanonical rate coefficients using RRKM expression.

        m_kfmc[i]    = SumOfStates / (plancksConst*ddos[i]) ;
    }

    // Test microcanonical rate coefficients.

//    if ( get_verbosity() ) 
        testMicroRateCoeffs() ;

    // Free work space.

    m_alloc.deallocate(ddosTS, MAXCELL) ;
    m_alloc.deallocate(ddos,   MAXCELL) ;
    
}

//
// Access microcanonical rate coefficients - cell values are averaged
// to give grain values. This code is similar to that in Molecule.cpp
// and this averaging should be done there. SHR 19/Sep/2004.
//
void Reaction::get_MicroRateCoeffs(vector<double> &kmc, int ngrn) {

    // Extract density of states of equilibrium molecule.

    double *ddos = m_alloc.allocate(MAXCELL) ; 
    m_Reactant->densityOfStates(ddos) ;

    int igsz = MAXCELL/ngrn ;

    // Check that there are enough cells.

    if (igsz < 1) {
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
        for (int j(0) ; j < igsz ; j++, idx1++ ) 
            smt += ddos[idx1] ;

        // Calculate average energy of the grain if it contains sum states.

        if ( smt > 0.0 ) {

            double smat(0.0) ;
            for (int j(0) ; j < igsz ; j++, idx3++ ) 
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

    m_alloc.deallocate(ddos, MAXCELL) ;
}

//
// Test the forward microcanonical rate coefficients.
//
void Reaction::testMicroRateCoeffs() {

    cout << endl << "Test of microcanonical rate coefficients" << endl << endl ;
  string comment("Microcanonical rate coefficients");
    TiXmlElement* list = WriteMainElement( m_pXmlEl, "me:microRateList", comment );

    // Allocate some work space for density of states.

    double *decll = m_alloc.allocate(MAXCELL) ;
    double *ddos  = m_alloc.allocate(MAXCELL) ;

    m_Reactant->cellEnergies(decll) ;
    m_Reactant->densityOfStates(ddos) ;

    for ( int n = 0 ; n < 29 ; ++n ) {

       double temp = 100.0*static_cast<double>(n + 2) ;
       double beta = 1.0/(0.695029*temp) ;

       double sm1 = 0.0 ;
       double sm2 = 0.0 ;
       double tmp = 0.0 ;
       for ( int i = 0 ; i < MAXCELL ; ++i ) {
           tmp  = ddos[i] * exp(-beta * decll[i]) ;
           sm1 += m_kfmc[i] * tmp ;
           sm2 +=             tmp ;
       }
       sm1 /= sm2 ; 
       formatFloat(cout, temp, 6, 7) ;
       formatFloat(cout, sm1,  6, 15) ;
       cout << endl ;
 
       //Add to XML document
       TiXmlElement* item = WriteElement(list, "me:microRate");
       WriteValueElement(item, "me:T",   temp, 6);
       WriteValueElement(item, "me:val", sm1,  6) ;
 
    }

    // Free work space.

    m_alloc.deallocate(decll, MAXCELL) ;
    m_alloc.deallocate(ddos,  MAXCELL) ;

}
}//namespace
