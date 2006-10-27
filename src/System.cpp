//-------------------------------------------------------------------------------------------
//
// System.cpp 
//
// Author: Struan Robertson 
// Date:   11/Feb/2003
//
// This file contains the implementation of the System class.
//
//-------------------------------------------------------------------------------------------
#include <iostream>
#include "System.h"
#include "Constants.h"
#include "formatfloat.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
System::System(): m_pMoleculeManager(0), m_pReactionManager(0) {

     m_pMoleculeManager = new MoleculeManager() ;

     m_pReactionManager = new ReactionManager(m_pMoleculeManager) ;

 }

System::~System() { }

//
// Parse an input data file.
//
bool System::parse(TiXmlElement* root)
{
  TiXmlElement* mollist = root->FirstChildElement("moleculeList");
  if(!mollist)
  {
    cerr << "No molecules have been specified" << endl;
    return false;
  }if(!m_pMoleculeManager->addmols(mollist))
      return false;;

  TiXmlElement* reaclist = root->FirstChildElement("reactionList");
  if(!reaclist)
  {
    cerr << "No reactions have been specified" << endl;
    return false;
  }if(!m_pReactionManager->addreactions(reaclist))
    return false;

  TiXmlElement* pnConditions = root->FirstChildElement("me:conditions");
  if(!pnConditions)
  {
    cerr << "No conditions have been specified" << endl;
    return false;
  }
  TiXmlElement* pnbathGas = pnConditions->FirstChildElement("me:bathGas");
  if(!pnbathGas || 
    ( !(m_pBathGasMolecule = dynamic_cast<BathGasMolecule*>(m_pMoleculeManager->find(pnbathGas->GetText())))))
  {
    cerr << "No bath gas specified" << endl;
    return false;
  }

  TiXmlElement* pnTemp = pnConditions->FirstChildElement("me:temperature");
  TiXmlElement* pnConc = pnConditions->FirstChildElement("me:conc");
  const char* Ttxt, *Ctxt;
  if(!pnTemp || !(Ttxt = pnTemp->GetText()) || !pnConc || !(Ctxt = pnConc->GetText()))
  {
    cerr << "The bath gas temperature or concentration has not been specified" << endl;
    return false;
  }
  istringstream Tss(Ttxt);
  Tss >> temp;
  istringstream Css(Ctxt);
  Css >> conc;
  return true;
}

//
// Begin calculation.
//
void System::calculate() 
{ 
	//
	// Collision parameters for Argon. (Temporary hardwiring SHR 18/May/2003).
	//
//double bthMass    = 39.95 ;
//double bthSigma   = 3.47 ;
//double bthEpsilon = 114 ;
  double bthMass    = m_pBathGasMolecule->getMass();
  double bthSigma   = m_pBathGasMolecule->getSigma();
  double bthEpsilon = m_pBathGasMolecule->getEpsilon();

//double temp = 300.0 ;
	double beta = 1.0/(boltzmann*temp) ;

//double conc = 1.0e8 ;

	// Build up collison matrix by looping over reactions.

	for (size_t i=0; i < m_pReactionManager->size() ; i++) {

		Reaction *reaction = (*m_pReactionManager)[i] ;

		reaction->CalcMicroRateCoeffs() ;

		// Work space to hold mircocanonical rates.

		vector<double> kmc(MAXGRN,0.0) ;

		reaction->get_MicroRateCoeffs(kmc, MAXGRN) ;

		CollidingMolecule *pmolecule = reaction->m_Reactant ;

		pmolecule->initCollisionOperator(beta) ;

		double omega = pmolecule->collisionFrequency(temp, conc, bthMass, bthSigma, bthEpsilon) ;

		// for (int i(0) ; i < MAXCELL ; i++) 
		//     kmc[i] /= omega ;

		pmolecule->diagCollisionOperator() ;

		// Calculate matrix elements

		double kinf = pmolecule->matrixElement(MAXGRN-1,MAXGRN-1,kmc,MAXGRN) ;

		cout << endl ;
		formatFloat(cout, kinf, 6, 15) ;

	}    

}
}//namespace
