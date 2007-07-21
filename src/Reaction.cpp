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
#include "SimpleRRKM.h"
#include "SimpleILT.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
	Reaction::Reaction(MoleculeManager *pMoleculeManager): m_pMoleculeManager(pMoleculeManager),
		m_Reactant(NULL), 
		m_Reactant2(NULL), 
		m_Product(NULL), 
		m_Product2(NULL),
		m_kfwd(0.0),
		m_kfmc(NULL),
		m_kfgrn(),
        m_pMicroRateCalculator(){}

	Reaction::~Reaction() 
	{
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

		//Read products ... if any.
		pMol2=NULL;
		PersistPtr ppProduct1 = ppReac->MoveTo("product");
		if (ppProduct1) {   
			pMol1 = GetMolRef(ppProduct1);

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
		}

		// Classify the reaction.

		if (m_Reactant && m_Product && !m_Reactant2 && !m_Product2)
			m_reactiontype = ISOMERIZATION ;  
		else if (m_Reactant && m_Product && m_Reactant2 && !m_Product2)
			m_reactiontype = ASSOCIATION ;  
		else if (m_Reactant)
			m_reactiontype = DISSOCIATION ; 
		else {
			m_reactiontype = ERROR_REACTION ; 
			cerr << "Unknown combination of reactants and products" << endl ;
			return false;
		}

		// Determine the method of MC rate coefficient calculation.

		bool bCalculator(false) ;
		const char* pMCRCMethodtxt = ppReac->ReadValue("me:MCRCMethod") ;
		if(pMCRCMethodtxt)
		{
			string strMCRCMethodtxt(pMCRCMethodtxt) ;
			if        (strMCRCMethodtxt == "Simple RRKM") {

				SimpleRRKM *pSimpleRRKM = new SimpleRRKM(m_pMoleculeManager, m_Reactant) ;
				if(pSimpleRRKM->Initialize(ppReac)){ 
					m_pMicroRateCalculator = auto_ptr<MicroRateCalculator>(pSimpleRRKM) ;
					bCalculator = true ;
				}

			} else if (strMCRCMethodtxt == "Simple ILT") {

				SimpleILT *pSimpleILT = new SimpleILT(m_pMoleculeManager, m_Reactant) ;
				if(pSimpleILT->Initialize(ppReac)){ 
					m_pMicroRateCalculator = auto_ptr<MicroRateCalculator>(pSimpleILT) ;
					bCalculator = true ;
				}

			} else {
				cerr << "Unknown method for the determination of Microcanonical rate coefficients" << endl;
			}
		}

		return bCalculator ;
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
		if(m_Reactant2 == NULL){                // Possible dissociation or isomerization.
			unimolecularspecies.push_back(m_Reactant) ;
		}

		if(m_Product && m_Product2 == NULL){	// Possible association or isomerization.
			unimolecularspecies.push_back(m_Product) ;
		}
	}


	//
	// Access microcanoincal rate coeffcients. 
	//
	void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc) {
		calcGrnAvrgMicroRateCoeffs();

		kmc = m_kfgrn ;
	}

	//
	// Calculate grain averaged microcanoincal rate coeffcients. 
	//
	void Reaction::calcGrnAvrgMicroRateCoeffs() {

		// Calculate microcanonical rate coefficients.
		if (m_kfmc.size()==0) 
			m_pMicroRateCalculator->calculateMicroRateCoeffs(m_kfmc) ;

		// Calculate Grain averages of microcanonicla rate coefficients.
		if (m_kfgrn.size()==0)
			grnAvrgMicroRateCoeffs() ;
	}

	//
	// Access microcanonical rate coefficients - cell values are averaged
	// to give grain values. This code is similar to that in Molecule.cpp
	// and this averaging should be done there. SHR 19/Sep/2004.
	//
	void Reaction::grnAvrgMicroRateCoeffs() {

		int ngrn = pSys->MAXGrn();
		m_kfgrn.resize(ngrn);

		// Extract density of states of equilibrium molecule.

		vector<double> ddos(pSys->MAXCell(),0.0) ; 
		m_Reactant->cellDensityOfStates(&ddos[0]) ;

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

				m_kfgrn[idx2] = smat/smt ;
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
	}

	//
	// Add microcanonical terms to collision operator
	//
	void Reaction::AddMicroRates(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

                // Calculate Microcanonical rate coefficients.

		calcGrnAvrgMicroRateCoeffs() ;

		// Add microcanonical rates to the collision operator.

		switch(m_reactiontype) {
			case ISOMERIZATION :

				AddIsomerReactionTerms(CollOptr, isomermap, rMeanOmega) ;

				break;
			case ASSOCIATION :

				AddAssocReactionTerms(CollOptr, isomermap, rMeanOmega) ;

				break;
			case DISSOCIATION :

				AddDissocReactionTerms(CollOptr, isomermap, rMeanOmega) ;

				break;
			default :
				cerr << "Unknown reaction type" << endl;
		}
	} 

	//
	// Add isomer reaction terms to collision matrix.
	//
	void Reaction::AddIsomerReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

		// Locate isomers in system matrix.

		const int rctLocation = isomermap[m_Reactant] ;
		const int pdtLocation = isomermap[m_Product] ;

		// Get densities of states for detailed balance.

		const int ngrn = pSys->MAXGrn();
		vector<double> rctDos(ngrn, 0.0) ;  
		vector<double> pdtDos(ngrn, 0.0) ;  

		m_Reactant->grnDensityOfStates(rctDos) ;
		m_Product->grnDensityOfStates(pdtDos) ;

		const int idx = m_Product->get_grnZpe() - m_Reactant->get_grnZpe() ;
		for ( int i = max(0,-idx) ; i < min(ngrn,(ngrn-idx)) ; ++i ) {
			int ll = i + idx ;
			int ii(rctLocation + ll) ;
			int jj(pdtLocation + i) ;
			(*CollOptr)[ii][ii] -= rMeanOmega * m_kfgrn[ll] ;                            // Forward loss reaction.
			(*CollOptr)[jj][jj] -= rMeanOmega * m_kfgrn[ll]*rctDos[ll]/pdtDos[i] ;       // Backward loss reaction from detailed balance.
			(*CollOptr)[ii][jj]  = rMeanOmega * m_kfgrn[ll]*sqrt(rctDos[ll]/pdtDos[i]) ; // Reactive gain.
			(*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                 // Reactive gain.
		}

	}

	//
	// Add (reversible) association reaction terms to collision matrix.
	//
	void Reaction::AddAssocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {
	}

	//
	// Add dissociation reaction terms to collision matrix.
	//
	void Reaction::AddDissocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

		// Locate reactant in system matrix.

		const int rctLocation = isomermap[m_Reactant] ;

		for ( int i = 0 ; i < m_Reactant->get_colloptrsize() ; ++i ) {
			int ii(rctLocation + i) ;
			(*CollOptr)[ii][ii] -= rMeanOmega * m_kfgrn[i] ;                            // Forward loss reaction.
		}
	}

}//namespace
