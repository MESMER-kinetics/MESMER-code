//-------------------------------------------------------------------------------------------
//
// SimpleRRKM.cpp 
//
// Author: Struan Robertson 
// Date:   1/Jul/2007
//
// This file contains the implementation of the Simple RRKM MC rate coefficient method.
//
//-------------------------------------------------------------------------------------------

#include <math.h>
#include "system.h"
#include "Constants.h"
#include "Persistence.h"
#include "SimpleRRKM.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{

	bool SimpleRRKM::Initialize(PersistPtr &ppReac) {

		// Is there a simple transition state?

		PersistPtr ppTransitionState = ppReac->MoveTo("me:transitionState") ;
		if (ppTransitionState) {
			PersistPtr ppmol = ppTransitionState->MoveTo("molecule");
			if(!ppmol){ 
				cerr << "Cannot locate Transition state: " << endl;
				return false;
			}
			const char* pRef = ppmol->ReadValue("ref");
			if(pRef)
				m_TransitionState = dynamic_cast<TransitionState*>(m_pMoleculeManager->find(pRef));
			if(m_TransitionState)
			{
				const char* pthreshtxt = ppReac->ReadValue("me:threshold",false);
				if(pthreshtxt)
				{
					stringstream ss(pthreshtxt);
					ss >> m_E0;
					return true ;
				}
			}
		}
		return false ;
	}

	//
	// Calculate the forward microcanonical rate coefficients, using simple RRKM theory.
	//
	void SimpleRRKM::calculateMicroRateCoeffs(vector<double> &kfmc) {

		double plancksConst = 1.0/2.998e+10 ; // Planck's constant in wavenumbers.

		// Allocate space to hold Micro-canonical rate coefficients.
		kfmc.resize(pSys->MAXCell());

		// Initialize microcanoincal rate coefficients.

		int i, j ;
		for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
			kfmc[i] = 0.0 ;
		}

		// Allocate some work space for density of states.

		vector<double> ddosTS(pSys->MAXCell(),0.0) ; // Transistion state density of states.
		vector<double> ddos(pSys->MAXCell(),0.0) ; // Density of states of equilibrim molecule.

		// Extract densities of states from molecules.

		m_Reactant->cellDensityOfStates(&ddos[0]) ;
		m_TransitionState->cellDensityOfStates(&ddosTS[0]) ;

		double SumOfStates  = 0.0 ;
		int thresholdEnergy = int(m_E0 * KCMLTOPCM) ;
		for (i = thresholdEnergy, j = 0 ; i < pSys->MAXCell() ; ++i, ++j ) {

			// Integrate transition state density of states.

			SumOfStates += ddosTS[j] ;

			// Calculate microcanonical rate coefficients using RRKM expression.

			kfmc[i] = SumOfStates / (plancksConst*ddos[i]) ;
		}

		// Test microcanonical rate coefficients.

		//    if ( get_verbosity() ) 
		testMicroRateCoeffs(kfmc, m_Reactant) ;

	}

}//namespace
