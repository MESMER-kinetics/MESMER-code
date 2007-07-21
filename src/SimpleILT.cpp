//-------------------------------------------------------------------------------------------
//
// SimpleILT.cpp 
//
// Author: Struan Robertson 
// Date:   1/Jul/2007
//
// This file contains the implementation of the Simple ILT MC rate coefficient method.
//
//-------------------------------------------------------------------------------------------

#include <math.h>
#include "system.h"
#include "Constants.h"
#include "Persistence.h"
#include "SimpleILT.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{

    bool SimpleILT::Initialize(PersistPtr &ppReac) {

        // Are there ILT parameters available?

        const char* pActEnetxt = ppReac->ReadValue("me:activationenergy",false);
        const char* pPreExptxt = ppReac->ReadValue("me:preexponential",false);
        if (pActEnetxt && pPreExptxt) {

            stringstream s1(pActEnetxt);
            s1 >> m_ActEne ;
            stringstream s2(pPreExptxt);
            s2 >> m_PreExp ;

            return true ;
        } else {

            return false ;
        }

    }

    //
    // Calculate the forward microcanonical rate coefficients, using simple ILT theory.
    //
    void SimpleILT::calculateMicroRateCoeffs(vector<double> &kfmc) {

        // Allocate space to hold Micro-canonical rate coefficients.
        kfmc.resize(pSys->MAXCell());

        // Initialize microcanoincal rate coefficients.

        int i ;
        for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
            kfmc[i] = 0.0 ;
        }

        // Allocate some work space for density of states.

        vector<double> ddos(pSys->MAXCell(),0.0) ; // Density of states of equilibrim molecule.

        // Extract densities of states from molecules.

        m_Reactant->cellDensityOfStates(&ddos[0]) ;

        // Conversion of EINF from kcal.mol^-1 to cm^-1

		int nEinf = int(m_ActEne*KCMLTOPCM) ;

        // Calculate microcanonical rate coefficients using simple ILT expression.

		for (i = nEinf ; i < pSys->MAXCell() ; ++i ) {
            kfmc[i] = m_PreExp*ddos[i-nEinf] / ddos[i] ;
        }

        // Test microcanonical rate coefficients.

        //    if ( get_verbosity() ) 
        testMicroRateCoeffs(kfmc, m_Reactant) ;

    }

}//namespace
