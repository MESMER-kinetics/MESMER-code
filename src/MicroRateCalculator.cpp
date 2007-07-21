//-------------------------------------------------------------------------------------------
//
// MicroRateCalculator.cpp
//
// Author: Struan Robertson 
// Date:   21/Jul/2007
//
// This file contains the implementation of the MicroRateCalculator class.
//
//-------------------------------------------------------------------------------------------

#include <math.h>
#include "system.h"
#include "Constants.h"
#include "Persistence.h"
#include "MicroRateCalculator.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
    //
    // Test the forward microcanonical rate coefficients.
    //
    void MicroRateCalculator::testMicroRateCoeffs(vector<double> &kfmc, CollidingMolecule *m_Reactant) const {

        cout << endl << "Test of microcanonical rate coefficients" << endl << endl ;
        string comment("Microcanonical rate coefficients");
        //		PersistPtr ppList = m_ppPersist->WriteMainElement("me:microRateList", comment );

        // Allocate some work space for density of states.

        vector<double> decll(pSys->MAXCell(),0.0) ; 
        vector<double> ddos(pSys->MAXCell(),0.0) ; 

        m_Reactant->cellEnergies(&decll[0]) ;
        m_Reactant->cellDensityOfStates(&ddos[0]) ;

        for ( int n = 0 ; n < 29 ; ++n ) {

            double temp = 100.0*static_cast<double>(n + 2) ;
            double beta = 1.0/(0.695029*temp) ;

            double sm1 = 0.0 ;
            double sm2 = 0.0 ;
            double tmp = 0.0 ;
            for ( int i = 0 ; i < pSys->MAXCell() ; ++i ) {
                tmp  = ddos[i] * exp(-beta * decll[i]) ;
                sm1 += kfmc[i] * tmp ;
                sm2 +=           tmp ;
            }
            sm1 /= sm2 ; 
            formatFloat(cout, temp, 6, 7) ;
            formatFloat(cout, sm1,  6, 15) ;
            cout << endl ;

            //Add to XML document
            //	PersistPtr ppItem = ppList->WriteElement("me:microRate");
            //	ppItem->WriteValueElement("me:T",   temp, 6);
            //	ppItem->WriteValueElement("me:val", sm1,  6) ;

        }
    }

}//namespace
