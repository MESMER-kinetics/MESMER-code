#include "math.h"
#include "System.h"

//
// Test the forward microcanonical rate coefficients.
//
namespace mesmer
{
  using namespace std;
  using namespace Constants ;

  bool MicroRateCalculator::testMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc, PersistPtr ppbase) const
  {
    System* pSys = pReact->GetSys();  

    vector<CollidingMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    CollidingMolecule * pReactant = unimolecularspecies[0];

    cout << endl << "Test of microcanonical rate coefficients" << endl << endl ;
    string comment("Microcanonical rate coefficients");
    PersistPtr ppList = ppbase->WriteMainElement("me:microRateList", comment );

    // Allocate some work space for density of states.

    vector<double> decll(pSys->MAXCell(),0.0) ; 
    vector<double> ddos(pSys->MAXCell(),0.0) ; 

    pReactant->cellEnergies(&decll[0]) ;
    pReactant->cellDensityOfStates(&ddos[0]) ;

    for(int i = 0 ; i < 29 ; ++i)
    {
      double Temperature = double(i+2)*100.0 ;
      double beta = 1.0/(boltzmann*Temperature) ;

      double sm1 = 0.0 ;
      double sm2 = 0.0 ;
      double tmp = 0.0 ;
      for ( int i = 0 ; i < pSys->MAXCell() ; ++i ) {
        tmp  = ddos[i] * exp(-beta * decll[i]) ;
        sm1 += kfmc[i] * tmp ;
        sm2 +=           tmp ;
      }
      sm1 /= sm2 ; 
      formatFloat(cout, Temperature, 6,  7) ;
      formatFloat(cout, sm1,         6, 15) ;
      cout << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->WriteElement("me:microRate");
      ppItem->WriteValueElement("me:T",   Temperature, 6) ;
      ppItem->WriteValueElement("me:val", sm1,         6) ;
    }
    return true;
  }
}//namespace