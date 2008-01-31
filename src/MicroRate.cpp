#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants ;

namespace mesmer
{

  bool MicroRateCalculator::testMicroRateCoeffs(Reaction*         pReact,
                                                PersistPtr        ppbase) const
  {
    vector<ModelledMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    ModelledMolecule * pReactant = unimolecularspecies[0];

    cinfo << "Test of microcanonical rate coefficients" << endl;

    string comment("Microcanonical rate coefficients");
    PersistPtr ppList = ppbase->XmlWriteMainElement("me:microRateList", comment );
    int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate some work space for density of states.

    vector<double> cellEne;
    vector<double> cellDOS;

    pReactant->getCellEnergies(cellEne) ;
    pReactant->getCellDensityOfStates(cellDOS) ;

    ctest << "\nMicrocanonical rate coefficients:\n{\n";
    for(int i = 0 ; i < 29 ; ++i)
    {
      double Temperature = double(i+2)*100.0 ;
      double beta = 1.0/(boltzmann_RCpK*Temperature) ;

      double sm1 = 0.0, sm2 = 0.0, tmp = 0.0;

      for ( int i = 0 ; i < MaximumCell ; ++i ) {
        tmp  = cellDOS[i] * exp(-beta * cellEne[i]) ;
        sm1 += pReact->m_CellKfmc[i] * tmp ;
        sm2 += tmp ;
      }
      sm1 /= sm2 ;
      formatFloat(ctest, Temperature, 6,  7) ;
      formatFloat(ctest, sm1,         6, 15) ;
      ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:microRate");
      ppItem->XmlWriteValueElement("me:T",   Temperature, 6) ;
      ppItem->XmlWriteValueElement("me:val", sm1,         6) ;
    }
    ctest << "}\n";
        
    return true;
  }
}//namespace
