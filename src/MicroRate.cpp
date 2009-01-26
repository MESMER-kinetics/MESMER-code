#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants ;

namespace mesmer
{

  bool MicroRateCalculator::testMicroRateCoeffs(
    Reaction*  pReact,
    PersistPtr ppbase) const
  {
    vector<Molecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    Molecule * pReactant = unimolecularspecies[0];

    string comment("Microcanonical rate coefficients");
    PersistPtr ppList = ppbase->XmlWriteMainElement("me:microRateList", comment );
    int MaximumGrain = (pReact->getEnv().MaxGrn - pReact->get_fluxFirstNonZeroIdx());

    // Allocate some work space for density of states.

    vector<double> grainEne;
    vector<double> grainDOS;
    const vector<double>& grainKfmc = pReact->get_GrainKfmc();
    pReactant->g_dos->getGrainEnergies(grainEne) ;
    pReactant->g_dos->getGrainDensityOfStates(grainDOS) ;

    ctest << "\nCanonical rate coefficients for " << pReact->getName() << ", calculated from microcanonical rates\n{\n";
    for(int i = 0 ; i < 29 ; ++i)
    {
      double Temperature = double(i+2)*100.0 ;
      double beta = 1.0/(boltzmann_RCpK*Temperature) ;

      double sm1 = 0.0, sm2 = 0.0, tmp = 0.0;

      for ( int i = 0 ; i < MaximumGrain ; ++i ) {
        tmp  = grainDOS[i] * exp(-beta * grainEne[i]) ;
        sm1 += grainKfmc[i] * tmp ;
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
  
  std::string MicroRateCalculator::getName() {return name;}

}//namespace
