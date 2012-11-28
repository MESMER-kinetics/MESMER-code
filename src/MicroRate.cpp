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
    if (!unimolecularspecies.size()){
      ctest << "\nNo microcanonical rate coefficients for " << pReact->getName() << endl;
      return true;
    }
    Molecule * pReactant = unimolecularspecies[0];

    string comment("Microcanonical rate coefficients");
    PersistPtr ppList = ppbase->XmlWriteMainElement("me:microRateList", comment );
    int MaximumGrain = (pReact->getEnv().MaxGrn - pReact->get_fluxFirstNonZeroIdx());

    // Allocate some work space for density of states.

    vector<double> grainEne;
    vector<double> grainDOS;
    const vector<double>& grainKfmc = pReact->get_GrainKfmc();
    pReactant->getDOS().getGrainEnergies(grainEne) ;
    pReactant->getDOS().getGrainDensityOfStates(grainDOS) ;

    ctest << "\nCanonical rate coefficients for " << pReact->getName() << ", calculated from microcanonical rates\n{\n";
    
		// get the Temperature from the MemserEnv beta value
		double MaxTemp = 1.0e+0/((pReact->getEnv().beta)*boltzmann_RCpK);
		double Temperature(0.0e+0);
		int num_temps(0);

		// calculate Canonical rate coefficients up to the Temperature of a particular MesmerEnv
		while(Temperature <= MaxTemp)
    {
      Temperature = double(num_temps+1)*100.0 ;
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

			++num_temps;
    }
    ctest << "}\n";

    return true;
  }
  
  //
  // This function retrieves the activation/threshold energy for an association reaction.
  //
  double MicroRateCalculator::get_ThresholdEnergy(Reaction* pReac) {

    if (!pReac->get_TransitionState()) {
      string s("No Transition State for ");
      throw (std::runtime_error(s + getID())); 
    }

    return (pReac->get_relative_TSZPE() - pReac->get_relative_rctZPE());
  }

//-----------------------------------------------------------------------------------------------
//
// ILT Utility methods
//

  //
  // Utility function to check for inconsistencies. 
  //
  bool MicroRateCalculator::ILTCheck(Reaction* pReac, PersistPtr ppReac)
  {
    // A few checks on features not allowed in ILT methods.
    
    if (pReac->get_TransitionState())
    {
      cerr << "Reaction " << pReac->getName() 
        << " uses ILT method, which should not have transition state."<<endl;
      return false;
    }
    const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling", optional) ;
    if(pTunnelingtxt)
    {
      cerr << "Tunneling parameter in Reaction " << pReac->getName() << " is invalid in ILT."<<endl;
      return false;
    }

    const char* pCrossingtxt = ppReac->XmlReadValue("me:crossing", optional) ;
    if(pCrossingtxt)
    {
      cerr << "Crossing parameter in Reaction " << pReac->getName() << " is invalid in ILT."<<endl;
      return false;
    }
    
    return true ;

  }

}//namespace
