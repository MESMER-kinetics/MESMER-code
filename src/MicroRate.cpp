#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants ;

namespace mesmer
{

  bool MicroRateCalculator::testMicroRateCoeffs(Reaction* pReact, PersistPtr ppbase) const
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
    size_t MaximumGrain = (pReact->getEnv().MaxGrn - pReact->get_fluxFirstNonZeroIdx());

    // Allocate some work space for density of states.

    vector<double> grainEne;
    vector<double> grainDOS;
    const vector<double>& grainKfmc = pReact->get_GrainKfmc();
    pReactant->getDOS().getGrainEnergies(grainEne) ;
    pReactant->getDOS().getGrainDensityOfStates(grainDOS) ;

    ctest << "\nCanonical (high pressure) rate coefficients for " << pReact->getName() << ", calculated from microcanonical rates\n{\n";
	ctest << pReact->TestRateCoeffHeader() << endl ;

	// Save the current value of excess concentration and set it to unity
	// to prevent division by zero for assocaiation type reactions.

	const double current_conc = pReact->get_concExcessReactant() ;
	pReact->set_concExcessReactant(1.0) ;

	// Save current value of beta.
	const double current_beta = pReact->getEnv().beta ;

    // Calculate Canonical rate coefficients up to the max. temperature givn by MesmerEnv.
	MesmerEnv &env = const_cast<MesmerEnv&>(pReact->getEnv()) ;
	double dTemp(100.0) ; // 100 K intervals.
    double Temp(0.0);
    size_t nTemp(size_t(pReact->getEnv().MaximumTemperature/dTemp)+1);
    for(size_t i(0) ; i < nTemp ; i++)
    {
      Temp += dTemp ;
      const double beta = 1.0/(boltzmann_RCpK*Temp) ;
      env.beta = beta ;

      double kf(0.0), qtot(0.0);
      for ( size_t i(0) ; i < MaximumGrain ; ++i ) {
        double tmp  = grainDOS[i] * exp(-beta * grainEne[i]) ;
        kf += grainKfmc[i] * tmp ;
        qtot += tmp ;
      }
      kf /= qtot ;
	  const double Keq = pReact->calcEquilibriumConstant() ;
	  const double kb  = kf/Keq ;
      formatFloat(ctest, Temp, 6,  7) ;
      formatFloat(ctest, kf,   6, 15) ;
      formatFloat(ctest, kb,   6, 15) ;
      formatFloat(ctest, Keq,  6, 15) ;
      ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:kinf");
      ppItem->XmlWriteValueElement("me:T",   Temp, 6) ;
      ppItem->XmlWriteValueElement("me:val", kf,   6) ;
      ppItem->XmlWriteValueElement("me:rev", kb,   6) ;
      ppItem->XmlWriteValueElement("me:Keq", Keq,  6) ;
    }
    ctest << "}\n";

	// Restore excess concentration value.
	pReact->set_concExcessReactant(current_conc) ;

	// Restore current value of beta.
	env.beta = current_beta ;

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
