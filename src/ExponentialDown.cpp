
#include "ExponentialDown.h"
#include "MolecularComponents.h"
#include "Molecule.h"
#include <cmath>

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance).
  ExponentialDown exponentialDown("ExponentialDown");
  //************************************************************

  double ExponentialDown::calculateTransitionProbability(double Ei, double Ej) {
    // return exp(-(Ei -Ej)/m_deltaEdown) ;
    return exp(-m_alpha*(Ei -Ej)) ;
  }

  bool ExponentialDown::ReadParameters(const gWellProperties* pgWellProperties) { 

    const double refTemp = pgWellProperties->getDeltaEdownRefTemp();
    const double dEdExp = pgWellProperties->getDeltaEdownExponent();
    const double dEdRef = pgWellProperties->getDeltaEdown();
    const double temperature = 1.0/(boltzmann_RCpK * pgWellProperties->getHost()->getEnv().beta);

    m_deltaEdown = dEdRef * pow((temperature/refTemp),dEdExp);
    m_alpha      = 1.0/m_deltaEdown;

    // issue a warning message and exit if delta_E_down is smaller than grain size.
    if (m_deltaEdown < double(pgWellProperties->getHost()->getEnv().GrainSize) && !pgWellProperties->getHost()->getFlags().allowSmallerDEDown){
      cerr << "Delta E down is smaller than grain size: the solution may not converge.";
      return false;
    }
    
    return true ; 

  }

}//namespace

