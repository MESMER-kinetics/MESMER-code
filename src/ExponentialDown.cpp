
#include "ExponentialDown.h"
#include "MolecularComponents.h"
#include "Molecule.h"

#include <cmath>

using namespace std ;

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

    PersistPtr ppPropList = pgWellProperties->get_MolProp() ;

    // The temperature dependence of <delta_E_down> is accounted for as:
    //
    // <delta_E_down>(T) = <delta_E_down>_ref * (T / refTemp)^dEdExp
    //
    // By default, dEdExp = 0, which means delta_E_down does not depend on temperature.
    // Reference temperature of <Delta E down>, refTemp, has default 298.

    double refTemp(298.);
    const char* pRefTemptxt  = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "referenceTemperature", optional );
    if(pRefTemptxt){
      stringstream s_temp(pRefTemptxt); s_temp >> refTemp;
    }

    double dEdExp(0.0);
    const char* pExponenttxt = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "exponent", optional);
    if(pExponenttxt){
      stringstream s_exp(pExponenttxt); s_exp >> dEdExp;
    }    

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

