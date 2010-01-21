
#include "ExponentialDown.h"
#include "MolecularComponents.h"
#include "Molecule.h"

#include <cmath>

using namespace std ;

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id. A new instance (from Clone())is now used in the calculations.
  ExponentialDown exponentialDown("ExponentialDown");
  //************************************************************

  ExponentialDown* ExponentialDown::Clone() {
    return new ExponentialDown(*this); //probably deleted in gWellProperties destructor
  }

  double ExponentialDown::calculateTransitionProbability(double Ei, double Ej) {
    // return exp(-(Ei -Ej)/m_deltaEdown) ;

    double deltaEDown = m_deltaEdown;
    if(m_dEdExp!=0.0) {
      const double temperature = 1.0/(boltzmann_RCpK * getParent()->getEnv().beta);
      deltaEDown = deltaEDown * pow((temperature/m_refTemp),m_dEdExp);
    }

    // issue a warning message and exit if delta_E_down is smaller than grain size.
    if (deltaEDown < double(getParent()->getEnv().GrainSize) && !getParent()->getFlags().allowSmallerDEDown)
      cerr << "Delta E down is smaller than grain size: the solution may not converge.";

    return exp(-(Ei -Ej)/deltaEDown) ;
  }

  bool ExponentialDown::ReadParameters(const Molecule* parent) { 

    setParent(parent);
    PersistPtr pp = parent->get_PersistentPointer();
    PersistPtr ppProp = pp->XmlMoveToProperty("me:deltaEDown");

    const char* txt = (parent->get_PersistentPointer())->XmlReadProperty("me:deltaEDown"); //required
    if(!txt)
      return false;
    istringstream idata(txt);
    double value(0.0);
    idata >> value;
    m_deltaEdown = value;

    double valueL   = ppProp->XmlReadDouble("lower", optional);
    double valueU   = ppProp->XmlReadDouble("upper", optional);
    double stepsize = ppProp->XmlReadDouble("stepsize" ,optional);
    if (valueL!=0.0 && valueU!=0.0){
      // make a range variable
      m_deltaEdown.set_range(valueL, valueU, stepsize, "deltaEdown");//incl parent in name?
      //Save PersistPtr of the XML source of this Rdouble
      RangeXmlPtrs.push_back(ppProp);
    }

    // The temperature dependence of <delta_E_down> is accounted for as:
    //
    // <delta_E_down>(T) = <delta_E_down>_ref * (T / refTemp)^dEdExp
    //
    // By default, dEdExp = 0, which means delta_E_down does not depend on temperature.
    // Reference temperature of <Delta E down>, refTemp, has default 298.

    m_refTemp = ppProp->XmlReadDouble("referenceTemperature", optional );
    if(m_refTemp==0.0)
      m_refTemp = 298.;

    m_dEdExp = ppProp->XmlReadDouble("exponent", optional); //defaults to 0.0

    return true ; 
  }

}//namespace

