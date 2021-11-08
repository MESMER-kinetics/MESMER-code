#include <iomanip>
#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants;

namespace mesmer
{

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
        << " uses ILT method, which should not have transition state." << endl;
      return false;
    }
    const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling", optional);
    if (pTunnelingtxt)
    {
      cerr << "Tunneling parameter in Reaction " << pReac->getName() << " is invalid in ILT." << endl;
      return false;
    }

    const char* pCrossingtxt = ppReac->XmlReadValue("me:crossing", optional);
    if (pCrossingtxt)
    {
      cerr << "Crossing parameter in Reaction " << pReac->getName() << " is invalid in ILT." << endl;
      return false;
    }

    return true;

  }

}//namespace
