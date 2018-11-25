#include <iomanip>
#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants;

namespace mesmer
{

  bool MicroRateCalculator::testMicroRateCoeffs(Reaction* pReact, PersistPtr ppbase) const
  {
    string comment("Canonical rate coefficients (calculated from microcanonical rate coefficients)");
    PersistPtr ppList = ppbase->XmlWriteMainElement("me:canonicalRateList", comment);

    stest << "\nCanonical (high pressure) rate coefficients for " << pReact->getName() << ", calculated from microcanonical rates\n{\n";
    //Number of reactants and products to set kf,kb and Keq units
    vector<Molecule*> vec;
    int nr = pReact->get_reactants(vec);
    int np = pReact->get_products(vec);
    stest << right << setw(7) << "T/K"
      << setw(20) << (nr == 2 ? "kf/cm3molecule-1s-1" : "kf/s-1")
      << setw(20) << (np == 2 ? "kb/cm3molecule-1s-1" : "kb/s-1")
      << setw(18) << ((np - nr) == 0 ? "Keq   " : ((np - nr) > 0 ? "Keq/moleculecm-3" : "Keq/cm3molecule-1"))
      << endl;

    // Save the current value of excess concentration and set it to unity
    // to prevent division by zero for assocaiation type reactions.
    const double current_conc = pReact->get_concExcessReactant();
    pReact->set_concExcessReactant(1.0);

    // Save current value of beta.
    const double current_beta = pReact->getEnv().beta;

    // Calculate Canonical rate coefficients up to the max. temperature givn by MesmerEnv.
    MesmerEnv &env = const_cast<MesmerEnv&>(pReact->getEnv());
    double dTemp = env.CRCTitl;
    double Temp  = env.CRCTMin;
    double TMax  = env.CRCTMax;
    vector<double> Coeffs;
    for (int j(0) ; Temp <= TMax; Temp += dTemp, j++) {
      env.beta = 1.0 / (boltzmann_RCpK*Temp);

      Coeffs.clear();
      pReact->HighPresRateCoeffs(&Coeffs);

      formatFloat(stest, Temp, 6, 7);
      formatFloat(stest, Coeffs[0], 6, 20);
      if (Coeffs.size()>1) // Output only forward rate if no ZPE has been provided.
      {
        formatFloat(stest, Coeffs[1], 6, 20);
        formatFloat(stest, Coeffs[2], 6, 18);
      }
      stest << endl;

      // Add to XML document.
      vector<Molecule*> vec;
      int nr = pReact->get_reactants(vec);
      int np = pReact->get_products(vec);
      PersistPtr ppItem = ppList->XmlWriteElement("me:kinf");
      PersistPtr pp = ppItem->XmlWriteValueElement("me:T", Temp, 6);
      if (j == 0) pp->XmlWriteAttribute("units", "K");
      pp = ppItem->XmlWriteValueElement("me:val", Coeffs[0], 6);
      if (j == 0) pp->XmlWriteAttribute("units", nr == 2 ? "cm3molecule-1s-1" : "s-1");
      if (Coeffs.size() > 1)
      {
        pp = ppItem->XmlWriteValueElement("me:rev", Coeffs[1], 6);
        if (j == 0) pp->XmlWriteAttribute("units", np == 2 ? "cm3molecule-1s-1" : "s-1");
        pp = ppItem->XmlWriteValueElement("me:Keq", Coeffs[2], 6);
        if (j == 0) pp->XmlWriteAttribute("units", ((np - nr) == 0 ? "" : ((np - nr) > 0 ? "moleculecm-3" : "cm3molecule-1")));
      }
    }
    stest << "}\n";

    // Restore excess concentration value.
    pReact->set_concExcessReactant(current_conc);

    // Restore current value of beta.
    env.beta = current_beta;

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
