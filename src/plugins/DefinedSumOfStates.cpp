//-------------------------------------------------------------------------------------------
//
// DefinedSumOfStates.cpp
//
// Author: Struan Robertson
// Date:   24/Mar/2012
//
// This class implements the defined sum of states model. This class exists to allow users to
// read in their own transition states sum of states, e.g. from FTST calculations or similar.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class DefinedSumOfStates : public MicroRateCalculator
  {
  public:

    // Constructor which registers with the list of MicroRateCalculators in the base class.
    DefinedSumOfStates(const std::string& id) : 
        MicroRateCalculator(id) {}

        virtual ~DefinedSumOfStates() {}
        virtual DefinedSumOfStates* Clone() { return new DefinedSumOfStates(*this); }

        virtual bool calculateMicroRateCoeffs(Reaction* pReac);

        virtual bool ReadParameters(Reaction* pReac) ;

  private:


  };

  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  DefinedSumOfStates theDefinedSumOfStates("DefinedSumOfStates");
  //************************************************************

  bool DefinedSumOfStates::ReadParameters(Reaction* pReact) {

    return true ; 
  }

  //
  // This method calculates the reaction flux. 
  //

  bool DefinedSumOfStates::calculateMicroRateCoeffs(Reaction* pReact)
  {

    PersistPtr ppReac = pReact->get_PersistentPointer();

    PersistPtr pp = ppReac->XmlMoveTo("me:SumOfStates") ;

    if (!pp) {
      throw(std::runtime_error("No sum of states define for reaction"+pReact->getName()))  ;
    }

    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "kJ/mol";

    while (pp = pp->XmlMoveTo("me:SumOfStatesPoint")) {
      double energy = pp->XmlReadDouble("energy", true);
      energy = getConvertedEnergy(units, energy);
      double angMomMag = pp->XmlReadDouble("angMomMag", true);
      double SumOfStates = pp->XmlReadDouble("angMomMag", true);
    }

    return true;
  }

}//namespace
