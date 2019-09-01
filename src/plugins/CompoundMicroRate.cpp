//-------------------------------------------------------------------------------------------
//
// CompoundMicroRate.cpp
//
// Author: Struan Robertson
// Date:    Aug/2019
//
// Calculates microcanonical rate coefficients as the sum of other microcanonical rates.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"
#include "../gDensityOfStates.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  class CompoundMicroRate : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    CompoundMicroRate(const char* id) : m_id(id) { Register(); }

    virtual ~CompoundMicroRate() {}
    virtual const char* getID()  { return m_id; }

    virtual CompoundMicroRate* Clone() { return new CompoundMicroRate(*this); }

    virtual bool ParseData(PersistPtr);

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

  private:
    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  CompoundMicroRate theCompoundMicroRate("CompoundMicroRate");
  //************************************************************

  bool CompoundMicroRate::calculateMicroCnlFlux(Reaction* pReact)
  {
    return true;
  }

  bool CompoundMicroRate::ParseData(PersistPtr)
  {
    bool status(true);

    return status;
  }

}//namespace
