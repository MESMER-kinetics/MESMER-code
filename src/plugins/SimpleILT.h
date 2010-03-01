#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

#include <vector>
#include <string>
#include "System.h"

namespace mesmer
{
  class SimpleILT : public MicroRateCalculator
  {
  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    SimpleILT(const std::string& id) : MicroRateCalculator(id){}

    virtual ~SimpleILT() {}
    virtual SimpleILT* Clone() { return new SimpleILT(*this); }

    virtual bool calculateMicroRateCoeffs(Reaction* pReact);

    virtual double get_ThresholdEnergy(Reaction* pReac) ;

  };
}//namespace

#endif // GUARD_SimpleILT_h