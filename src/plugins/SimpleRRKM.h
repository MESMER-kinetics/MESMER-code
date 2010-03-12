#ifndef GUARD_SimpleRRKM_h
#define GUARD_SimpleRRKM_h

#include <vector>
#include <string>
#include "../System.h"

namespace mesmer
{
  class SimpleRRKM : public MicroRateCalculator
  {
  public:
  
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    SimpleRRKM(const std::string& id) : MicroRateCalculator(id){}
  
    virtual ~SimpleRRKM() {}

    virtual SimpleRRKM* Clone() { return new SimpleRRKM(*this); }

    virtual bool ReadParameters(Reaction* pReac);
  
    virtual bool calculateMicroRateCoeffs(Reaction* pReact);
  };
}//namespace

#endif // GUARD_SimpleRRKM_h
