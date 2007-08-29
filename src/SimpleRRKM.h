#ifndef GUARD_SimpleRRKM_h
#define GUARD_SimpleRRKM_h

#include <vector>
#include <string>
#include "System.h"

namespace mesmer
{

class SimpleRRKM : public MicroRateCalculator
{
public:

  ///Constructor which registers with the list of MicroRateCalculators in the base class
  SimpleRRKM(const std::string& id) : MicroRateCalculator(id){}

  virtual ~SimpleRRKM() {}

  virtual bool calculateMicroRateCoeffs(Reaction* pReact, std::vector<double> &kfmc);
};

}//namespace

#endif // GUARD_SimpleRRKM_h