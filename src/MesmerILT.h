#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

#include <vector>
#include <string>
#include "System.h"

namespace mesmer
{
  class MesmerILT : public MicroRateCalculator
  {
  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    MesmerILT(const std::string& id) : MicroRateCalculator(id){}

    virtual ~MesmerILT() {}

    virtual bool calculateMicroRateCoeffs(Reaction* pReact, std::vector<double> &kfmc);
  };
}//namespace

#endif // GUARD_SimpleILT_h
