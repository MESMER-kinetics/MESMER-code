#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

#include "ClassicalRotor.h"

namespace mesmer
{
  class MesmerILT : public MicroRateCalculator
  {
  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    MesmerILT(const std::string& id) : MicroRateCalculator(id) {}

    virtual ~MesmerILT() {}

    virtual bool calculateMicroRateCoeffs(Reaction* pReact) ;

  private:

    bool calculateAssociationMicroRates(Reaction* pReact);
    bool calculateUnimolecularMicroRates(Reaction* pReact);

  };
}//namespace

#endif // GUARD_SimpleILT_h
