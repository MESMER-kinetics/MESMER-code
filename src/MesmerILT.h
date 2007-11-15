#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

#include "System.h"

namespace mesmer
{
  class MesmerILT : public MicroRateCalculator
  {
  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    MesmerILT(const std::string& id) : MicroRateCalculator(id){}

    virtual ~MesmerILT() {}

    virtual bool calculateMicroRateCoeffs(Reaction* pReact, std::vector<double> &cellKfmc, const MesmerEnv &mEnv);

    // provide a function to define particular counts of the convoluted DOS of two molecules
    virtual bool countDimerCellDOS(SuperMolecule* rcts, const MesmerEnv &mEnv);

  };
}//namespace

#endif // GUARD_SimpleILT_h
