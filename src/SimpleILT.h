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

    virtual bool calculateMicroRateCoeffs(Reaction* pReact, std::vector<double> &cellKfmc, const MesmerEnv &mEnv);
      
    // provide a function to define particular counts of the convoluted DOS of two molecules
    bool countDimerCellDOS(ModelledMolecule* mol1, ModelledMolecule* mol2, std::vector<double> &dimerCellDOS, const MesmerEnv &mEnv);

  };
}//namespace

#endif // GUARD_SimpleILT_h
