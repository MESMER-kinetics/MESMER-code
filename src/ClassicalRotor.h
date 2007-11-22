#ifndef GUARD_ClassicalRotor_h
#define GUARD_ClassicalRotor_h

#include "System.h"

namespace mesmer
{
  class ClassicalRotor : public DensityOfStatesCalculator
  {
    // provide a function to define particular counts of the convoluted DOS of two molecules
    bool countDimerCellDOS(SuperMolecule* rcts, const MesmerEnv &mEnv); // convolutes the DOS of two molecules

    // provide a function to define particular counts of the DOS of a molecule
    bool countMonomerCellDOS(ModelledMolecule* mol, const MesmerEnv &mEnv);

  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    ClassicalRotor(const std::string& id) : DensityOfStatesCalculator(id){}

    virtual ~ClassicalRotor() {}

    virtual bool countCellDOS(ModelledMolecule* mol, const MesmerEnv &mEnv);

  };
}//namespace

#endif // GUARD_ClassicalRotor_h
