#ifndef GUARD_ClassicalRotor_h
#define GUARD_ClassicalRotor_h

#include "DensityOfStates.h"
#include "MolecularComponents.h"

namespace mesmer
{
  class ClassicalRotor : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell, PersistPtr ppDOSC);

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    ClassicalRotor(const std::string& id) : DensityOfStatesCalculator(id){}

    virtual ~ClassicalRotor() {}

  } ;

}//namespace

#endif // GUARD_ClassicalRotor_h
