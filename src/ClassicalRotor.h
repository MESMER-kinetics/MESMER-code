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
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell, PersistPtr ppDOSC=NULL);

    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class calculates a complete DOS: it is not an extra class. 
    ClassicalRotor(const std::string& id) : DensityOfStatesCalculator(id, false){}

    virtual ~ClassicalRotor() {}

  } ;

}//namespace

#endif // GUARD_ClassicalRotor_h
