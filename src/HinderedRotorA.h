#ifndef GUARD_HinderedRotorA_h
#define GUARD_HinderedRotorA_h

#include "DensityOfStates.h"
#include "MolecularComponents.h"

namespace mesmer
{
  class HinderedRotorA : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell,
      PersistPtr ppDOSC=NULL);

    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class is an extra DOS class: there needs to be a non-extra DensityOfStatesCalculator class
    HinderedRotorA(const std::string& id) : DensityOfStatesCalculator(id, true){}

    virtual ~HinderedRotorA() {}

  } ;

}//namespace

#endif // GUARD_HinderedRotorA_h
