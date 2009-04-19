#ifndef GUARD_HinderedRotorInterpolation_h
#define GUARD_HinderedRotorInterpolation_h

#include "DensityOfStates.h"
#include "MolecularComponents.h"

namespace mesmer
{
  class HinderedRotorInterpolation : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell, PersistPtr ppDOSC);

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    HinderedRotorInterpolation(const std::string& id) : DensityOfStatesCalculator(id){}

    virtual ~HinderedRotorInterpolation() {}

  } ;

}//namespace

#endif // GUARD_HinderedRotorInterpolation_h
