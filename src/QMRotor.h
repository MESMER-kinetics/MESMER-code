#ifndef GUARD_QMRotor_h
#define GUARD_QMRotor_h

#include "DensityOfStates.h"
#include "MolecularComponents.h"

namespace mesmer
{
  class QMRotor : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell, PersistPtr ppDOSC);

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    QMRotor(const std::string& id) : DensityOfStatesCalculator(id){}

    virtual ~QMRotor() {}

  } ;

}//namespace

#endif // GUARD_QMRotor_h
