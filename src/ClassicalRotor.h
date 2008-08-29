#ifndef GUARD_ClassicalRotor_h
#define GUARD_ClassicalRotor_h

#include "System.h"

namespace mesmer
{
  class ClassicalRotor : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the convolved DOS of two molecules.
    virtual bool countDimerCellDOS(ModelledMolecule* p_mol1, ModelledMolecule*  p_mol2, std::vector<double>& rctsCellDOS) ; 

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(ModelledMolecule* mol);

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    ClassicalRotor(const std::string& id) : DensityOfStatesCalculator(id){}

    virtual ~ClassicalRotor() {}

  } ;

}//namespace

#endif // GUARD_ClassicalRotor_h
