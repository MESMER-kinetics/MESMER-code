#ifndef GUARD_HinderedRotorA_h
#define GUARD_HinderedRotorA_h

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

namespace mesmer
{
  class HinderedRotorA : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell);

    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    //needs to be specified.
    HinderedRotorA(const std::string& id) : DensityOfStatesCalculator(id, true){}

    virtual ~HinderedRotorA() {}
    virtual HinderedRotorA* Clone() { return new HinderedRotorA(*this); }

    //Read data from XML and store in this instance.
    virtual bool ReadParameters(Molecule* pMol, PersistPtr ppDOSC=NULL);

  private:
    std::string m_bondID;
    double      m_barrier;
    int         m_periodicity;
    double      m_vibFreq;
  } ;

}//namespace

#endif // GUARD_HinderedRotorA_h
