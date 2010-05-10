#ifndef GUARD_HinderedRotorA_h
#define GUARD_HinderedRotorA_h

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

namespace mesmer
{
  class HinderedRotorA : public DensityOfStatesCalculator
  {
  public:
    //Read data from XML. Some is stored hear and some in a MolecularComponent class.
    virtual bool ReadParameters(Molecule* pMol, PersistPtr ppDOSC=NULL);

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell);

    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    //needs to be specified.
    HinderedRotorA(const std::string& id) : DensityOfStatesCalculator(id, true){}

    virtual ~HinderedRotorA() {}
    virtual HinderedRotorA* Clone() { return new HinderedRotorA(*this); }

  private:
    std::string m_bondID;
    double      m_barrier;
    int         m_periodicity;
    double      m_vibFreq;
    double      m_reducedMomentInertia;
  } ;

}//namespace

#endif // GUARD_HinderedRotorA_h
