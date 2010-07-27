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

    // Provide a function to calculate contribution to canonical partition function.
    // (Mostly for testing purposes.)
    virtual double canPrtnFnCntrb(const double beta) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorA(const std::string& id) : DensityOfStatesCalculator(id, true),
      m_bondID(),
      m_reducedMomentInertia(0.0),
      m_periodicity(1),
      m_potentialCosCoeff(),
      m_expansion(4) {}

    virtual ~HinderedRotorA() {}
    virtual HinderedRotorA* Clone() { return new HinderedRotorA(*this); }

  private:

    // Calculation of the modifed Bessel function, Io(x), for real x.
    double ModifiedBessalFuncion(const double x) const ;

    // Calculate cosine coefficients from potential data points.
    void CosineFourierCoeffs(vector<double> &angle, vector<double> &potential) ;

    std::string m_bondID;

    double m_reducedMomentInertia;
    int    m_periodicity;

	vector<double> m_potentialCosCoeff ;	// The cosine coefficients of the hindered rotor potential.

    size_t m_expansion ; // Number of coefficients in the cosine expansion.

  } ;

}//namespace

#endif // GUARD_HinderedRotorA_h
