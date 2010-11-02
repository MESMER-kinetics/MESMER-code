#ifndef GUARD_HinderedRotorQM1D_h
#define GUARD_HinderedRotorQM1D_h

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

namespace mesmer
{
  class HinderedRotorQM1D : public DensityOfStatesCalculator
  {
  public:
    //Read data from XML. Some is stored hear and some in a MolecularComponent class.
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell);

    // Provide a function to calculate contribution to canonical partition function.
    // (Mostly for testing purposes.)
    virtual double canPrtnFnCntrb(const double beta) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorQM1D(const std::string& id) : DensityOfStatesCalculator(id, true),
      m_bondID(),
      m_reducedMomentInertia(0.0),
      m_periodicity(1),
      m_potentialCosCoeff(),
      m_expansion(4),
	  m_energyLevels() {}

    virtual ~HinderedRotorQM1D() {}
    virtual HinderedRotorQM1D* Clone() { return new HinderedRotorQM1D(*this); }

  private:

    // Calculate cosine coefficients from potential data points.
    void CosineFourierCoeffs(vector<double> &angle, vector<double> &potential) ;

    std::string m_bondID;

    double m_reducedMomentInertia;
    int    m_periodicity;

	vector<double> m_potentialCosCoeff ; // The cosine coefficients of the hindered rotor potential.

    size_t m_expansion ;                 // Number of coefficients in the cosine expansion.

	vector<double> m_energyLevels ;	     // The energies of the hindered rotor states.

  } ;

}//namespace

#endif // GUARD_HinderedRotorQM1D_h
