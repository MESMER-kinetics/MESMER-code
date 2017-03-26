#ifndef GUARD_Potential1D_h
#define GUARD_Potential1D_h

//-------------------------------------------------------------------------------------------
//
// Potential1D.h
//
// Author: Struan Robertson
// Date:   25/Mar/2017
//
// Definition of a the Potential1D class. The principal use of this class is to express the 
// potential of a single anharmonic mode such that its states can be calculated. There are
// two types of potential: an analytic potential expressed as a polynomial and a numerical
// potential that is represented as a spline. These potentials expose the the same abstract
// base class.
//
//-------------------------------------------------------------------------------------------

#include "XMLPersist.h"
#include <vector>
#include "../Spline.h"

namespace mesmer
{

  // Abstract base class.

  class Potential1D
  {
  public:

		Potential1D(size_t npoints) : m_units("kJ/mol"), m_npoints(101), m_minx(0.0), m_maxx(10.0) { m_npoints = npoints; };
		virtual ~Potential1D() {};

    virtual void InitializePotential(PersistPtr pp) = 0;

    virtual void calculatePotential(std::vector<double> &abscissa, std::vector<double> &potential) const = 0;

		double get_characteristicLength() const { return (m_maxx - m_minx) ; }

  protected:

    std::string m_units; // Units that potential energy is expressed in.
    size_t m_npoints;    // Number of absissca points for which the potential is required.
    double m_minx;       // The maximum value of the abscissa to use.
    double m_maxx;       // The maximum value of the abscissa to use.
  };

  // Analytical potential based on a polynomial.

  class AnalyticalPotential : public Potential1D
  {
  public:

    AnalyticalPotential(size_t npoints) : Potential1D(npoints), m_Coeff() {};
    ~AnalyticalPotential() {};

    virtual void InitializePotential(PersistPtr pp);

    virtual void calculatePotential(std::vector<double> &abscissa, std::vector<double> &potential) const;

  private:

    std::vector<double> m_Coeff; // The polynomial coefficients of the vibration potential.
  };

  // Numerical potential based on a spline.

  class NumericalPotential : public Potential1D
  {
  public:

    NumericalPotential(size_t npoints) : Potential1D(npoints) {};
    ~NumericalPotential() {};

    virtual void InitializePotential(PersistPtr pp);

    virtual void calculatePotential(std::vector<double> &abscissa, std::vector<double> &potential) const;

  private:

		Spline m_spline;
	
	};

}  //namespace

#endif // GUARD_Potential1D_h
