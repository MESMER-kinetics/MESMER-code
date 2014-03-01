#ifndef GUARD_HinderedRotorUtils_h
#define GUARD_HinderedRotorUtils_h

//-------------------------------------------------------------------------------------------
//
// HinderedRotorUtils.h
//
// Author: Struan Robertson
// Date:   1/Mar/2014
//
// Definition of a utility class that is inherited by the hindered rotor methods. 
//
//-------------------------------------------------------------------------------------------

#include "../MolecularComponents.h"
#include <vector>

namespace mesmer
{
  class gStructure ;

  class HinderedRotorUtils
  {
  public:

	HinderedRotorUtils() :  m_potentialCosCoeff(),
	  m_potentialSinCoeff(),
	  m_expansion(4),
	  m_useSinTerms(false) {} ; 
	~HinderedRotorUtils() {} ; 

  protected:

	std::vector<double> m_potentialCosCoeff ; // The cosine coefficients of the hindered rotor potential.
	std::vector<double> m_potentialSinCoeff ; // The sine coefficients of the hindered rotor potential.

	size_t m_expansion ;                 // Number of coefficients in the cosine expansion.

	bool m_useSinTerms ;                 // If true sine terms are used in the representation of the potential.

	// Calculate the reduced moment of inertia.
	double reducedMomentInertia(gStructure& gs, pair<string,string>& bondats, std::vector<double>& mode) ;

	// Calculate the Fourier coefficients from potential data points.
	void FourierCoeffs(vector<double> &angle, vector<double> &potential) ;

	// Calculate potential.
	double CalculatePotential(double angle) const ;

	// Calculate gradient.
	double CalculateGradient(double angle) const ;

  } ;

}  //namespace

#endif // GUARD_HinderedRotorUtils_h
