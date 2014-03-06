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

	std::vector<double>& get_PotentialCosCoeff() {return m_potentialCosCoeff ; } ;
	std::vector<double>& get_PotentialSinCoeff() {return m_potentialSinCoeff ; } ;
	void set_Expansion(size_t expansion) {m_expansion = expansion ; } ;
	void set_UseSinTerms(bool useSinTerms) {m_useSinTerms = useSinTerms ; } ;
	size_t get_Expansion() const {return m_expansion ; } ;
	bool get_UseSinTerms() const {return m_useSinTerms ; } ;

	// Calculate the reduced moment of inertia.
	double reducedMomentInertia(gStructure& gs, pair<string,string>& bondats, std::vector<double>& mode) ;

	// Calculate the Fourier coefficients from potential data points.
	void FourierCoeffs(vector<double> &angle, vector<double> &potential) ;

	// Calculate potential.
	double CalculatePotential(double angle) const ;

	// Calculate gradient.
	double CalculateGradient(double angle) const ;

  private:

	std::vector<double> m_potentialCosCoeff ; // The cosine coefficients of the hindered rotor potential.
	std::vector<double> m_potentialSinCoeff ; // The sine coefficients of the hindered rotor potential.

	size_t m_expansion ;                 // Number of coefficients in the cosine expansion.

	bool m_useSinTerms ;                 // If true sine terms are used in the representation of the potential.

  } ;

}  //namespace

#endif // GUARD_HinderedRotorUtils_h
