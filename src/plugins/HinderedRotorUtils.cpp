//-------------------------------------------------------------------------------------------
//
// HinderedRotorUtils.cpp
//
// Author: Struan Robertson
// Date:   1/Mar/2014
//
// Implementation of a utility class that is inherited by the hindered rotor methods. 
//
//-------------------------------------------------------------------------------------------

#include "HinderedRotorUtils.h"
#include "../vector3.h"

namespace mesmer
{
  using namespace std;
  using OpenBabel::vector3;

  //
  // Calculate cosine coefficients from potential data points.
  //
  void HinderedRotorUtils::PotentialFourierCoeffs(vector<double> &angle, vector<double> &potential)
  {
	size_t ndata = potential.size() ;

	// Locate the potential minimum and shift to that minimum.

	double vmin(potential[0]), amin(angle[0]) ;
	for (size_t i(1); i < ndata; ++i) {
	  if (potential[i] < vmin){
		vmin = potential[i] ;
		amin = angle[i] ;
	  }
	}

	for (size_t i(0); i < ndata; ++i) {
	  potential[i] -= vmin ;
	  angle[i]     -= amin ;
	  angle[i]     *= M_PI/180. ;
	}

    FourierCosCoeffs(angle, potential, m_potentialCosCoeff, m_expansion) ;
	if (m_useSinTerms) {
      FourierSinCoeffs(angle, potential, m_potentialSinCoeff, m_expansion) ;
	} else {
	  for(size_t k(0); k < m_expansion; ++k) {
		m_potentialSinCoeff.push_back(0.0) ;
	  }
	}

	// Test potential

	ctest << "          Angle         Potential          Series\n";
	for (size_t i(0); i < ndata; ++i) {
	  double clcPtnl = CalculatePotential(angle[i]) ;
	  ctest << formatFloat(angle[i], 6, 15) << ", " <<  formatFloat(potential[i], 6, 15) << ", " <<  formatFloat(clcPtnl, 6, 15) <<'\n' ;
	}
	ctest << endl ;

	return ;
  }

  // Calculate potential.
  double HinderedRotorUtils::CalculatePotential(double angle) const {

	if (m_potentialCosCoeff.size() == 0)
	  return 0.0 ;

	double sum(0.0) ;
	for(size_t k(0); k < m_potentialCosCoeff.size(); ++k) {
	  double nTheta = double(k) * angle;
	  sum += m_potentialCosCoeff[k] * cos(nTheta) + m_potentialSinCoeff[k] * sin(nTheta);
	}

	return sum ;
  }

  // Calculate potential gradient.
  double HinderedRotorUtils::CalculateGradient(double angle) const {

	if (m_potentialCosCoeff.size() == 0)
	  return 0.0 ;

	double sum(0.0) ;
	for(size_t k(0); k < m_potentialCosCoeff.size(); ++k) {
	  double nTheta = double(k) * angle;
	  sum += double(k)*(-m_potentialCosCoeff[k]*sin(nTheta) + m_potentialSinCoeff[k]*cos(nTheta)) ;
	}

	return sum ;
  }

}//namespace

