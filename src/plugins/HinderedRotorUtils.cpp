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
#include "vector3.h"

namespace mesmer
{
  using namespace std;
  using OpenBabel::vector3;

  //
  // Calculate the reduced moment of inertia.
  //
  double HinderedRotorUtils::reducedMomentInertia(gStructure& gs, pair<string,string>& bondats, vector<double>& mode) {

    vector3 coords1 = gs.GetAtomCoords(bondats.first);
    vector3 coords2 = gs.GetAtomCoords(bondats.second);

    // Calculate moment of inertia about bond axis of atoms on one side of bond...
    vector<string> atomset;
    atomset.push_back(bondats.second); //will not look beyond this atom on the other side of the bond
    gs.GetAttachedAtoms(atomset, bondats.first);
    atomset.erase(atomset.begin()); //the other side of the bond is not in this set
    double mm1 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);
    gs.CalcInternalRotVec(atomset, coords1, coords2, mode) ;

    //...and the other side of the bond
    atomset.clear();
    atomset.push_back(bondats.first);
    gs.GetAttachedAtoms(atomset, bondats.second);
    atomset.erase(atomset.begin());
    double mm2 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);
    gs.CalcInternalRotVec(atomset, coords2, coords1, mode) ;

    /*
    Is the reduced moment of inertia needed about the bond axis or, separately for the set of
    atoms on each side of the bond, about a parallel axis through their centre of mass?
    See:
    http://www.ccl.net/chemistry/resources/messages/2001/03/21.005-dir/index.html
    http://www.ccl.net/chemistry/resources/messages/2001/03/31.002-dir/index.html
    The bond axis is used here.
    */

    return mm1 * mm2 / ( mm1 + mm2 ); //units a.u.*Angstrom*Angstrom

  }

  //
  // Calculate cosine coefficients from potential data points.
  //
  void HinderedRotorUtils::FourierCoeffs(vector<double> &angle, vector<double> &potential)
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

	// Determine the cosine coefficients.

	for(size_t k(0); k < m_expansion; ++k) {
	  double sum(0.0) ;
	  for(size_t i(0); i < ndata; ++i) {
		double nTheta = double(k) * angle[i];
		sum += potential[i] * cos(nTheta);
	  }
	  m_potentialCosCoeff.push_back(2.0*sum/double(ndata)) ;
	}
	m_potentialCosCoeff[0] /= 2.0 ;

	// Determine the sine coefficients.

	if (m_useSinTerms) {
	  for(size_t k(0); k < m_expansion; ++k) {
		double sum(0.0) ;
		for(size_t i(0); i < ndata; ++i) {
		  double nTheta = double(k) * angle[i];
		  sum += potential[i] * sin(nTheta);
		}
		m_potentialSinCoeff.push_back(2.0*sum/double(ndata)) ;
	  }
	  m_potentialSinCoeff[0] = 0.0 ;
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

