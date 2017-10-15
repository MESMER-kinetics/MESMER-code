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

	// Read potential parameters
	void HinderedRotorUtils::ReadPotentialParameters(PersistPtr ppHR, const string& bondID, vector<double>& potentialCosCoeff, vector<double>& potentialSinCoeff) {

		PersistPtr pp = ppHR->XmlMoveTo("me:HinderedRotorPotential");

		if (pp) {

			const char* p = pp->XmlReadValue("format", true);
			string format(p);

			p = pp->XmlReadValue("units", optional);
			string units = p ? p : "kJ/mol";

			if (format == "analytical") {

				// Analytical potential.

				vector<int> indicies;
				vector<double> coefficients;
				int maxIndex(0);
				while (pp = pp->XmlMoveTo("me:PotentialPoint"))
				{
					int index = pp->XmlReadInteger("index", optional);
					indicies.push_back(index);
					maxIndex = max(maxIndex, index);

					double coefficient = pp->XmlReadDouble("coefficient", optional);
					if (IsNan(coefficient))
						coefficient = 0.0;
					coefficient = getConvertedEnergy(units, coefficient);
					coefficients.push_back(coefficient);
				}

				// As coefficients can be supplied in any order, they are sorted here.
				potentialCosCoeff.resize(++maxIndex);
				potentialSinCoeff.resize(++maxIndex);
				for (size_t i(0); i < coefficients.size(); i++) {
					potentialCosCoeff[indicies[i]] = coefficients[i];
					potentialSinCoeff[i] = 0.0;
				}

				// Shift potential so that lowest minimum is at zero.

				ShiftPotential();

			}
			else if (format == "numerical") {

				// Numerical potential.

				vector<double> potential;
				vector<double> angle;
				set_Expansion(pp->XmlReadInteger("expansionSize", optional));

				// Check if sine terms are to be used.
				set_UseSinTerms(pp->XmlReadBoolean("useSineTerms") || pp->XmlReadBoolean("UseSineTerms"));

				while (pp = pp->XmlMoveTo("me:PotentialPoint"))
				{
					double anglePoint = pp->XmlReadDouble("angle", optional);
					if (IsNan(anglePoint))
						anglePoint = 0.0;
					angle.push_back(anglePoint);

					double potentialPoint = pp->XmlReadDouble("potential", optional);
					if (IsNan(potentialPoint))
						potentialPoint = 0.0;
					potentialPoint = getConvertedEnergy(units, potentialPoint);
					potential.push_back(potentialPoint);
				}

				PotentialFourierCoeffs(angle, potential);

			}
			else {

				// Unknown format.

				cinfo << "Unknown hindering potential format for " << bondID << ", assuming free rotor." << endl;
				potentialCosCoeff.push_back(0.0);

			}

		}

	}

	// Shift potential to origin.
	void HinderedRotorUtils::ShiftPotential() {

		// A coarse search for minima is done over intervals of 1 degree a minima 
		// being located when the gradient changes sign from -ve to +ve. Then a 
		// Newton-Raphson type iteration is applied to get a better estimate of
		// the location of the minimum. Finally, the potential is shifted.

		double minPotential(1.e10);
		double anga(0.0), grda(0.0);
		size_t nDegrees(360);
		for (size_t i(0); i <= nDegrees; ++i) {
			double angb = double(i) * M_PI / 180.;
			double grdb = CalculateGradient(angb);
			if (grdb < 0.0) {
				anga = angb;
				grda = grdb;
			}
			else if (grda < 0.0) {
				double angc(0.0), grdc(1.0), tol(1.e-05);
				int n(0);
				while (fabs(grdc) > tol && n < 10) {
					angc = (grdb*anga - grda*angb) / (grdb - grda);
					grdc = CalculateGradient(angc);
					if (grdc > 0.0) {
						anga = angc;
						grda = grdc;
					}
					else {
						angb = angc;
						grdb = grdc;
					}
					n++;
				}
				double potential = CalculatePotential(angc);
				minPotential = min(minPotential, potential);
			}
		}
	}

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

	// Update the potential and and configuration phase difference.

	m_phase += amin ;

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

