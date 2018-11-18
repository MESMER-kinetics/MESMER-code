#ifndef GUARD_MesmerMath_h
#define GUARD_MesmerMath_h

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "Matrix.h"
#include "MesmerPrecision.h"

using namespace mesmer;

//
// Some basic utility functions:
//
inline int nint(double x) { return (x>0)? int(x+0.5):int(x-0.5) ; } ;

template<class T>
inline const T SQR(const T a) {return a*a;}

// This routine is copied from the following source and modified for purpose to used as a template:
//  ggm.cpp -- computation of ggm function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
// Returns ggm function of argument 'x'.
//
// NOTE: Returns 1e308 if argument is a negative integer or 0,
//      or if argument exceeds 171.
//

template <class T>
const T MesmerGamma(const T& x)
{
  T T_PI = acos(-1.);
  int i,k,m;
  T ga(0.0),gr(0.0),r(0.0),z(0.0);

  static T g[] = {
	1.0,
	0.5772156649015329,
	-0.6558780715202538,
	-0.420026350340952e-1,
	0.1665386113822915,
	-0.421977345555443e-1,
	-0.9621971527877e-2,
	0.7218943246663e-2,
	-0.11651675918591e-2,
	-0.2152416741149e-3,
	0.1280502823882e-3,
	-0.201348547807e-4,
	-0.12504934821e-5,
	0.1133027232e-5,
	-0.2056338417e-6,
	0.6116095e-8,
	0.50020075e-8,
	-0.11812746e-8,
	0.1043427e-9,
	0.77823e-11,
	-0.36968e-11,
	0.51e-12,
	-0.206e-13,
	-0.54e-14,
	0.14e-14
  };

  if (x > 171.0) 
	return 1e308;    // This value is an overflow flag.

  if ((to_double(x) - double(int(to_double(x)))) < 1.e-14) {
	if (x > 0.0) {
	  ga = 1.0;               // use factorial
	  for (i=2; i<x; i++) 
		ga *= i;
	}
	else 
	  ga = 1e308;
  } else {
	if (abs(x) > 1.0) {
	  z = abs(x);
	  m = (int)to_double(z);
	  r = 1.0;
	  for (k=1;k<=m;++k) 
		r *= (z-k);
	  z -= m;
	}
	else 
	  z = x;

	gr = g[24];
	for (k=23;k>=0;--k) 
	  gr = gr*z+g[k];

	ga = 1.0/(gr*z);
	if (abs(x) > 1.0) {
	  ga *= r;
	  if (x < 0.0) 
		ga = -T_PI/(x*ga*sin(T_PI*x));
	}
  }
  return ga;
}

// Convolves rovibrational DOSs
void Convolution(const std::vector<double> &f1,
  const std::vector<double> &f2,
  std::vector<double> &conv,
  const int n = 0);

// Convolves rovibrational DOSs
void FastLaplaceConvolution(const std::vector<double> &data, const std::vector<double> &respns, std::vector<double> &convolution);

void getCellEnergies(const size_t &cellNumber, const double &cellSize, std::vector<double>& cellEne);

bool convertToFourierCoefficients(const size_t expansion, vector<double>& ak, vector<double>& bk, double& a0, vector<double> pesEnes);

double getEnergyFromFourierCoefficients(double theta, const vector<double> ak, const vector<double> bk, const double a0);

// airy function used for WKB transmission probabilities, which is approximate for deep tunnelling
// but accurate in the shallow tunnelling regime

void airy(double x, double& ai, double& aip, double& bi, double& bip);

// airy2 and its associated functions are accurate over the entire tunnelling regime, 
// including deep tunnelling

void airy2(const double x, double &ai);

void bessik(const double x, const double xnu, double &ri, double &rk, double &rip, double &rkp);

void bessjy(const double x, const double xnu, double &rj, double &ry, double &rjp, double &ryp);

void beschb(const double x, double &gam1, double &gam2, double &gampl, double &gammi);

double chebev(const double a, const double b, std::vector<double> &c, const int m, const double x);

double ModifiedBessalFunction(const double x) ;

double ChiSquaredPrbFn(double x, double a) ;

void nrerror(std::string message);

// end airy2 functions

// Calculate the Cosine Fourier coefficients for general data.
void FourierCosCoeffs(vector<double> &angle, vector<double> &value, vector<double> &cosCoeff, size_t expansion) ;

// Calculate the Sine Fourier coefficients for general data.
void FourierSinCoeffs(vector<double> &angle, vector<double> &value, vector<double> &sinCoeff, size_t expansion) ;

// cumulative distribution function and its inverse
double RationalApproximation(double t) ;

double NormalCDFInverse(double p) ;

// Utility function to fit a polynomial to a suppled data set.
vector<double> FitPoly(size_t order,
  vector<double>::const_iterator xstart,
  vector<double>::const_iterator xend,
  vector<double>::const_iterator ystart);

#endif // GUARD_MesmerMath_h
