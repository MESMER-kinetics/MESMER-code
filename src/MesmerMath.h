#ifndef GUARD_MesmerMath_h
#define GUARD_MesmerMath_h


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
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "Matrix.h"

using namespace mesmer;

template <class T>
const T MesmerGamma(const T& x)
{
	T T_PI = acos(-1.);
	int i,k,m;
	T ga,gr,r,z;

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

	if (x > 171.0) return 1e308;    // This value is an overflow flag.
	if (to_double(x) == (int) to_double(x)) {
		if (x > 0.0) {
			ga = 1.0;               // use factorial
			for (i=2;i<x;--i) ga *= i;
		}
		else ga = 1e308;
	}
	else {
		if (abs(x) > 1.0) {
			z = abs(x);
			m = (int)to_double(z);
			r = 1.0;
			for (k=1;k<=m;++k) r *= (z-k);
			z -= m;
		}
		else z = x;
		gr = g[24];
		for (k=23;k>=0;--k) gr = gr*z+g[k];
		ga = 1.0/(gr*z);
		if (abs(x) > 1.0) {
			ga *= r;
			if (x < 0.0) ga = -T_PI/(x*ga*sin(T_PI*x));
		}
	}
	return ga;
}
template<class T>
inline const T SQR(const T a) {return a*a;}


//convolutes rovibrational DOSs
void Convolution(const std::vector<double> &f1,
								 const std::vector<double> &f2,
								 std::vector<double> &conv,
								 const int n = 0);

//convolutes rovibrational DOSs
void FastLaplaceConvolution(const std::vector<double> &data, const std::vector<double> &respns, std::vector<double> &convolution);

void getCellEnergies(int cellNumber, std::vector<double>& cellEne);

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

void nrerror(std::string message);

// end airy2 functions

//---------------------------------------------------------------------------
// MULTIPLY PORTIONS OF TWO MATRICES: matrix1 and matrix2
// If the complete multiplication procedure is (M'= A^T M B), then M is the target matrix. 
// A and B have to be square matrices, but M does not have to be a square matrix.
// This function is either multiplying the (A^T M)or (M B) part.
// The rule here is simple: if matrix1 is     transposed, it is A, and therefore matrix2 is M; 
//                          if matrix1 is NOT transposed, it is M, and therefore matrix2 is B.
// Both matrix1 and matrix2 are constant matrices.
// matrix3 is the location where the result of the multiplication is placed.
template<class T> // The first two numbers are the indices of x and y, the second two numbers are the sizes of the matrix.
bool matrices_multiplication
(const Matrix<T>* matrix1, const size_t m1x, const size_t m1y, const size_t n1x, const size_t n1y,
 const Matrix<T>* matrix2, const size_t m2x, const size_t m2y, const size_t n2x, const size_t n2y,
 Matrix<T>* const matrix3, const bool transposeM1){
	 // matrix3 is a constant pointer but variable T
	 //-----------------------------  1st check if the matrices match up and valid  -----------------------------
	 if (n1x < 1 || n1y < 1 || n2x < 1 || n2y < 1) return false;

	 //-----------------------------  2nd define the dimension of the output matrix  -----------------------------
	 // This procedure requires first to make sure that the non-target matrix is a square matrix, and whether the summation
	 // indices are of the same size.
	 size_t n3x(0), n3y(0);

	 if (n1y != n2x) return false;
	 else if (!transposeM1 && n2x == n2y){ n3x = n1x; n3y = n1y; }
	 else if ( transposeM1 && n1x == n1y){ n3x = n2x; n3y = n2y; }
	 else return false;

	 Matrix<T> mpMul((n3x > n3y) ? n3x : n3y); // whichever is larger

	 //-----------------------------  3rd check if the indices go out of range  -----------------------------
	 if (m1x + n1x - 1 > matrix1->size() || m1y + n1y - 1 > matrix1->size()) return false;
	 if (m2x + n2x - 1 > matrix2->size() || m2y + n2y - 1 > matrix2->size()) return false;

	 //-----------------------------  4th do the multiplication  -----------------------------
	 if (transposeM1){ // When the second matrix is the target, the first matrix has to be transposed.
		 for(size_t i(0); i < n1x; ++i){
			 for(size_t j(0); j < n2y; ++j){
				 T sm = 0.0;
				 const size_t i1x = i + m1x;
				 const size_t j2y = j + m2y;
				 for(size_t k(0); k< n1y; ++k){
					 sm += (*matrix1)[k + m1y][i1x] * (*matrix2)[k + m2x][j2y];
				 }
				 mpMul[i][j] = sm;
			 }
		 }
	 }
	 else{
		 for(size_t i(0); i < n1x; ++i){
			 for(size_t j(0); j < n2y; ++j){
				 T sm = 0.0;
				 const size_t i1x = i + m1x;
				 const size_t j2y = j + m2y;
				 for(size_t k(0); k< n1y; ++k){
					 sm += (*matrix1)[i1x][k + m1y] * (*matrix2)[k + m2x][j2y];
				 }
				 mpMul[i][j] = sm;
			 }
		 }
	 }

	 //-----------------------------  5th copy the temporary matrix to the target position  -----------------------------
	 for(size_t i(0); i < n3x; ++i){
		 for(size_t j(0); j < n3y; ++j){
			 (*matrix3)[i][j] = mpMul[i][j];
		 }
	 }

	 return true;
}

#endif // GUARD_MesmerMath_h
