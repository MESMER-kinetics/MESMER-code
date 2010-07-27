#ifndef GUARD_oneDimensionalFourierGridHamiltonian_h
#define GUARD_oneDimensionalFourierGridHamiltonian_h

#include <cmath>
#include "dMatrix.h"

namespace mesmer{

	int oneDimensionalFourierGridHamiltonian(const double imu, // reduced moment of inertia
																					 double (*vsub)(double, const std::vector<double>, const std::vector<double>, const double),
																					 std::vector<double>& eValues,
																					 dMatrix& eVectors,
																					 const int nx,
																					 const std::vector<double> ak,
																					 const std::vector<double> bk,
																					 const double a0
																					 );
};

#endif