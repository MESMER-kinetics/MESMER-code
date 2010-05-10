#ifndef GUARD_oneDimensionalFourierGridHamiltonian_h
#define GUARD_oneDimensionalFourierGridHamiltonian_h

#include <cmath>
#include "dMatrix.h"

namespace mesmer{

	int oneDimensionalFourierGridHamiltonian(const double imu, // reduced moment of inertia
																					 double (*vsub)(double),
																					 std::vector<double>& eigenvalues,
																					 dMatrix& eigenvectors,
																					 const int numberGridPoint
																					 );
};

#endif