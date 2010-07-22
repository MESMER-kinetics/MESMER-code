//------------------------------------------------------------------------------------------
//    Translated to C and modified to fit the usage of MESMER by Chi-Hsiu Liang
//
//    =========================
//    Original file information
//    FGHEVEN
//    F. Gogtas, G.G. Balint-Kurti and C.C. Marston,
//    QCPE Program No. 647, (1993).
//    =========================
//    
//    The following information may be different from 
//    This program solves one dimensional Schrodinger equation for
//    bound state eigenvalues and eigenfunctions corresponding to a
//    potential V(x).
//
//    The Fourier Grid Hamiltonian method is used in the program.
//    This method is fully described in:
//    C.C. Marston and G.G.Balint-Kurti, J.Chem. Phys.,91,3571(1989).
//    and
//    G.G. Balint-Kurti, R.N. Dixon and C.C. Marston, Internat. Rev.
//    Phys. Chem.,11, 317 (1992).
//
//    The program uses an even number of grid points. The Hamiltonian
//    matrix elements are calculated using analytic formula described
//    in the second of the above references.
//
//    The analytical Hamiltonian expression given in the reference
//    contains a small error.  The formula should read:
//
//               H(i,j) = {(h**2)/(4*m*(L**2)} *
//                         [(N-1)(N-2)/6 + (N/2)] + V(Xi),   if i=j
//
//               H(i,j) = {[(-1)**(i-j)] / m } *
//                         { h/[2*L*sin(M_PI*(i-j)/N)]}**2 ,     if i#j
//
//    The eigenvalues of the Hamiltonian matrix which lie below the
//    asymptotic (large x) value of V(x) correspond to the bound state
//    energies. The corresponding eigenvectors are the representation
//    of the bound state wavefunctions on a regular one dimensional grid.
//
//    The dimensions of the arrays generally depend on the number of grid
//    points used.  This number numberGridPoint is set in a PARAMETER statement.
//    The following arrays represent  :
//
//       xa  : X-coordinates i.e. position on a grid.
//
//       ar  : The Hamiltonian Matrix
//
//       wch : The eigenvalues for respective energy levels.
//
//       zr  : The eigenvectors (X,Y); where
//                                     x : Wavefunction.
//                                     y : Energy level.
//    The folowing data constans represent :
//       imu                  : Moments of inertia of the fragments (amu * Angstrom^2)
//       zl                   : Grid length (radian 2 pi)
//       dx                   : Grid spacings (radian)
//
//    All quantities in au (unless otherwise stated)
//------------------------------------------------------------------------------------------
#include "OneDimensionalFGH.h"
#include "Constants.h"

using namespace Constants;

namespace mesmer{

	int oneDimensionalFourierGridHamiltonian(const double imu, // reduced moment of inertia
																					 double (*vsub)(double, const std::vector<double>, const std::vector<double>, const double),
																					 std::vector<double>& eValues,
																					 dMatrix& eVectors,
																					 const int nx,
																					 const std::vector<double> ak,
																					 const std::vector<double> bk,
																					 const double a0
																					 ){
    const double psq = M_PI * M_PI;
		const double zl = 2.0 * M_PI; // The total length of the system
		const double const1 = psq / (imu * (zl * zl)); // c^2 * (2 * pi)^2 / (4 * mu * zl^2)
		const double const2 = const1 * ((nx + 2) * (nx - 2) / 6.0 + 1.0);
		const double darg = M_PI / double(nx);

		dMatrix ar(nx, nx); // Hamiltonian matrix

		double xa = 0.0; // initial value
		const double dx = 2.0 * M_PI / double(nx); // stepsize
    ctest << "Energies in Hartree:\n{\n";
		for (int i(0); i < nx; ++i){
			for (int j(0); j <= i; ++j){
				double ijD = double(i) - double(j);
				if (!ijD){
					ar[i][j] = const2;
				}
				else{
					const double ratio = 1.0 / sin(darg * double(ijD));
					ar[i][j]=(int(ijD) % 2 ? 1.0 : -1.0) * const1 * ratio * ratio;
				}
			}
			// Find the potential value at x
			// Add the potential value when the kronecker delta function
			// equals to one, i.e. when i and j are equal
      const double vv = (*vsub)(xa, ak, bk, a0);
      ctest << xa << "\t" << vv << endl;
			ar[i][i] += vv;
			xa += dx;
		}
    ctest << "}\n";

		// Now fill out Hamiltonian matrix.
		for (int i(0); i < nx; ++i){
			for (int j(0); j <= i; ++j){
				ar[j][i] = ar[i][j];
			}
		}

		// Now call eigenvalue solvers.
		eValues.clear();
		eValues.resize(nx, 0.0);
		ar.diagonalize(&eValues[0]);

		if (eVectors[0] != NULL){
			for (int i(0); i < nx; ++i){
				for (int j(0); j < nx; ++j){
					eVectors[i][j] = ar[i][j];
				}
			}
		}

    // output the wavefucntions data
    if (1){
      ctest << "Rotational eigenvalues:\n{\n";
      for (unsigned int i(0); i < nx; ++i){
        ctest << eValues[i] * Hartree_in_RC << endl;
      }
      ctest << "}\n\nWavefunctions:\n{\n";
      for (unsigned int i(0); i < nx; ++i){
        for (unsigned int j(0); j < nx; ++j){
          ctest << eVectors[j][i] << "\t";
        }
        ctest << endl;
      }
    }
		return 0;
	}
};
