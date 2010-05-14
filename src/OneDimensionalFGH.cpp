//------------------------------------------------------------------------------------------
//    Translated to C and import to MESMER by Chi-Hsiu Liang
//
//    FGHEVEN
//    F. Gogtas, G.G. Balint-Kurti and C.C. Marston,
//    QCPE Program No. 647, (1993).
//
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
//       XA  : X-coordinates i.e. position on a grid.
//
//       AR  : The Hamiltonian Matrix
//
//       WCH : The eigenvalues for respective energy levels.
//
//       ZR  : The eigenvectors (X,Y); where
//                                     X : Wavefunction.
//                                     Y : Energy level.
//    The folowing data constans represent :
//       ima, imb             : Moments of inertia of the fragments
//       gridLength           : Grid length (360 degree or of radian 2 pi, default to be degree)
//       gridSpacing          : Grid spacings (degree or radian)
//
//    All quantities in au (unless otherwise stated)
//------------------------------------------------------------------------------------------
#include "OneDimensionalFGH.h"

namespace mesmer{

	int oneDimensionalFourierGridHamiltonian(const double imu, // reduced moment of inertia
																					 double (*vsub)(double, const std::vector<double>, const std::vector<double>, const double),
																					 std::vector<double>& eigenvalues,
																					 dMatrix& eigenvectors,
																					 const int numberGridPoint,
																					 const std::vector<double> ak,
																					 const std::vector<double> bk,
																					 const double a0
																					 ){
		const double psq = M_PI * M_PI;
		const double gridLength = 1.0; // temporary assigned value
		const double const1 = psq / (imu * (gridLength * gridLength));
		const double nfact1 = (numberGridPoint - 1) * (numberGridPoint - 2);
		const double nfact2 = (numberGridPoint - 2) / 2;
		const double const2 = const1 * (double(nfact1) / 6.0 + 1.0 + double(nfact2));
		const double darg = M_PI / double(numberGridPoint);

		dMatrix hamiltonian(numberGridPoint, numberGridPoint);

		double x = 0.0; // temporary
		const double dx(2.0 * M_PI / double(numberGridPoint)); // temporary
		for (int i(0); i < numberGridPoint; ++i){
			for (int j(0); j <= i; ++j){
				double ijD = (i - j);
				if (!ijD){
					hamiltonian[i][j] = const2;
				}
				else{
					const double ratio=1.0 / sin(darg * double(ijD));
					hamiltonian[i][j]=(int(ijD) % 2 ? 1.0 : -1.0) * const1 * ratio * ratio;
				}
			}
			// Find the potential value at x
			// Add the potential value when the kronecker delta function
			// equals to one, i.e. when i and j are equal
      const double temp = (*vsub)(x, ak, bk, a0);
			hamiltonian[i][i] += temp;
			x += dx;
		}

		// Now fill out Hamiltonian matrix.
		for (int i(0); i < numberGridPoint; ++i){
			for (int j(0); j <= i; ++j){
				hamiltonian[j][i] = hamiltonian[i][j];
			}
		}

		// Now call eigenvalue solvers.
		eigenvalues.clear();
		eigenvalues.resize(numberGridPoint, 0.0);
		hamiltonian.diagonalize(&eigenvalues[0]);

		if (eigenvectors[0] != NULL){
			for (int i(0); i < numberGridPoint; ++i){
				for (int j(0); j < numberGridPoint; ++j){
					eigenvectors[i][j] = hamiltonian[i][j];
				}
			}
		}

		return 0;
	}
};
