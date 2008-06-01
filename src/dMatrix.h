#ifndef GUARD_dMatrix_h
#define GUARD_dMatrix_h

//-------------------------------------------------------------------------------------------
//
// dMatrix.h
//
// Author: Struan Robertson
// Date:   30/Mar/2003
//
// This header file contains the declaration of the dMatrix class.  This class inherits from
// Matrix and wraps calls to LAPACK functions.
//
//-------------------------------------------------------------------------------------------
#include "MesmerPrecision.h"
#include "Matrix.h"
#include <string>

namespace mesmer
{
	class dMatrix : public Matrix<double> {

	public:

		// Constructor
		dMatrix(int n) : Matrix<double>(n, 0.0) { } ;

		//
		// Wrapped call to EISPACK routine to diagonalise matrix.
		//
		void diagonalize(double *rr) {

			int size ;
			size = static_cast<int>(m_msize) ;

			//  Allocate memory for work array
			double *work = new double[size] ;
			double *rrProxy = new double[size] ;

			tred2(m_matrix, size, rrProxy, work) ;
			tqli(rrProxy, work, size, m_matrix) ;

			for (int i = 0; i < size; ++i){
				rr[i] = to_double(rrProxy[i]);
			}

			delete [] work ;
			delete [] rrProxy ;

		}

		// 
		// Solve a set of linear equations with a single right hand side.
		//
		void solveLinearEquationSet(double *rr) {

			int size ;
			size = static_cast<int>(m_msize) ;

			//  Allocate memory for work array
			int *indx = new int[size] ;

			double d ;
			ludcmp(m_matrix, size, indx, d) ;
			lubksb(m_matrix, size, indx, rr) ;

			delete [] indx ;

		};

	private:

		//
		// EISPACK methods for diagonalizing matrix.
		//
		void    tred2   (double **a, int n, double *d, double *e) ;
		void    tqli    (double *d, double *e, int n, double **z) ;
		double  pythag  (double a, double b) ;

		//
		// NR LU methods for linear equation solving.
		//
		void ludcmp(double **a,  int n, int *indx, double d) ;
		void lubksb(double **a,  int n, int *indx, double* b) ;

	} ;

}//namespacer mesmer


#endif // GUARD_dMatrix_h
