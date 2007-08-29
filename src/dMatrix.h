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

#include "Matrix.h"
#include <string>

class dMatrix : public Matrix<double> {

public:

  // Constructor
  dMatrix(int n) : Matrix<double>(n) { } ;

  // Wrapped call to LAPACK routine to diagonalise matrix.
  void diagonalize(double *rr) {

    int size ;
    size = static_cast<int>(m_msize) ;

    //  Allocate memory for work array
    double *work = new double[size] ;

    tred2(m_matrix, size, rr, work) ;
    tqli(rr, work, size, m_matrix) ;

    delete [] work ;
  }

private:

  //
  // EISPACK methods for diagonalizing matrix.
  //
  void    tred2   (double **a, int n, double *d, double *e) ;
  void    tqli    (double *d, double *e, int n, double **z) ;
  double  pythag  (double a, double b) ;

} ;


#endif // GUARD_dMatrix_h
