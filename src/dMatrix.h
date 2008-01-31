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
class dMatrix : public Matrix<MesmerHP> {

public:

  // Constructor
  dMatrix(int n) : Matrix<MesmerHP>(n, 0.0) { } ;

  // Wrapped call to LAPACK routine to diagonalise matrix.
  void diagonalize(double *rr) {

    int size ;
    size = static_cast<int>(m_msize) ;

    //  Allocate memory for work array
    MesmerHP *work = new MesmerHP[size] ;
    MesmerHP *rrProxy = new MesmerHP[size] ;
    
    tred2(m_matrix, size, rrProxy, work) ;
    tqli(rrProxy, work, size, m_matrix) ;

    for (int i = 0; i < size; ++i){
      rr[i] = to_double(rrProxy[i]);
    }
    delete [] work ;
  }

private:

  //
  // EISPACK methods for diagonalizing matrix.
  //
  void    tred2   (MesmerHP **a, int n, MesmerHP *d, MesmerHP *e) ;
  void    tqli    (MesmerHP *d, MesmerHP *e, int n, MesmerHP **z) ;
  MesmerHP  pythag  (MesmerHP a, MesmerHP b) ;

} ;

}//namespacer mesmer


#endif // GUARD_dMatrix_h
