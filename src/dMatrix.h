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
#include "TMatrix.h"
#include "MesmerPrecision.h"

namespace mesmer
{
  // double version of Matrix
  class dMatrix : public TMatrix<double>{
    public: dMatrix(int n, const double init = 0.0) : TMatrix<double>(n, init) { } ;
  };

  // double-double version of Matrix
  class ldMatrix : public TMatrix<long double>{
    public: ldMatrix(int n, const long double init = 0.0) : TMatrix<long double>(n, init) { } ;
  };

  // double-double version of Matrix
  class ddMatrix : public TMatrix<dd_real>{
    public: ddMatrix(int n, const dd_real init = 0.0) : TMatrix<dd_real>(n, init) { } ;
  };

  // quad-double version of Matrix
  class qdMatrix : public TMatrix<qd_real>{
    public: qdMatrix(int n, const qd_real init = 0.0) : TMatrix<qd_real>(n, init) { } ;
  };
}//namespacer mesmer


#endif // GUARD_dMatrix_h
