#ifndef GUARD_dMatrix_h
#define GUARD_dMatrix_h

//-------------------------------------------------------------------------------------------
//
// dMatrix.h
//
// Author: Struan Robertson
// Date:   30/Mar/2003
//
// This header file contains the declaration of the dMatrix class.  This
// class inherits from Matrix and wraps calls to EISPACK functions.
//
//-------------------------------------------------------------------------------------------
#include "TMatrix.h"
#include "MesmerPrecision.h"
#include "logExp.h"

using namespace logExpGroup;

namespace mesmer
{
  // double version of Matrix
  class dMatrix : public TMatrix<double>{
    public: dMatrix(size_t n, const double init = 0.0) : TMatrix<double>(n, init) { } ;
  };

  // double version of Matrix
  class lpdMatrix : public TMatrix<logExp>{
    public: lpdMatrix(size_t n, const double init = 0.0) : TMatrix<logExp>(n, init) { } ;
  };

}//namespacer mesmer


#endif // GUARD_dMatrix_h
