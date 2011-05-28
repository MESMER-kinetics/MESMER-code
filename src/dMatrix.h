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

namespace mesmer
{
  // double version of Matrix
  typedef TMatrix<double> dMatrix;

  // double-double version of Matrix
  typedef TMatrix<long double> ldMatrix;

  // double-double version of Matrix
  typedef TMatrix<dd_real> ddMatrix;

  // quad-double version of Matrix
  typedef TMatrix<qd_real> qdMatrix;

  // complex version of Matrix
  typedef TMatrix<complex<double> > cMatrix;

}//namespacer mesmer


#endif // GUARD_dMatrix_h
