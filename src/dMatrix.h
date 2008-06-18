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
		public: dMatrix(int n) : TMatrix<double>(n) { } ;
  };

  // double-double version of Matrix
	class ldMatrix : public TMatrix<long double>{
		public: ldMatrix(int n) : TMatrix<long double>(n) { } ;
	};

  // double-double version of Matrix
	class ddMatrix : public TMatrix<dd_real>{
		public: ddMatrix(int n) : TMatrix<dd_real>(n) { } ;
	};

  // quad-double version of Matrix
	class qdMatrix : public TMatrix<qd_real>{
		public: qdMatrix(int n) : TMatrix<qd_real>(n) { } ;
	};
}//namespacer mesmer


#endif // GUARD_dMatrix_h
