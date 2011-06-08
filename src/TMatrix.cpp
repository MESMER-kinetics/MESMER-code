//-------------------------------------------------------------------------------------------
//
// TMatrix.cpp
//
// Author: Struan Robertson
// Date:   7/Apr/2011
//
// This contains the implementation of the the TMatrix class. As TMatrix is a template class
// this file will largely contain particular specializations. 
//
//-------------------------------------------------------------------------------------------
#include "TMatrix.h"

namespace mesmer
{
  //
  // Template specialization for complex numbers.
  //
  template<>
  void TMatrix<complex<double> >::diagonalize(complex<double> *rr)  {
    throw std::runtime_error("The direct diagonalization of complex matrices is not yet implemented.");
  }


}//namespacer mesmer

