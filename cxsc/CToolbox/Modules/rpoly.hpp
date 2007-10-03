//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// This program/module is free software for non-commercial use. For details
// on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
// This program/module is distributed WITHOUT ANY WARRANTY. For details,
// see the "Disclaimer / Legal Matters" of the book (page iv).
//
//============================================================================
//----------------------------------------------------------------------------
// File: rpoly (header)
// Purpose: Definition of the class for the representation of a real
//    polynomial by its coefficients.
// Class RPolynomial:
//    Deg()        : to get the degree of the polynomial
//    RPolynomial(): constructors
//    operator []  : component access
//    operator >>  : input operator for data of type RPolynomial
//    operator <<  : output operator for data of type RPolynomial
//----------------------------------------------------------------------------
#ifndef __RPOLY_HPP
#define __RPOLY_HPP

#include <rvector.hpp>     // Real vector arithmetic

using namespace cxsc;
using namespace std;

class RPolynomial {
  private:
    rvector coeff;
  public:
    RPolynomial ( int );
    RPolynomial ( RPolynomial& );
    real& operator[] ( int i ) { return coeff[i]; }

    friend int Deg ( RPolynomial& );
    friend istream& operator>> ( istream&, RPolynomial& );
    friend ostream& operator<< ( ostream&, RPolynomial );
};
#endif
