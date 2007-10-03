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
// File: cpoly (header)
// Purpose: Definition of the class for the representation of a complex
//    polynomial by its coefficients.
// Class CPolynomial:
//    Deg()        : to get the degree of the polynomial
//    CPolynomial(): constructors
//    operator []  : component access
//    operator >>  : input operator for data of type CPolynomial
//    operator <<  : output operator for data of type CPolynomial
//----------------------------------------------------------------------------
#ifndef __CPOLY_HPP
#define __CPOLY_HPP

#include <cvector.hpp>     // Complex vector arithmetic

using namespace cxsc;
using namespace std;

class CPolynomial {
  private:
    cvector coeff;
  public:
    CPolynomial ( int );
    CPolynomial ( CPolynomial& );
    complex& operator[] ( int i ) { return coeff[i]; }

    friend int Deg ( CPolynomial& );
    friend istream& operator>> ( istream&, CPolynomial& );
    friend ostream& operator<< ( ostream&, CPolynomial );
};
#endif
