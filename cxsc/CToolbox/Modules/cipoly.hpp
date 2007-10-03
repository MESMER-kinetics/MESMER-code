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
// File: cipoly (header)
// Purpose: Definition of the class for the representation of a complex
//    interval polynomial by its coefficients.
// Class CIPolynomial:

//    Deg()         : to get the degree of the polynomial
//    CIPolynomial(): constructors
//    operator []   : component access
//    in()          : Contained-in-the-interior relation for two
//                    complex interval polynomials
//    Blow()        : Epsilon inflation
//    operator >>   : input operator for data of type CIPolynomial
//    operator <<   : output operator for data of type CIPolynomial
//----------------------------------------------------------------------------
#ifndef __CIPOLY_HPP
#define __CIPOLY_HPP

#include <civector.hpp>     // Complex interval vector arithmetic

using namespace cxsc;
using namespace std;

class CIPolynomial {
  private:
    civector coeff;
  public:
    CIPolynomial ( int );
    CIPolynomial ( const CIPolynomial& );
    cinterval& operator[] ( int i ) { return coeff[i]; }
    const cinterval& operator[] ( int i ) const { return coeff[i]; }

    friend int Deg ( const CIPolynomial& );
    friend int in ( const CIPolynomial&, const CIPolynomial& );
    friend CIPolynomial Blow ( CIPolynomial, const real );
    friend istream& operator>> ( istream&, CIPolynomial& );
    friend ostream& operator<< ( ostream&, CIPolynomial );
};
#endif
