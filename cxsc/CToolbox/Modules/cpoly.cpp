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
// File: cpoly (implementation)
// Purpose: Definition of the class for the representation of a complex
//    polynomial by its coefficients.
// Class CPolynomial:
//    Deg()        : to get the degree of the polynomial
//    CPolynomial(): constructors
//    operator []  : component access
//    operator >>  : input operator for data of type CPolynomial
//    operator <<  : output operator for data of type CPolynomial
//----------------------------------------------------------------------------
#include <cpoly.hpp>

using namespace cxsc;
using namespace std;

int Deg ( CPolynomial& p ) { return Ub(p.coeff); }

CPolynomial::CPolynomial( int n )
{
  Resize(coeff,0,n);
  coeff = 0.0;
}

CPolynomial::CPolynomial( CPolynomial& p )
{
  Resize(coeff,0,Deg(p));
  coeff = p.coeff;
}

istream& operator>> ( istream& in, CPolynomial& p )
{
  cout << "  x^0 * ";  in >> p[0];
  for (int i = 1; i <= Deg(p); i++)
    { cout << "+ x^" << i << " * ";  in >> p[i]; }
  return in;
}

ostream& operator<< ( ostream& out, CPolynomial p )
{
  int PolyIsZero, n = Deg(p);

  PolyIsZero = 1;
  for (int i = 0; i <= n; i++) {
    if (p[i] == 0.0) continue;
    if (PolyIsZero)
      out << "  ";
    else
      out << "+ ";
    out << p[i] << " * x^" << i << endl;
    PolyIsZero = 0;
  }
  if (PolyIsZero) out << "  0 (= zero polynomial)" << endl;
  return out;
}
