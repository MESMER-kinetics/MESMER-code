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
// File: cipoly (implementation)
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
#include <ci_util.hpp>     // Complex interval utility functions
#include <cipoly.hpp>

using namespace cxsc;
using namespace std;

int Deg ( const CIPolynomial& p ) { return Ub(p.coeff); }

CIPolynomial::CIPolynomial( int n )
{
  Resize(coeff,0,n);
  coeff = 0.0;
}

CIPolynomial::CIPolynomial( const CIPolynomial& p )
{
  Resize(coeff,0,Deg(p));
  coeff = p.coeff;
}

int in ( const CIPolynomial& p, const CIPolynomial& q )               // Contained-in-the-
{                                                         // interior relation
  int i, incl = 1;                                        //------------------

  for (i = 0; (incl == 1) && (i <= Deg(p)); i++)
    incl = in(p[i],q[i]);
  return incl;
}

CIPolynomial Blow ( CIPolynomial p, const real eps )      // Epsilon inflation
{                                                         //------------------
  for (int i = 0; i <= Deg(p); i++)
    p[i] = Blow(p[i],eps);
  return p;
}

istream& operator>> ( istream& in, CIPolynomial& p )
{
  cout << "  x^0 * ";  in >> p[0];
  for (int i = 1; i <= Deg(p); i++)
    { cout << "+ x^" << i << " * ";  in >> p[i]; }
  return in;
}

ostream& operator<< ( ostream& out, CIPolynomial p )
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
