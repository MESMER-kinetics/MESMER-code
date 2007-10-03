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
// Example: Evaluation of a real polynomial with maximum accuracy
//----------------------------------------------------------------------------
#include <rpeval.hpp>     // Polynomial evaluation


using namespace cxsc;
using namespace std;


int main ( )
{
  int      Err, No, n;
  real     t, y;
  interval yy;

  do {
    cout << "Enter the degree of the polynomial (>=0): ";
    cin  >> n; cout << endl;
  } while (n < 0);

  RPolynomial p(n);

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Enter the coefficients of p in increasing order:" << endl;
  cin  >> p; cout << endl;
  cout << "Enter the argument t for evaluation:" << endl;
  cin  >> t; cout << endl;

  RPolyEval(p,t,y,yy,No,Err);

  if (!Err) {
    cout << "Polynomial:" << endl << p <<endl;
    cout << "Floating-point evaluation of p(t) using Horner's scheme:"
         << endl << " " << y << endl << endl;
    cout << "Verified inclusion of p(t):"
         << endl << yy << endl << endl;
    cout << "Number of iterations needed: " << No << endl;
  }
  else
    cout << RPolyEvalErrMsg(Err) << endl;

  return 0;
}
