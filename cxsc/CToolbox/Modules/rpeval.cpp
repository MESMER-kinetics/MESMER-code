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
// File: rpeval (implementation)
// Purpose: Evaluation of a real polynomial p with maximum accuracy.
// Method: Transformation of Horner's scheme for evaluating p(t) to a linear
//    system of equations A(t)*x = p, where A(t) is a bidiagonal, Toeplitz
//    matrix and x = (x_n,...,x_0). By iterative refinement the floating-
//    point approximation is improved and the exact solution x is enclosed in
//    an interval vector [x].
// Global functions:
//    RPolyEval()      : computes an enclosure for the exact value p(t) with
//                       maximum accuracy
//    RPolyEvalErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#include <string.h>         // String handling
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mvi_util.hpp>     // Interval matrix/vector utility functions
#include <rpeval.hpp>

using namespace cxsc;
using namespace std;

const int kmax = 10;        // Maximum number of iteration steps

//---------------------------------------------------------------------------
// Error codes used in this module. In the comments below p[0],..., p[n] are
// the coefficients of a polynomial p.
//---------------------------------------------------------------------------
const int
  NoError  = 0,    // No error occurred
  ItFailed = 1;    // Maximum number of iterations exceeded

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* RPolyEvalErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    if (Err == ItFailed)
      sprintf(Hlp,"Maximum number of iterations (=%d) exceeded",kmax);
    else
      strcpy(Hlp,"Code not defined");
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
}

//----------------------------------------------------------------------------
// Purpose: Determination of p(t) (a polynomial p with argument t) with
//    maximum accuracy.
// Parameters:
//    In : 'p'  : represents a real polynomial by its coefficients
//         't'  : specifies the point of evaluation
//    Out: 'z'  : floating-point approximation computed by Horner's scheme
//         'zz' : enclosure of p(t)
//         'k'  : number of iterations needed
//         'Err': error flag
// Description: The polynomial 'p' is defined by its coefficients and, 't'
//    denotes the evaluation point. Horner's scheme for evaluating p is equi-
//    valent to computing the solution of a linear system of equations
//    A(t)*x = p, where A(t) is a bidiagonal, Toeplitz matrix and
//    x =(x_n,...,x_0). The component x_0 of x is equal to p(t). The solution
//    x is enclosed in an interval vector [x] by iterative refinement. The
//    first element of [x] is an enclosure of p(t). It is returned in the
//    variable 'zz'. A floating-point approximation computed by Horner's
//    scheme is returned in 'z'. The number of iterations is returned in 'k'
//    and the state of success is returned in 'Err'.
// Remark: The polynomial's coefficients and the argument are assumed to be
//    exact floating-point numbers!
//----------------------------------------------------------------------------
void RPolyEval ( RPolynomial p,  real t, real& z,
                 interval&   zz, int& k, int&  Err )
{
  int           i, j, n = Deg(p);            // The j-th column of the matrix
  ivector       rr(0,n), yy(0,n);            // 'y' is used to store the j-th
  rvector       x(0,n);                      // correction for the solution of
  rmatrix       y(0,n,0,kmax);               // A(t)*y = p.
  dotprecision  Accu;                        //-------------------------------
  idotprecision IAccu;

  k = 0;                                     // Initialization
  Err = NoError;                             //---------------

  if (n == 0)                                // Constant polynomial
    { z = p[0]; zz = p[0]; }                 //--------------------

  else if (n == 1) {                         // Linear polynomial
    Accu = p[0];                             //------------------
    accumulate(Accu,t,p[1]);
    rnd(Accu,z);
    rnd(Accu,zz);
  }
  else { // n > 1
    x = 0.0; y = 0.0;                       // Initialization x = 0 and y = 0
                                            //--------------------------------

    x[n] = p[n];                            // Computation of a first approxi-
    for (i = n-1; i >= 0; i--)              // mation using Horner's scheme
      x[i] = x[i+1]*t + p[i];               //--------------------------------
    z = x[0];

    y[Col(0)] = x;                          // Iterative refinement for the
    do {                                    // solution of A*x = p
                                            //--------------------------------
      if (k > 0)                    // If a residual was computed, its middle
        y[Col(k)] = mid(yy);        // is stored as the next correction of 'y'

      yy[n] = 0.0;                          // Computation of the residual [r]
      for (i = n-1; i >= 0; i--) {          // and evaluation of the interval
        Accu = p[i];                        // system A*[y] = [r]
        for(j = 0; j <= k; j++) Accu -= y[i][j];
        for(j = 0; j <= k; j++) accumulate(Accu,t,y[i+1][j]);
        rnd(Accu,rr[i]);
        yy[i] = yy[i+1]*t + rr[i];
      }

      IAccu = yy[0];                                 // Determination of a new
      for (j = 0; j <= k; j++) IAccu += y[0][j];     // enclosure [z] of p(t)
      rnd(IAccu,zz);

      k++;
    } while ( !(UlpAcc(zz,1) || (k > kmax)) );

    if ( !UlpAcc(zz,1) ) Err = ItFailed;
  }
}
