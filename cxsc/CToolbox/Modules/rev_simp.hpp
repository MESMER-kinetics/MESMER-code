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
// File: rev_simp (header)
// Purpose: Determine an optimal value 'z', an optimal basic index set 'v',
//    and an optimal solution vector 'x' for a linear programming problem
//    P = (A,b,c) given in the standard form:
//                   ( z = c^t * x = max! )
//         (LP)      (     A * x = b      )
//                   (       x >= 0       ).
// Global function:
//    RevSimplex()      : determines the values z, v, x for the linear
//                        programming problem P = (A,b,c)
//    RevSimplexErrMsg(): to get an error message text
//----------------------------------------------------------------------------
#ifndef __REV_SIMP_HPP
#define __REV_SIMP_HPP

#include <rmatrix.hpp>      // Real matrix/vector arithmetic
#include <intvector.hpp>    // Integer vector type
#include <mv_util.hpp>      // Matrix/vector utility functions
                            // (needed for output of integer vectors)

using namespace cxsc;
using namespace std;

extern char* RevSimplexErrMsg( int );
extern void  RevSimplex      ( rmatrix, rvector, rvector,
                               rvector&, intvector&, real&, int& );
#endif
