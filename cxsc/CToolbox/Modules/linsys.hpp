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
// File: linsys (header)
// Purpose: Computation of a verified solution of a square linear system of
//    equations A*x = b with full real matrix A and real right-hand side b.
// Global functions:
//    LinSolve()      : to get a verified enclosure of the solution (two
//                      versions)
//    LinSolveErrMsg(): to get an error message text
//----------------------------------------------------------------------------
#ifndef __LINSYS_HPP
#define __LINSYS_HPP

#include <rmatrix.hpp>     // Real matrix/vector arithmetic
#include <ivector.hpp>     // Interval vector arithmetic

using namespace cxsc;
using namespace std;

extern char* LinSolveErrMsg ( int );
extern void  LinSolve ( const rmatrix&, const rvector&, ivector&, int& );
extern void  LinSolve ( const rmatrix&, const rvector&, ivector&, real&, int& );
#endif
