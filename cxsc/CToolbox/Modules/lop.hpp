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
// File: lop (header)
// Purpose: Determine enclosures for the optimal value 'z_opt', for the set
//    of optimal basic index sets 'V_opt', and for the set of solution
//    vectors 'X_opt' for a linear programming problem P = (A,b,c) given in
//    the standard form:
//               ( z = c^t * x = max! )
//         (LP)  (     A * x = b      )
//               (      x >= 0        )
//    with an initial optimal basic index set.
// Global functions:
//    LinOpt()      : determines the enclosures of z_opt, V_opt, and X_opt
//                    for LP P = (A,b,c)
//    LinOptErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __LOP_HPP
#define __LOP_HPP

#include <intmatrix.hpp>    // Integer matrix/vector types
#include <rmatrix.hpp>      // Real matrix/vector arithmetic
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mv_util.hpp>      // Real matrix/vector utility functions
#include <mvi_util.hpp>     // Interval matrix/vector utility functions

using namespace cxsc;
using namespace std;

extern char* LinOptErrMsg ( int );
extern void  LinOpt ( rmatrix, rvector, rvector, intvector,
                      interval&, intmatrix&, imatrix&, int&, int& );
#endif
