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
// File: gop (header)
// Purpose: Computing enclosures for all global minimizers and for the global
//    minimum value of a twice continuously differentiable multi-dimensional,
//    scalar valued function, assuming that the global minimum is a stationa-
//    ry point. If it is a boundary point of the search area with gradient of
//    the function being different from zero, the method fails in its form
//    presented here.
// Global functions:
//    AllGOp()      : computes enclosures for all zeros
//    AllGOpErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __GOP_HPP
#define __GOP_HPP

#include <intvector.hpp>    // Integer vector type
#include <interval.hpp>     // Interval arithmetic
#include <hess_ari.hpp>     // Hessian differentiation arithmetic

// Additional header files to get access to interval standard functions and
// interval utility functions which often are needed for implementing user
// defined functions of type 'HessType'
//--------------------------------------------------------------------------
#include <i_util.hpp>       // Interval utility functions
#include <imath.hpp>        // Interval mathematical functions

using namespace cxsc;
using namespace std;

const int MaxCountGOp = 10000; // Maximum count of result components

extern char* AllGOpErrMsg ( int );
extern void  AllGOp ( HTscalar_FctPtr, ivector, real, imatrix&, intvector&,
                      int&, interval&, int&, int = MaxCountGOp );
#endif
