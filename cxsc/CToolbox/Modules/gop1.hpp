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
// File: gop1 (header)
// Purpose: Computing enclosures for all global minimizers and for the global
//    minimum value of a twice continuously differentiable one-dimensional,
//    scalar valued function.
// Global functions:
//    AllGOp1()      : computes enclosures for all global optimizers
//    AllGOp1ErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __GOP1_HPP
#define __GOP1_HPP

#include <interval.hpp>      // Interval arithmetic
#include <intvector.hpp>     // Integer vector type
#include <ivector.hpp>       // Interval vector arithmetic
#include <ddf_ari.hpp>       // Differentiation arithmetic

using namespace cxsc;
using namespace std;

const int MaxCountGOp1 = 10000;  // Maximum count of result components

extern char* AllGOp1ErrMsg ( int );
extern void  AllGOp1 ( ddf_FctPtr, interval, real, ivector&, intvector&,
                       int&, interval&, int&, int = MaxCountGOp1 );
#endif
