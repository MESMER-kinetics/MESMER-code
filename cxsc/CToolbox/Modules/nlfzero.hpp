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
// File: nlfzero (header)
// Purpose: Computing enclosures for all zeros of a continuously
//    differentiable one-dimensional, real-valued function.
// Global functions:
//    AllZeros()      : computes enclosures for all zeros
//    AllZerosErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __NLFZERO_HPP
#define __NLFZERO_HPP

#include <intvector.hpp>     // Integer vector type
#include <ivector.hpp>       // Interval vector arithmetic
#include <mvi_util.hpp>      // Interval matrix/vector utility functions
#include <ddf_ari.hpp>       // Differentiation arithmetic

using namespace cxsc;
using namespace std;

const int MaxCount = 10000;  // Maximum count of result components

extern char* AllZerosErrMsg ( int );
extern void  AllZeros ( ddf_FctPtr, interval, real, ivector&, intvector&,
                        int&, int&, int = MaxCount );
#endif
