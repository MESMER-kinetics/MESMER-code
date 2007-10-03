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
// File: nlinsys (header)
// Purpose: Computing enclosures for all solutions of systems of nonlinear
//    equations given by continuously differentiable functions.
// Global functions:
//    AllNLSS()      : computes enclosures for all solutions
//    AllNLSSErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __NLINSYS_HPP
#define __NLINSYS_HPP

#include <intvector.hpp>    // Integer vector type
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mv_util.hpp>      // Real matrix/vector utility functions
#include <mvi_util.hpp>     // Interval matrix/vector utility functions
#include <grad_ari.hpp>     // Gradient differentiation arithmetic


using namespace cxsc;
using namespace std;


const int MaxCountNLSS = 10000; // Maximum count of result components

extern char* AllNLSSErrMsg ( int );
extern void  AllNLSS ( GTvector_FctPtr, ivector, real, imatrix&, intvector&,
                       int&, int&, int = MaxCountNLSS );
#endif
