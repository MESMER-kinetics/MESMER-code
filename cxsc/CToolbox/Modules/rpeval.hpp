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
// File: rpeval (header)
// Purpose: Evaluation of a real polynomial p with maximum accuracy.
// Global functions:
//    RPolyEval()      : computes an enclosure for the exact value p(t) with
//                       maximum accuracy
//    RPolyEvalErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __RPEVAL_HPP
#define __RPEVAL_HPP

#include <interval.hpp>     // Interval arithmetic
#include <rpoly.hpp>        // Real polynomial type

using namespace cxsc;
using namespace std;

extern char* RPolyEvalErrMsg ( int );
extern void  RPolyEval ( RPolynomial, real, real&, interval&, int&, int& );
#endif
