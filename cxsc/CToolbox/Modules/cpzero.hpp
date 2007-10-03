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
// File: cpzero (header)
// Purpose: Determination and enclosure of a root of a complex polynomial,
//    and of the coefficients of the deflated polynomial.
// Global functions:
//    CPolyZero()      : computes an enclosure for a root and for the
//                       deflated complex polynomial
//    CPolyZeroErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __CPZERO_HPP
#define __CPZERO_HPP

#include <cinterval.hpp>  // Complex interval arithmetic
#include <cpoly.hpp>      // Complex polynomial type
#include <cipoly.hpp>     // Complex interval polynomial type

using namespace cxsc;
using namespace std;

extern char* CPolyZeroErrMsg ( int );
extern void  CPolyZero ( CPolynomial, complex,
                         CIPolynomial&, cinterval&, int& );
#endif
