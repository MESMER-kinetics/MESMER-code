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
// File: i_util (implementation)
// Purpose: Utilities of type 'interval'
// Global functions:
//    in()      : contained-in relation for a real
//    in()      : contained-in-the-interior relation for two intervals
//    rnd()     : to round a dotprecision argument to an interval
//    Blow()    : epsilon inflation
//    Disjoint(): test for disjointedness of two intervals
//    AbsMin()  : smallest absolute value of an interval
//    AbsMax()  : greatest absolute value of an interval
//    RelDiam() : relative diameter of an interval
//    UlpAcc()  : to check whether the width of an interval is less than a
//                certain number of ulps (ulp = unit in the last place of
//                the mantissa)
//    Power()   : exponentiation by an integer for intervals
//----------------------------------------------------------------------------
#include <imath.hpp>       // Interval mathematical functions
#include <imatrix.hpp>     // Interval matrix/vector arithmetic
#include <i_util.hpp>

// obsolete file

