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
// File: mvi_util (implementation)
// Purpose: Utilities of type 'ivector' and 'imatrix'.
// Global functions:
//    in()        : contained-in-the-interior relation for two interval
//                  vectors and for an integer in an interval vector
//    Blow()      : epsilon inflation for an interval vector
//    Disjoint()  : test for disjointedness of two interval vectors
//    Id()        : identity interval matrix
//    Zero()      : check for a zero interval matrix
//    mid()       : midpoint vector and matrix
//    MaxRelDiam(): to get the maximum of the relative diameters of the
//                  components of an interval vector.
//    UlpAcc()    : to check whether all components of an interval vector
//                  have width less than a certain number of ulps
//                  (ulp = unit in the last place of the mantissa).
//    VecLen()    : length of an interval vector
//    RowLen()    : length of the rows of an interval matrix
//    ColLen()    : length of the columns of an interval matrix
//    DoubleSize(): for doubling the size of an interval vector or matrix
//----------------------------------------------------------------------------
#include <mvi_util.hpp>
