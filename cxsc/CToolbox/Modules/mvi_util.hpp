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
// File: mvi_util (header)
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
#ifndef __MVI_UTIL_HPP
#define __MVI_UTIL_HPP

#include <i_util.hpp>      // Interval utility functions
#include <rmatrix.hpp>     // Real matrix/vector arithmetic
#include <imatrix.hpp>     // Interval matrix/vector arithmetic
#include <mv_util.hpp>     // Utility functions for real matrices, real
 
//using namespace cxsc;
//using namespace std;
                          // vectors, and integer vectors

//extern int      in         ( const ivector&, const ivector& );
//extern int      in         ( int,      ivector& );
//extern ivector  Blow       ( const ivector&, real );
//extern int      Disjoint   ( ivector&, ivector& );
//extern imatrix  Id         ( imatrix& );
//extern int      Zero       ( ivector& );
//extern rmatrix  mid        ( imatrix& );
//extern rvector  mid        ( ivector& );
//extern real     MaxRelDiam ( const ivector& );
//extern real     MaxRelDiam ( const ivector_slice& );
//extern real     MaxRelDiam ( const imatrix_subv& );
//extern int      UlpAcc     ( ivector&, int );
//extern int      VecLen     ( ivector& );
//extern int      RowLen     ( imatrix& );
//extern int      ColLen     ( imatrix& );
//extern void     DoubleSize ( ivector& );
//extern void     DoubleSize ( imatrix& );
#endif
