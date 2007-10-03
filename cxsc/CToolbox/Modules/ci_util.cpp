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
// File: ci_util (implementation)
// Purpose: Utilities of type 'cinterval'
// Global functions:
//    in()  : Contained-in-the-interior relation for two complex intervals
//    Blow(): Epsilon inflation
//----------------------------------------------------------------------------
#include <i_util.hpp>      // Interval utility functions
#include <ci_util.hpp>

// obsolete file

/*
using namespace cxsc;
using namespace std;

int in ( cinterval& x, cinterval& y )    // Contained-in-the-interior relation
{                                        //-----------------------------------
  return ( in(Re(x),Re(y)) && in(Im(x),Im(y)) );
}

cinterval Blow ( cinterval x, const real& eps )           // Epsilon inflation
{                                                         //------------------
  return _cinterval(Blow(Re(x),eps),Blow(Im(x),eps));
}
*/
