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
// File: expreval (header)
// Purpose: Computation of an enclosure of the value of a real arithmetic
//    expression composed of the operations +, -, *, /, and ^, where ^
//    denotes exponentiation by an integer.
// Type:
//    Stagg_FctPtr        : pointer for a function of type 'Staggered'
// Class Staggered (staggered data representation):
//    Staggered()         : constructors
//    operators =         : assignment of a real argument
//    operators +, -, *, /: both operands of type 'Staggered' or one of
//                          type 'Staggered' and one of type 'real'
//    Power()             : argument of type 'Staggered', exponent of type
//                          'integer' (Exponentiation by an integer)
//    Eval()              : main function for evaluation of an expression
//    EvalErrMsg()        : to get an error message text
// Class StaggArray (arrays of type 'Staggered'):
//    StaggArray()        : constructors
//    ~StaggArray()       : destructor
//    operator []         : component access
//----------------------------------------------------------------------------
#ifndef __EXPREVAL_HPP
#define __EXPREVAL_HPP

#include <interval.hpp>     // Interval arithmetic
#include <i_util.hpp>       // Interval utility functions
#include <rvector.hpp>      // Real vector arithmetic

using namespace cxsc;
using namespace std;

extern char* EvalErrMsg ( int );

class Staggered;
class StaggArray;

typedef Staggered  (*Stagg_FctPtr)(StaggArray&);

class Staggered {
  private:
    rvector   Val;
    interval  Err;

    friend void InitEntry       ( real );
    friend void UpdateError     ( interval );
    friend void UpdateStaggComp ( int );

  public:
    Staggered ( );
    Staggered ( const Staggered& );

    Staggered& operator= ( const real& );
    Staggered& operator= ( const Staggered& );

    friend Staggered  operator+ ( const Staggered& , const Staggered& );
    friend Staggered  operator- ( const Staggered& , const Staggered& );
    friend Staggered  operator* ( const Staggered& , const Staggered& );
    friend Staggered  operator/ ( const Staggered& , const Staggered& );

    friend Staggered  operator+ ( const Staggered& , const real& );
    friend Staggered  operator- ( const Staggered& , const real& );
    friend Staggered  operator* ( const Staggered& , const real& );
    friend Staggered  operator/ ( const Staggered& , const real& );

    friend Staggered  operator+ ( const real& , const Staggered& );
    friend Staggered  operator- ( const real& , const Staggered& );
    friend Staggered  operator* ( const real& , const Staggered& );
    friend Staggered  operator/ ( const real& , const Staggered& );

    friend Staggered Power ( const Staggered& , const int );

    friend void Eval ( Stagg_FctPtr, rvector, real, real&, interval&, int&,
                                                                      int& );

}; // class Staggered

class StaggArray {
  private:
    Staggered *SA;
    int       Dim;

  public:
    StaggArray  ( );
    StaggArray  ( int );
    ~StaggArray ( );
    StaggArray  ( StaggArray& );

    Staggered& operator[] ( int );
}; // class StaggArray
#endif



