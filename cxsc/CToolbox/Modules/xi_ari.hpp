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
// File: xi_ari (header)
// Purpose: Definition of an extended interval arithmetic which allows the
//    division by an interval containing zero.
// Global type:
//    KindType         : component type of extended intervals
// Global function:
//    EmptyIntval()    : empty set as irregular interval
// Class xinterval:
//    operators %, -, &: operators of extended interval arithmetic
//    operator =       : assignment operator
// Updates:
//    04.03.1996 : Modification of extended interval operations and kind type
//----------------------------------------------------------------------------
#ifndef __XI_ARI_HPP
#define __XI_ARI_HPP

#include <interval.hpp>     // Interval arithmetic
#include <ivector.hpp>      // Interval vector arithmetic

using namespace cxsc;
using namespace std;

typedef enum { Finite, PlusInfty, MinusInfty, Double, Empty } KindType;

extern interval EmptyIntval ( );   // Irregular (empty) interval

class xinterval {                  // Extended intervals according
  private:                         // to the above definition
    KindType  kind;                //-----------------------------
    real      inf, sup;

  public:
    xinterval ( );
    xinterval ( const KindType&, const real&, const real& );
    xinterval ( const xinterval& );

    xinterval& operator= ( const xinterval& );

    friend xinterval operator- ( const real&, const xinterval& );
    friend xinterval operator% ( const interval&, const interval& );
    friend ivector   operator& ( const interval&, const xinterval& );
};

xinterval operator- ( const real&, const xinterval& );
xinterval operator% ( const interval&, const interval& );
ivector   operator& ( const interval&, const xinterval& );
#endif
