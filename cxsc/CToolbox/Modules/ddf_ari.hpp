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
// File: ddf_ari (header)
// Purpose: Definition of an interval differentiation arithmetic which allows
//    function evaluation with automatic differentiation up to second order.
// Type:
//    ddf_FctPtr              : pointer for a function of type 'DerivType'
// Class DerivType:
//    DerivType()             : constructors
//    operators +, -, *, /    : operators of diff. arithmetic
//    operator =              : assignment operator
//    DerivConst()
//    DerivVar()              : to define derivative constants/variables
//    fValue()
//    dfValue()
//    ddfValue()              : to get function and derivative values
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEval()                 : to compute function value only
//    dfEval()                : to compute function and first derivative
//                              value
//    ddfEval()               : to compute function, first, and second
//                              derivative value
//----------------------------------------------------------------------------
#ifndef __DDF_ARI_HPP
#define __DDF_ARI_HPP

#include <interval.hpp>     // Interval arithmetic
#include <i_util.hpp>       // Interval utility functions

using namespace cxsc;
using namespace std;

class DerivType;

typedef DerivType (*ddf_FctPtr)(const DerivType&);

class DerivType {
  private:
    interval f, df, ddf;

  public:
    DerivType ( );
    DerivType ( const interval&, const interval&, const interval& );
    DerivType ( const DerivType& );

    DerivType& operator= ( const DerivType& );

    friend DerivType DerivConst ( const real& );
    friend DerivType DerivConst ( const interval& );
    friend DerivType DerivVar   ( const real& );
    friend DerivType DerivVar   ( const interval& );

//    friend interval fValue   ( const DerivType& );
//    friend interval dfValue  ( const DerivType& );
//    friend interval ddfValue ( const DerivType& );

    friend inline const interval fValue   ( const DerivType& u ) { return u.f;   }  
    friend inline const interval dfValue  ( const DerivType& u ) { return u.df;  }  
    friend inline const interval ddfValue ( const DerivType& u ) { return u.ddf; }  

    friend DerivType operator+ ( DerivType& );
    friend DerivType operator- ( const DerivType& );

    friend DerivType operator+ ( const DerivType&, const DerivType& );
    friend DerivType operator- ( const DerivType&, const DerivType& );
    friend DerivType operator* ( const DerivType&, const DerivType& );
    friend DerivType operator/ ( const DerivType&, const DerivType& );

    friend DerivType operator+ ( const DerivType&, const interval& );
    friend DerivType operator- ( const DerivType&, const interval& );
    friend DerivType operator/ ( const DerivType&, const interval& );
    friend DerivType operator* ( const DerivType&, const interval& );

    friend DerivType operator+ ( const interval&, const DerivType& );
    friend DerivType operator- ( const interval&, const DerivType& );
    friend DerivType operator* ( const interval&, const DerivType& );
    friend DerivType operator/ ( const interval&, const DerivType& );

    friend DerivType operator+ ( const DerivType&, const real& );
    friend DerivType operator- ( const DerivType&, const real& );
    friend DerivType operator* ( const DerivType&, const real& );
    friend DerivType operator/ ( const DerivType&, const real& );

    friend DerivType operator+ ( const real&, const DerivType& );
    friend DerivType operator- ( const real&, const DerivType& );
    friend DerivType operator* ( const real&, const DerivType& );
    friend DerivType operator/ ( const real&, const DerivType& );

    friend DerivType sqr   ( const DerivType& );
    friend DerivType power ( const DerivType&, int );
    friend DerivType sqrt  ( const DerivType& );
    friend DerivType exp   ( const DerivType& );
    friend DerivType ln    ( const DerivType& );

    friend DerivType sin    ( const DerivType& );
    friend DerivType cos    ( const DerivType& );
    friend DerivType tan    ( const DerivType& );
    friend DerivType cot    ( const DerivType& );
    friend DerivType asin   ( const DerivType& );
    friend DerivType acos   ( const DerivType& );
    friend DerivType atan   ( const DerivType& );
    friend DerivType acot   ( const DerivType& );

    friend DerivType sinh   ( const DerivType& );
    friend DerivType cosh   ( const DerivType& );
    friend DerivType tanh   ( const DerivType& );
    friend DerivType coth   ( const DerivType& );
    friend DerivType asinh  ( const DerivType& );
    friend DerivType acosh  ( const DerivType& );
    friend DerivType atanh  ( const DerivType& );
    friend DerivType acoth  ( const DerivType& );

    friend void fEval  ( ddf_FctPtr f, interval, interval&);
    friend void dfEval ( ddf_FctPtr f, interval, interval&, interval& );
    friend void ddfEval( ddf_FctPtr f, interval,  interval&, interval&,
                                                             interval& );
};

DerivType DerivConst ( const real& );
DerivType DerivConst ( const interval& );
DerivType DerivVar   ( const real& );
DerivType DerivVar   ( const interval& );
#endif
