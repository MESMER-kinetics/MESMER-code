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
//
//       Modified for C-XSC++ (Version 0.9.1) Stefan Wedner May 2000
//
//----------------------------------------------------------------------------
//
// File: grad_ari (header)
// Purpose: Definition of a multi-dimensional interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to first order (i.e. gradient or Jacobian matrix).
// Types:
//    GTscalar_FctPtr         : pointer for a scalar valued function
//    GTvector_FctPtr         : pointer for a vector valued function
// Class GradType:
//    GradType()              : constructors
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for arguments of types
//                              GradType, interval, and real
//    GradVar()               : to define GradType variables
//    fValue(), gradValue()   : to get function and derivative values
//    operators +, -, *, /    : operators of diff. arithmetic
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEvalG()                : to compute function value only
//    fgEvalG()               : to compute function and first derivative
//                              value (gradient)
// Class GTvector:
//    GTvector()              : constructors
//    ~GTvector()             : destructor
//    Dim()                   : to get the actual dimension
//    operator =              : assignment operator
//    operator []             : component access
//    fValue(), JacValue()    : to get function and derivative values
//    fEvalJ()                : to compute function value only
//    fJEvalJ()               : to compute function and first derivative
//                              value (Jacobian matrix)
//----------------------------------------------------------------------------
#ifndef __GRAD_ARI_HPP
#define __GRAD_ARI_HPP

#include <interval.hpp>     // Interval arithmetic
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mvi_util.hpp>     // Interval matrix/vector utility functions

using namespace cxsc;
using namespace std;

class GradType;
class GTvector;

typedef GradType (*GTscalar_FctPtr)(const GTvector&);
typedef GTvector (*GTvector_FctPtr)(const GTvector&);

class GradType {
    int     nmax;
    ivector g;

    friend void TestSize ( const GradType&, const GradType&, const char* );

  public:
    GradType ( );
    explicit GradType ( int );
    GradType ( const GradType& );

    friend void Resize ( GradType&, int );
    int  Dim ( ) const { return nmax; };

    interval& operator[] ( int );

    GradType& operator= ( const GradType& );
    GradType& operator= ( const interval& );
    GradType& operator= ( const real& );


    friend GTvector GradVar ( const ivector& );
    friend GTvector GradVar ( const rvector& );

    friend interval fValue    ( const GradType& );
    friend ivector  gradValue ( const GradType& );

    friend GradType operator+ ( GradType&);
    friend GradType operator- ( const GradType&);
    friend GradType operator+ ( const GradType&, const GradType&);
    friend GradType operator- ( const GradType&, const GradType&);
    friend GradType operator* ( const GradType&, const GradType&);
    friend GradType operator/ ( const GradType&, const GradType&);
    friend GradType operator+ ( const GradType&, const interval&);
    friend GradType operator- ( const GradType&, const interval&);
    friend GradType operator* ( const GradType&, const interval&);
    friend GradType operator/ ( const GradType&, const interval&);
    friend GradType operator+ ( const interval&, const GradType&);
    friend GradType operator- ( const interval&, const GradType&);
    friend GradType operator* ( const interval&, const GradType&);
    friend GradType operator/ ( const interval&, const GradType&);
    friend GradType operator+ ( const GradType&, const real&);
    friend GradType operator- ( const GradType&, const real&);
    friend GradType operator* ( const GradType&, const real&);
    friend GradType operator/ ( const GradType&, const real&);
    friend GradType operator+ ( const real&, const GradType&);
    friend GradType operator- ( const real&, const GradType&);
    friend GradType operator* ( const real&, const GradType&);
    friend GradType operator/ ( const real&, const GradType&);
    friend GradType sqr   ( const GradType& );
    friend GradType power ( const GradType&, const int );
    friend GradType sqrt  ( const GradType& );
    friend GradType exp   ( const GradType& );
    friend GradType ln    ( const GradType& );
    friend GradType sin   ( const GradType& );
    friend GradType cos   ( const GradType& );
    friend GradType tan   ( const GradType& );
    friend GradType cot   ( const GradType& );
    friend GradType asin  ( const GradType& );
    friend GradType acos  ( const GradType& );
    friend GradType atan  ( const GradType& );
    friend GradType acot  ( const GradType& );
    friend GradType sinh  ( const GradType& );
    friend GradType cosh  ( const GradType& );
    friend GradType tanh  ( const GradType& );
    friend GradType coth  ( const GradType& );
    friend GradType asinh ( const GradType& );
    friend GradType acosh ( const GradType& );
    friend GradType atanh ( const GradType& );
    friend GradType acoth ( const GradType& );

    friend void fEvalG   ( GTscalar_FctPtr, ivector, interval& );
    friend void fgEvalG  ( GTscalar_FctPtr, ivector, interval&, ivector& );
};

class GTvector {
    int      nmax;
    GradType *gt;

    friend void TestSize ( const GTvector&, const GTvector&, const char* );

  public:
    explicit GTvector  ( int );
    GTvector  ( const GTvector& );
    ~GTvector ( );

    int Dim ( ) const { return nmax; };

    GTvector& operator=  ( const GTvector& );

    GradType& operator[] ( int );
    const GradType& operator[] ( int ) const;

    friend ivector  fValue   ( const GTvector& );
    friend imatrix  JacValue ( const GTvector& );

    friend void fEvalJ  ( GTvector_FctPtr, ivector, ivector& );
    friend void fJEvalJ ( GTvector_FctPtr, ivector, ivector&, imatrix& );
};
#endif
