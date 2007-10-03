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
// File: hess_ari (header)
// Purpose: Definition of a multi-dimensional interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to second order (i.e. Hessian matrix).
// Types:
//    HTscalar_FctPtr         : pointer for a scalar valued function
//    HTvector_FctPtr         : pointer for a vector valued function
// Class LowTriMatrix:
//    LowTriMatrix()          : constructors
//    ~LowTriMatrix()         : destructor
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for triangular matrices
//                              and real values
//    operator []             : component access
// Class HessType:
//    HessType()              : constructors
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for arguments of types
//                              HessType, interval, and real
//    Dim()                   : returns the dimension of the Hessian matrix
//    HessVar()               : to define HessType variables
//    fValue(), gradValue(),
//    hessValue()             : to get function and derivative values
//    operators +, -, *, /    : operators of diff. arithmetic
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEvalH()                : to compute function value only
//    fgEvalH()               : to compute function and first derivative
//                              value (gradient)
//    fghEvalH()              : to compute function value, gradient, and
//                              Hessian matrix value
// Class HTvector:
//    HTvector()              : constructors
//    ~HTvector()             : destructor
//    operator =              : assignment operator
//    operator []             : component access
//    Dim()                   : to get the actual dimension
//    fValue(), JacValue()    : to get function and derivative values
//    fEvalJ()                : to compute function value only
//    fJEvalJ()               : to compute function and first derivative
//                              value (Jacobian)
//----------------------------------------------------------------------------
#ifndef __HESS_ARI_HPP
#define __HESS_ARI_HPP

#include <interval.hpp>     // Interval arithmetic
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mvi_util.hpp>     // Interval matrix/vector utility functions

using namespace cxsc;
using namespace std;

class LowTriMatrix {
    int      dim;
    ivector *rows;

    void init ( int );
    friend void Resize ( LowTriMatrix&, int );

  public:
    LowTriMatrix ( );
    explicit LowTriMatrix ( int );
    LowTriMatrix ( LowTriMatrix& );
    ~LowTriMatrix ( );

    LowTriMatrix& operator= ( const LowTriMatrix& );
    LowTriMatrix& operator= ( const real& );
    ivector& operator[] ( int );
    const ivector& operator[] ( int ) const;
};

class HessType;
class HTvector;

typedef HessType (*HTscalar_FctPtr)(const HTvector&);
typedef HTvector (*HTvector_FctPtr)(const HTvector&);

class HessType {
    int              nmax;
    interval         f;
    ivector          g;
    LowTriMatrix h;

    friend void Resize ( HessType&, int );
    friend void TestSize ( const HessType&, const HessType&, const char* );

  public:
    HessType ( );
    explicit HessType ( int );
    HessType ( const HessType& );

    HessType& operator= ( const HessType& );
    HessType& operator= ( const interval& );
    HessType& operator= ( const real& );

    int  Dim ( ) const { return nmax; };

    friend HTvector HessVar ( const ivector& );
    friend HTvector HessVar ( const rvector& );

    friend interval fValue    ( const HessType& );
    friend ivector  gradValue ( const HessType& );
    friend imatrix  hessValue ( const HessType& );

    friend HessType operator+ ( HessType&);
    friend HessType operator- ( const HessType&);
    friend HessType operator+ ( const HessType&, const HessType&);
    friend HessType operator- ( const HessType&, const HessType&);
    friend HessType operator* ( const HessType&, const HessType&);
    friend HessType operator/ ( const HessType&, const HessType&);
    friend HessType operator+ ( const HessType&, const interval&);
    friend HessType operator- ( const HessType&, const interval&);
    friend HessType operator* ( const HessType&, const interval&);
    friend HessType operator/ ( const HessType&, const interval&);
    friend HessType operator+ ( const interval&, const HessType&);
    friend HessType operator- ( const interval&, const HessType&);
    friend HessType operator* ( const interval&, const HessType&);
    friend HessType operator/ ( const interval&, const HessType&);
    friend HessType operator+ ( const HessType&, const real&);
    friend HessType operator- ( const HessType&, const real&);
    friend HessType operator* ( const HessType&, const real&);
    friend HessType operator/ ( const HessType&, const real&);
    friend HessType operator+ ( const real&, const HessType&);
    friend HessType operator- ( const real&, const HessType&);
    friend HessType operator* ( const real&, const HessType&);
    friend HessType operator/ ( const real&, const HessType&);

    friend HessType sqr   ( const HessType& );
    friend HessType power ( const HessType&, const int );
    friend HessType sqrt  ( const HessType& );
    friend HessType exp   ( const HessType& );
    friend HessType ln    ( const HessType& );
    friend HessType sin   ( const HessType& );
    friend HessType cos   ( const HessType& );
    friend HessType tan   ( const HessType& );
    friend HessType cot   ( const HessType& );
    friend HessType asin  ( const HessType& );
    friend HessType acos  ( const HessType& );
    friend HessType atan  ( const HessType& );
    friend HessType acot  ( const HessType& );
    friend HessType sinh  ( const HessType& );
    friend HessType cosh  ( const HessType& );
    friend HessType tanh  ( const HessType& );
    friend HessType coth  ( const HessType& );
    friend HessType asinh ( const HessType& );
    friend HessType acosh ( const HessType& );
    friend HessType atanh ( const HessType& );
    friend HessType acoth ( const HessType& );

    friend void fEvalH   ( HTscalar_FctPtr, ivector, interval& );
    friend void fgEvalH  ( HTscalar_FctPtr, ivector, interval&, ivector& );
    friend void fghEvalH ( HTscalar_FctPtr, ivector, interval&, ivector&,
                                                                imatrix& );
};

class HTvector {
    int       nmax;
    HessType *ht;

    friend void TestSize ( const HTvector&, const HTvector&, const char* );

  public:
    explicit HTvector  ( int );
    HTvector  ( const HTvector& );
    ~HTvector ( );

    HTvector& operator=  ( const HTvector& );
    HessType& operator[] ( int );
    const HessType& operator[] ( int ) const;

    int Dim ( ) const { return nmax; };

    friend ivector  fValue   ( const HTvector& );
    friend imatrix  JacValue ( const HTvector& );

    friend void fEvalJ  ( HTvector_FctPtr, ivector, ivector& );
    friend void fJEvalJ ( HTvector_FctPtr, ivector, ivector&, imatrix& );
};
#endif
