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
// File: lst_ari (header)
// Purpose: Definition of a list arithmetic used in connection with an
//    interval bisection method in global optimization for storing pairs of
//    an interval vector and a real value.
// Class Pair:
//    Pair()        : constructors
//    ~Pair()       : destructor
//    operator =    : assignment operator
//    _Pair()       : for explicit generating a pair
//    GetInt()      : get the interval vector component
//    GetFyi()      : get the function value
// Global type:
//    PairPtr       : list of pairs
// Global constant:
//    EmptyList     : the empty list (NULL pointer)
// Global functions and operators:
//    operator +    : adding a new element to a list
//    Next(), Head(): access functions for lists
//    Length()      : access function to length of list
//    FreeAll()     : free complete list
//    MultiDelete() : deletes several elements in a list
//    DelHead()     : deletes the first element of a list
//    ListToMatrix(): transfer a list into a dynamic interval matrix
//----------------------------------------------------------------------------
#ifndef __LST_ARI_HPP
#define __LST_ARI_HPP

#include <intvector.hpp>     // Integer vector type
#include <imatrix.hpp>      // Interval matrix/vector arithmetic

using namespace cxsc;
using namespace std;

// #define EmptyList NULL      // The empty list
#define EmptyList 0         // The empty list
                            //---------------

struct                      // Structure used for a list of pairs
  PairElmt;                 //-----------------------------------

typedef                     // Pointer to a list of pairs
  PairElmt* PairPtr;        //---------------------------

class Pair {                // Pair of an interval vector 'intv'
  private:                  // and the real value 'inf(f(intv))'
    ivector intv;           //----------------------------------
    real    fyi;
  public:
    Pair ( ) { };
    Pair ( const Pair& );

    Pair& operator= ( const Pair& );

    friend Pair    _Pair  ( const ivector&, real );
    friend Pair    _Pair  ( const ivector_slice&, real );
    friend Pair    _Pair  ( const imatrix_subv&, real );
    friend ivector GetInt ( const Pair& );
    friend real    GetFyi ( const Pair& );
};

Pair    _Pair  ( const ivector&, real );
Pair    _Pair  ( const ivector_slice&, real );
Pair    _Pair  ( const imatrix_subv&, real );
ivector GetInt ( const Pair& );
real    GetFyi ( const Pair& );

extern void    MultiDelete  ( PairPtr&, const real& );
extern void    ListToMatrix ( PairPtr, imatrix&, intvector& );
extern PairPtr operator+    ( PairPtr, const Pair& );

extern void    FreeAll      ( PairPtr& );
extern int     Length       ( PairPtr );
extern Pair    Head         ( PairPtr );
extern void    DelHead      ( PairPtr& );
#endif
