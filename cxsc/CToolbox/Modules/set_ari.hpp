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
// File: set_ari (header)
// Purpose: Implementation of basic operations for sets of integers (also
//    called indices).
// Class IndexSet:
//    IndexSet()   : constructors
//    ~IndexSet()  : destructor
//    operator =   : assignment operator for sets
//    operator +   : unification of two sets or of a set and an index
//    operator -   : difference of two sets or of a set and an index
//    operator []  : returns an element of the set
//    operator ==  : test for equality of two sets
// Global functions:
//    Size()       : returns the actual number of elements in the set
//    Complement() : returns the complementary set
//    SetToVector(): assignment of an index set to an integer vector
//----------------------------------------------------------------------------
#ifndef __SET_ARI_HPP
#define __SET_ARI_HPP

#include <intvector.hpp>     // Integer vector type


using namespace cxsc;
using namespace std;


class IndexSet {
  int   MaxSize;
  char* Index;

  public:
    ~IndexSet ( );
    IndexSet  ( int = 0, char = '\0' );
    IndexSet  ( const IndexSet& );
    IndexSet& operator= ( const IndexSet& );
    IndexSet  operator+ ( const IndexSet& ) const;
    IndexSet  operator- ( const IndexSet& ) const;
    IndexSet  operator+ ( int ) const;
    IndexSet  operator- ( int ) const;
    int       operator[] ( int ) const;
    int       operator== ( const IndexSet& );

  friend int      Size        ( const IndexSet& );
  friend IndexSet Complement  ( const IndexSet& );
  friend void     SetToVector ( const IndexSet&, const intmatrix_subv& );
};
#endif
