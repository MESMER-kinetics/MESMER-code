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
// File: lop_ari (header)
// Purpose: Definition of a linearly linked list as an abstract data type for
//    the representation of a list of integer sets.
// Global type:
//    BaseList : representation of a linearly linked list of index sets
// Global constant:
//    EmptyList: the empty list (NULL pointer)
// Global functions:
//    FreeAll(): free complete list
//    extract(): extracts submatrix/subvector depending on the
//               actual index set
//    in()     : returns TRUE, if index set is in list
//    select() : selects first element of a list
//    insert() : inserts index set at the head of a list
//    del()    : deletes index set from list
//    append() : appends 2nd list to end of 1st list
//    remove() : removes elements of 2nd list from 1st list
//----------------------------------------------------------------------------
#ifndef __LOP_ARI_HPP
#define __LOP_ARI_HPP

#include <set_ari.hpp>     // Index set handling
#include <rmatrix.hpp>     // Real matrix/vector arithmetic

using namespace cxsc;
using namespace std;

#define EmptyList NULL           // The empty list
                                 //---------------

struct                           // Structure used for a list of sets
  BaseListElement;               //----------------------------------

typedef                          // Pointer to a list of integer sets
  BaseListElement* BaseList;     //----------------------------------

extern void     FreeAll ( BaseList& );
extern rmatrix  extract ( rmatrix&, const IndexSet& );
extern rvector  extract ( rvector&, const IndexSet& );
extern int      in      ( const IndexSet&, BaseList );
extern IndexSet select  ( BaseList );
extern void     insert  ( BaseList&, const IndexSet& );
extern void     del     ( BaseList&, const IndexSet& );
extern void     append  ( BaseList&, BaseList& );
extern void     remove  ( BaseList&, BaseList );
#endif
