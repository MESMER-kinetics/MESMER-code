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
// File: mv_util (implementation)
// Purpose: Utilities of type 'intvector', 'intmatrix', 'rvector', and
//    'rmatrix'.
// Global functions:
//    VecLen()    : length of an integer or real vector
//    RowLen()    : length of the rows of an integer or real matrix
//    ColLen()    : length of the columns of an integer or real matrix
//    Id()        : identity real matrix
//    transp()    : transposed of a real matrix
//    DoubleSize(): for doubling the size of an integer vector or matrix
//    operator << : output of integer vector
//----------------------------------------------------------------------------
#include <mv_util.hpp>
/*
using namespace cxsc;
using namespace std;

int VecLen ( const intvector& v )          // Length of an integer vector
  { return Ub(v)-Lb(v)+1; }                //----------------------------

int RowLen ( const intmatrix& A )    // Length of the rows of a integer matrix
  { return Ub(A,2)-Lb(A,2)+1; }      //---------------------------------------

int ColLen ( const intmatrix& A )// Length of the columns of an integer matrix
  { return Ub(A,1)-Lb(A,1)+1; }  //-------------------------------------------

int VecLen ( const rvector& v )                     // Length of a real vector
  { return Ub(v)-Lb(v)+1; }                         //------------------------

int RowLen ( const rmatrix& A )         // Length of the rows of a real matrix
  { return Ub(A,2)-Lb(A,2)+1; }         //------------------------------------

int ColLen ( const rmatrix& A )      // Length of the columns of a real matrix
  { return Ub(A,1)-Lb(A,1)+1; }      //---------------------------------------

rmatrix Id ( const rmatrix& A )                        // Real identity matrix
{                                                      //---------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  rmatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = (i==j) ? 1.0 : 0.0;
  return B;
}

rmatrix transp ( rmatrix& A )                       // Transposed matrix
{                                                         //------------------
  int      n;
  rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
  return res;
}

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------
void DoubleSize ( intvector& x )
{
  int n = Lb(x);
  Resize(x,n,2*Ub(x)-n+1);
}

void DoubleSize ( intmatrix& A )
{
  int n = Lb(A,1);
  Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
}

ostream& operator<< ( ostream& os, intvector& v )         // Output of integer
{                                                         // vectors
  int i, newline = (Ub(v)-Lb(v) > 15);                    //------------------

  for (i = Lb(v); i <= Ub(v); i++) {
    os << v[i] << ' ';
    if (newline) os << endl;
  }
  if (!newline) os << endl;
  return os;
}
*/
