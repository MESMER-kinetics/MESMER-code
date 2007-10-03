/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.2.1)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2007 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: imatrix.cpp,v 1.11 2007/05/16 14:58:58 cxsc Exp $ */

#define _CXSC_CPP

#include "imatrix.hpp"
#include "matrix.inl"
#include "imatrix.inl"
#include "ivecrmat.inl"

namespace cxsc {
imatrix Id ( imatrix& A )                          // Interval identity matrix
{                                                  //-------------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  imatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = _interval( (i==j) ? 1.0 : 0.0 );
  return B;
}

rmatrix mid ( imatrix& A )                              // Matrix of midpoints
{                                                       //--------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  rmatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = mid(A[i][j]);
  return B;
}

real MaxRelDiam ( const imatrix_subv& v )                    // Maximum relative diameter
{                                                 //--------------------------
  real r;
  int  i, l=Lb(v), u=Ub(v);

  r = 0.0;
  for (i = l; i <= u; i++)
    if (RelDiam(v[i]) > r) r = RelDiam(v[i]);
  return r;
}

int RowLen ( imatrix& A )           // Length of the rows of a interval matrix
  { return Ub(A,2)-Lb(A,2)+1; }     //----------------------------------------

int ColLen ( imatrix& A )        // Length of the columns of a interval matrix
  { return Ub(A,1)-Lb(A,1)+1; }  //-------------------------------------------

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------

void DoubleSize ( imatrix& A )
{
  int n = Lb(A,1);
  Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
}
} // namespace cxsc

