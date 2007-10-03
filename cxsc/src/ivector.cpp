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

/* CVS $Id: ivector.cpp,v 1.11 2007/05/16 14:58:59 cxsc Exp $ */

#define _CXSC_CPP

#include "ivector.hpp"
#include "vector.inl"
#include "ivector.inl"

namespace cxsc {
int in ( const ivector& x, const ivector& y )          // Inner inclusion for two ivectors
{                                          //---------------------------------
  int i = Lb(x), n = Ub(x), ok = 1;

  while (ok && i <= n) { ok = in(x[i],y[i]); i++; }
  return ok;
}

int in ( int x, ivector& y )     // Inner inclusion of an integer in a ivector
{                                //-------------------------------------------
  int    i = Lb(y), n = Ub(y), ok = 1;
  double d = x;

  while (ok && i <= n) { ok = in(d,y[i]); i++; }
  return ok;
}

ivector Blow ( const ivector& x, real eps )                     // Epsilon deflation
{                                                         //------------------
  int     i;
  ivector h(Lb(x),Ub(x));

  for (i = Lb(x); i <= Ub(x); i++) h[i] = Blow(x[i],eps);
  return h;
}

int Disjoint ( ivector& a, ivector& b )             // Test for disjointedness
{                                                   //------------------------
  int al = Lb(a), au = Ub(a), bl = Lb(b);
  int disjointed, i, d;

  d = bl - al;

  i = al;
  disjointed = 0;
  do {
    if (Disjoint(a[i],b[i+d]))
      disjointed = 1;
    else
      i++;
  } while ( !(disjointed || i > au) );
  return disjointed;
}

int Zero ( ivector& x )                               // Check for zero vector
{                                                     //----------------------
  int i, ok;
  for (i = Lb(x), ok = 1; i <= Ub(x) && ok; i++) ok = (x[i] == 0.0);
  return ok;
}

rvector mid ( ivector& v )                              // Vector of midpoints
{                                                       //--------------------
  int i;
  int l = Lb(v), u = Ub(v);
  rvector w(l,u);

  for (i = l; i <= u; i++) w[i] = mid(v[i]);
  return w;
}

real MaxRelDiam ( const ivector& v )                    // Maximum relative diameter
{                                                 //--------------------------
  real r;
  int  i, l=Lb(v), u=Ub(v);

  r = 0.0;
  for (i = l; i <= u; i++)
    if (RelDiam(v[i]) > r) r = RelDiam(v[i]);
  return r;
}

real MaxRelDiam ( const ivector_slice& v )                    // Maximum relative diameter
{                                                 //--------------------------
  real r;
  int  i, l=Lb(v), u=Ub(v);

  r = 0.0;
  for (i = l; i <= u; i++)
    if (RelDiam(v[i]) > r) r = RelDiam(v[i]);
  return r;
}

//----------------------------------------------------------------------------
// Checks if the diameter of the interval vector 'v' is less or equal to 'n'
// ulps. Ulp is an abbreviation for: units in the last place of the mantissa.
//----------------------------------------------------------------------------
int UlpAcc ( ivector& v, int n )
{
  int k, upper;

  for (k = Lb(v), upper = Ub(v); (k < upper) && UlpAcc(v[k],n); k++);
  return UlpAcc(v[k],n);
}

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------
void DoubleSize ( ivector& x )
{
  int n = Lb(x);
  Resize(x,n,2*Ub(x)-n+1);
}

} // namespace cxsc

