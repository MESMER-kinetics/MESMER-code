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

/* CVS $Id: cimatrix.cpp,v 1.11 2007/05/16 14:58:55 cxsc Exp $ */

#define _CXSC_CPP

#include "cimatrix.hpp"
#include "matrix.inl"
#include "cimatrix.inl"
#include "civecrmat.inl"
#include "civeccmat.inl"
#include "civecimat.inl"

// they should also be included here
#include "cvecimat.inl"
#include "iveccmat.inl"
#include "cmatimat.inl"

namespace cxsc {
int RowLen ( cimatrix& A )       // Length of the rows of a cinterval matrix
  { return Ub(A,2)-Lb(A,2)+1; }  //------------------------------------------

int ColLen ( cimatrix& A )       // Length of the columns of a cinterval matrix
  { return Ub(A,1)-Lb(A,1)+1; }  //---------------------------------------------

} // namespace cxsc

