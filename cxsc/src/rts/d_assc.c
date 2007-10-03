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

/* CVS $Id: d_assc.c,v 1.9 2007/05/16 14:59:26 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_assc.c                              */
/*                                                              */
/*      Entries         : void d_assc(a,b)                      */
/*                        d_otpc *a, b     ;                    */
/*                                                              */
/*      Arguments       : a = complex dotprecision variable     */
/*                            assigned to                       */
/*                        b = complex dotprecision value to be  */
/*                            assigned                          */
/*                                                              */
/*      Description     : Assignment of intervall dotprecision  */
/*                        variables                             */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local void d_assc(d_otpc *a,d_otpc b)
#else
local void d_assc(a,b)

d_otpc *a;
d_otpc b;
#endif
{
     E_TPUSH("d_assc")
	d_ass(&(a->RE ),b.RE );
     d_ass(&(a->IM ),b.IM );
     E_TPOPP("d_assc")

     return;
}

