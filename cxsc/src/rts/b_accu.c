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

/* CVS $Id: b_accu.c,v 1.9 2007/05/16 14:59:19 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_accu.c                              */
/*                                                              */
/*      Description     : Definition of global accumulators.    */
/*                                                              */
/*      Note            : Accumulators are used by old-fashioned*/
/*                        dot product routines only.            */
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

static a_btyp o_accu[4*A_LENGTH];

local dotprecision b_acrl = &o_accu[0];
local dotprecision b_acru = &o_accu[A_LENGTH];
local dotprecision b_acil = &o_accu[2*A_LENGTH];
local dotprecision b_aciu = &o_accu[3*A_LENGTH];

#if VAX_VMS_C
#ifdef LINT_ARGS
local void b_accu(void)
#else
local void b_accu()
#endif
        {
        }
#endif

