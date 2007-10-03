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

/* CVS $Id: l_rmath.hpp,v 1.18 2007/05/16 14:59:01 cxsc Exp $ */

#ifndef l_rmath_hpp_included
#define l_rmath_hpp_included

#include "l_real.hpp"

namespace cxsc {

//! Calculates \f$ x^2  \f$
inline l_real sqr     (const l_real&) throw(); // Sqr(x)
//! Calculates \f$ \sqrt{x}  \f$
       l_real sqrt    (const l_real&) throw(ERROR_LREAL_STD_FKT_OUT_OF_DEF);
                                               // Sqrt(x)
//! Calculates \f$ \sqrt{x^2+y^2} \f$
       l_real sqrtx2y2(const l_real&, const l_real&) throw(); 
                                               // Sqrt(x^2+y^2)
//! Calculates \f$ \sqrt{1+x^2} \f$
       l_real sqrt1px2(const l_real&) throw(); // Sqrt(1+x^2)
//! Calculates \f$ \sqrt{(x+1)-1} \f$
inline l_real sqrtp1m1(const l_real &) throw(); 
//! Calculates \f$ \sqrt{x^2-1} \f$
inline l_real sqrtx2m1(const l_real &) throw();
//! Calculates \f$ \sqrt{1-x^2} \f$
inline l_real sqrt1mx2(const l_real &) throw();
//! Calculates \f$ \exp(x)-1 \f$
inline l_real expm1   (const l_real &x) throw();
//! Calculates \f$ \exp(-x^2) \f$
inline l_real expmx2  (const l_real&) throw();
//! Calculates \f$ \ln{\sqrt{x^2+y^2}} \f$
inline l_real ln_sqrtx2y2(const l_real& x, const l_real& y) throw();
//! Calculates \f$ \arccos(1+x) \f$
inline l_real acoshp1 (const l_real& x);

// inline l_real sqrt   (const l_real &, int);  // Sqrt(x, n)
// inline l_real sin    (const l_real&) throw();        // Sin(x)
// inline l_real cos    (const l_real&) throw();        // Cos(x)
// inline l_real tan    (const l_real&) throw();        // Tan(x)
// inline l_real cot    (const l_real&) throw();        // Cot(x)
// inline l_real asin   (const l_real&);        // ASin(x)
// inline l_real acos   (const l_real&);        // ACos(x)
// inline l_real atan   (const l_real&);        // ATan(x) 
// inline l_real acot   (const l_real&);        // ACot(x)
// inline l_real exp    (const l_real&) throw();        // Exp(x)
// inline l_real ln     (const l_real&);        // Ln(x)
// inline l_real sinh   (const l_real&) throw();        // Sinh(x)
// inline l_real cosh   (const l_real&) throw();        // Cosh(x)
// inline l_real tanh   (const l_real&) throw();        // Tanh(x) 
// inline l_real coth   (const l_real&) throw();        // Coth(x)
// inline l_real asinh  (const l_real&);        // ASinh(x)
// inline l_real acosh  (const l_real&);        // ACosh(x)
// inline l_real atanh  (const l_real&);        // ATanh(x)
// inline l_real acoth  (const l_real&);        // ACoth(x)

//! Calculates \f$ x^y \f$
inline l_real pow    (const l_real&, const l_real&); // Pow(x,y)
//! Calculates \f$ x^n \f$
l_real power         (const l_real&, const int);     // Power(x,n)

} // namespace cxsc 

#include "l_rmath.inl"
#endif

