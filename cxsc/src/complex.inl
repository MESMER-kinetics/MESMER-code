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

/* CVS $Id: complex.inl,v 1.12 2007/05/16 14:58:56 cxsc Exp $ */


namespace cxsc {
// ---- Constructors ----------------------------------------------

inline complex & complex::operator= (const real & r) throw()
{
   re=r;im=0;
   return *this;
}

      // ---- Std.Operators ---------------------------------------
inline complex operator -(const complex &a) throw () 
{
   return complex(-a.re,-a.im);
}

inline complex operator +(const complex &a) throw ()
{
   return a;
}

inline complex operator +(const complex &a,const complex &b) throw()
{
   return complex(a.re+b.re,a.im+b.im);
}

inline complex operator -(const complex &a,const complex &b) throw()
{
   return complex(a.re-b.re,a.im-b.im);
}

inline complex & operator +=(complex &a, const complex &b) throw() { return a=a+b; }
inline complex & operator -=(complex &a, const complex &b) throw() { return a=a-b; }
inline complex & operator *=(complex &a, const complex &b) throw() { return a=a*b; }
inline complex & operator /=(complex &a, const complex &b) throw() { return a=a/b; }

inline complex operator +(const complex & a,const real & b) throw() 
{ 
   return complex(a.re+b,a.im);
}

inline complex operator +(const real & a,const complex & b) throw()
{
   return complex(a+b.re,b.im);
}

inline complex operator -(const complex & a,const real & b) throw()
{
   return complex(a.re-b,a.im);
}

inline complex operator -(const real & a,const complex & b) throw()
{
   return complex(a-b.re,-b.im);
}

inline complex operator *(const complex & a,const real & b) throw()
{
//   return a*_complex(b);
     return complex(a.re*b,a.im*b);  // Blomquist, 07.11.02;
}

inline complex operator *(const real & a,const complex & b) throw()
{
//   return _complex(a)*b;
     return complex(a*b.re, a*b.im);  // Blomquist, 07.11.02;
}

inline complex operator /(const complex & a,const real & b) throw()
{
//   return a/_complex(b);
     return complex(a.re/b, a.im/b);  // Blomquist, 07.11.02;
}

inline complex operator /(const real & a,const complex & b) throw()
{
   return _complex(a)/b;
}

inline complex & operator +=(complex & a, const real & b) throw() { return a=a+b; }
inline complex & operator -=(complex & a, const real & b) throw() { return a=a-b; }
inline complex & operator *=(complex & a, const real & b) throw() { return a=a*b; }
inline complex & operator /=(complex & a, const real & b) throw() { return a=a/b; }

      // ---- Comp.Operat.  ---------------------------------------
inline bool operator!  (const complex & a)                    throw() { return !a.re && !a.im; }
inline bool operator== (const complex & a, const complex & b) throw() { return a.re==b.re && a.im==b.im; }
inline bool operator!= (const complex & a, const complex & b) throw() { return a.re!=b.re || a.im!=b.im; }
inline bool operator== (const complex & a, const real & b)    throw() { return !a.im && a.re==b; }
inline bool operator== (const real & a, const complex & b)    throw() { return !b.im && a==b.re; }
inline bool operator!= (const complex & a, const real & b)    throw() { return !!a.im || a.re!=b; }
inline bool operator!= (const real & a, const complex & b)    throw() { return !!b.im || a!=b.re; }

      // ---- Others   -------------------------------------------

inline complex conj(const complex & a) throw() { return complex(a.re,-a.im); }


// ----------- Directed Rounding, Blomquist -------------------------------
// ------------------------------------------------------------------------

   // -------------------- addition --------------------------------

inline complex addd(const complex& a, const complex& b) throw()
{ return complex(addd(a.re,b.re), addd(a.im,b.im)); }

inline complex addu(const complex& a, const complex& b) throw()
{ return complex(addu(a.re,b.re), addu(a.im,b.im)); }

inline complex addd(const complex& a, const real& b) throw()
{ return complex(addd(a.re,b), a.im); }

inline complex addu(const complex& a, const real& b) throw()
{ return complex(addu(a.re,b), a.im); }

inline complex addd(const real& a, const complex& b) throw()
{ return complex(addd(a,b.re), b.im); }

inline complex addu(const real& a, const complex& b) throw()
{ return complex(addu(a,b.re), b.im); }
   // ----------------- subtraction: ----------------------------

inline complex subd(const complex& a, const complex& b) throw()
{ return complex(subd(a.re,b.re), subd(a.im,b.im)); }

inline complex subu(const complex& a, const complex& b) throw()
{ return complex(subu(a.re,b.re), subu(a.im,b.im)); }

inline complex subd(const complex& a, const real& b) throw()
{ return complex(subd(a.re,b), a.im); }

inline complex subu(const complex& a, const real& b) throw()
{ return complex(subu(a.re,b), a.im); }

inline complex subd(const real& a, const complex& b) throw()
{ return complex(subd(a,b.re), -b.im); }

inline complex subu(const real& a, const complex& b) throw()
{ return complex(subu(a,b.re), -b.im); }

   // --------------- multiplikation ------------------------

inline complex muld(const complex &a, const real &b) throw()
{ return complex( muld(a.re,b), muld(a.im,b) ); }

inline complex mulu(const complex &a, const real &b) throw()
{ return complex( mulu(a.re,b), mulu(a.im,b) ); }

inline complex muld(const real &a, const complex &b) throw()
{ return complex( muld(a,b.re), muld(a,b.im) ); }

inline complex mulu(const real &a, const complex &b) throw()
{ return complex( mulu(a,b.re), mulu(a,b.im) ); }

   // -------------- division ---------------------------------

inline complex divd(const complex &a, const real &b) throw()
{ return complex( divd(a.re,b), divd(a.im,b) ); }

inline complex divu(const complex &a, const real &b) throw()
{ return complex( divu(a.re,b), divu(a.im,b) ); }

inline complex divd(const real &a, const complex &b) throw()
{ return divd(_complex(a),b); }

inline complex divu(const real &a, const complex &b) throw()
{ return divu(_complex(a),b); }
} // namespace cxsc

