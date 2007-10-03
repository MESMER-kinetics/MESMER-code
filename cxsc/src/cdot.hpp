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

/* CVS $Id: cdot.hpp,v 1.15 2007/05/16 14:58:54 cxsc Exp $ */

#ifndef CXSC_CDOT_HPP_INCLUDED
#define CXSC_CDOT_HPP_INCLUDED

#include <iostream>
#include <string>
#include "dot.hpp"
#include "complex.hpp"
#include "l_real.hpp" 
#include "l_complex.hpp"
namespace cxsc {
class l_complex;


//! The Data Type cdotprecision
/*!
The data types dotprecision, idotprecision, cdotprecision and cidotprecision are based on the scalar data types real, interval, complex, and
cinterval, respectively. Variables of these data types permit the exact representation of products of two arbitrary numbers of the
corresponding scalar base type and the exact summation of an arbitrary number of such products in a dotprecision accumulator,
i.e. in a fixed-point format of suitable size.

\sa cxsc::dotprecision
*/
class cdotprecision
{  
   private:
      // ---- Datenelemente ---------------------------------------
      dotprecision re,im;

   public:
      // ---- Constructors  ---------------------------------------
      //! Constructor of class cdotprecision
      cdotprecision() {}
      //! Constructor of class cdotprecision
      inline cdotprecision(const dotprecision&, const dotprecision&);
      //! Constructor of class cdotprecision
      inline cdotprecision(const real &, const real &);
      //! Constructor of class cdotprecision
      inline cdotprecision(const cdotprecision&);
      //! Constructor of class cdotprecision
      inline cdotprecision(const l_real &, const l_real &);
      //! Constructor of class cdotprecision
      inline cdotprecision(const l_complex &);
//    explicit inline cdotprecision(int n)                   : re(n),im(0) { }
//    explicit inline cdotprecision(int a,int b)             : re(a),im(b) { }
//    explicit inline cdotprecision(long n)                  : re(n),im(0) { }
//    explicit inline cdotprecision(long a,long b)           : re(a),im(b) { }
//    explicit inline cdotprecision(const float &d)          : re(d),im(0) { }
//    explicit inline cdotprecision(const float &a, const float &b) : re(a),im(b) { }
//    explicit inline cdotprecision(const double &d)         : re(d),im(0) { }
//    explicit inline cdotprecision(const double &a, const double &b) : re(a),im(b) { }
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const real         &r)   : re(r), im(0) { }
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const complex      &c)   : re(Re(c)),im(Im(c)) { }
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const dotprecision &r)   : re(r), im(0) { }
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const l_real &r)         : re(r), im(0) { }

//    inline cdotprecision& operator= (const double & a)       { re=a;im=0; return *this; }
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const real & a)
                                      { re=a;im=0; return *this; }
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const complex & a) 
                                      { re=Re(a),im=Im(a); return *this; }
//    inline cdotprecision& operator= (const int & a)          { re=a;im=0; return *this; }
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const dotprecision & a)
                                      { re=a;im=0; return *this; }
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const cdotprecision & a)
                                      { re=a.re,im=a.im; return *this; }
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const l_real & a) 
                                      { re=a;im=0; return *this; }

      // ---- Destruktor    ----
      // ~cdotprecision() {} unnoetig

      // ---- Typwandlungen ----
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const dotprecision&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const real&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const l_real&); 
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const complex&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const dotprecision&, const dotprecision&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const real&,const real&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const l_real&, const l_real&); 
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const l_complex&);

      // ---- Ausgabefunkt. ---------------------------------------
      //! Implementation of standard input method
      friend std::istream& operator >> (std::istream& s, cdotprecision& a)       throw();
      //! Implementation of standard output method
      friend std::ostream& operator << (std::ostream& s, const cdotprecision& a) throw();
      //! Implementation of standard input method
      friend std::string&  operator >> (std::string&  s, cdotprecision& a)       throw();
      //! Implementation of standard output method
      friend std::string&  operator << (std::string&  s, const cdotprecision& a) throw();
      //! Implementation of standard input method
      friend void          operator >> (const std::string &s,cdotprecision& a)   throw();
      //! Implementation of standard input method
      friend void          operator >> (const char *s       ,cdotprecision& a)   throw();

      // ---- Standardfunkt ---- (arithmetische Operatoren)
      //! Implementation of standard algebraic negative sign operation
      friend inline cdotprecision operator -(const cdotprecision &) throw();
      //! Implementation of standard algebraic positive sign operation
      friend inline cdotprecision operator +(const cdotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const cdotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const complex &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const complex &,const cdotprecision &) throw();
      
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const dotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const dotprecision &,const cdotprecision &) throw();
      
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const real &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const real &,const cdotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const l_real &, const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const l_real &, const cdotprecision &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const l_real &) throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const dotprecision &) throw();

      // ---- Vergleichsop. ----

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const dotprecision &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const dotprecision &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const dotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const complex &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const complex &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const complex &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const complex &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const real &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const real &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const real &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const real &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard equality operation
      friend inline bool operator ==(const l_real &, const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const l_real &, const cdotprecision &) throw();

      //! Implementation of standard negation operation
      friend inline bool operator !(const cdotprecision&) throw();
      
      // ---- Funktionen    ----
      
      //! Returns the real part of the complex dotprecision value
      friend inline dotprecision & Re(cdotprecision& a);   // { return a.re; }
      //! Returns the imaginary part of the complex dotprecision value
      friend inline dotprecision & Im(cdotprecision& a)  throw(); // { return a.im; }
      //! Returns the real part of the complex dotprecision value
      friend inline const dotprecision & Re(const cdotprecision& a);  // { return a.re; }
      //! Returns the imaginary part of the complex dotprecision value
      friend inline const dotprecision & Im(const cdotprecision& a)  throw(); // { return a.im; }
      //! Sets the real part of a complex dotprecision value
      friend inline cdotprecision& SetRe (cdotprecision& a, const dotprecision& b)  throw();
      //! Sets the imaginary part of a complex dotprecision value
      friend inline cdotprecision& SetIm (cdotprecision& a, const dotprecision& b) throw();
      
      //! Returns the conjugated complex dotprecision value
      friend inline cdotprecision conj(const cdotprecision &a) throw();

      //! Converting the exact complex dotprecision value with one rounding into a complex value
      friend void rnd(const cdotprecision &,complex &,rndtype = RND_NEXT) throw();
      //! Converting the exact complex dotprecision value with one rounding into the nearest lower and upper complex value
      friend void rnd(const cdotprecision &,complex &,complex &) throw();
      //! Converting the exact complex dotprecision value with one rounding into a complex value
      friend complex rnd(const cdotprecision &,rndtype = RND_NEXT) throw();
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend        void accumulate  (cdotprecision&, const complex&, const complex&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cdotprecision&, const complex&, const real&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cdotprecision&, const real&, const complex&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cdotprecision&, const real&, const real&) throw();
};

inline cdotprecision _cdotprecision(const l_complex&);

inline dotprecision & Re(cdotprecision& a)   { return a.re; }
inline dotprecision & Im(cdotprecision& a)  throw() { return a.im; }
inline const dotprecision & Re(const cdotprecision& a) { return a.re; }
inline const dotprecision & Im(const cdotprecision& a)  throw() { return a.im; }

// ---------------- cdotprecision +(-)= l_complex ---------------------------- 
inline cdotprecision & operator += (cdotprecision &cd, 
                                           const l_complex &lc) throw() 
	{
	    cd = cd + _cdotprecision(lc); return cd;
	}
inline cdotprecision & operator -= (cdotprecision &cd, 
                                           const l_complex &lc) throw() 
	{
	    cd = cd - _cdotprecision(lc); return cd;
	}
	

//----------------------------------------------------------------------
// global verfuegbare CDotprecision Akku's   

#define MAXCDOTAKKU     (MAXDOTAKKU / 2)
extern cdotprecision cdotakku[MAXCDOTAKKU];
  
//----------------------------------------------------------------------

} // namespace cxsc 

#include "cdot.inl"

#endif 

