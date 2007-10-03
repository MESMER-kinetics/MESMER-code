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

/* CVS $Id: dot.hpp,v 1.14 2007/05/16 14:58:57 cxsc Exp $ */

#ifndef dot_hpp_included
#define dot_hpp_included

#include <iostream>
#include <string> 

#include "compiler.h"
#include "RtsTyp.h"

namespace cxsc {


// ---- RTS - Definitionen ----
    

  // Form des dotprec. Akkus:
  //      Pointer auf a_btyp
  //  ist Pointer auf long
typedef  d_otpr   Dotprecision;

// ----------------------------
} // namespace cxsc 

#include "real.hpp"
#include "ioflags.hpp"

namespace cxsc {

class dotprecision;
class idotprecision;


//----------------------------------------------------------------------
// global verfgbare Dotprecision Variablen
//
//  dotakku[0..3] stehen fuer Matrix, Langzahl u.a. Pakete zur
//                Verfuegung
//  dotakku[4]    wird in den skalaren Paketen intern verwendet

#define MAXDOTAKKU      5
extern dotprecision dotakku[MAXDOTAKKU];


//!The Data Type dotprecision
/*!
The data types dotprecision, idotprecision, cdotprecision and cidotprecision are based on the scalar data types real, interval, complex, and
cinterval, respectively. Variables of these data types permit the exact representation of products of two arbitrary numbers of the
corresponding scalar base type and the exact summation of an arbitrary number of such products in a dotprecision accumulator,
i.e. in a fixed-point format of suitable size.

In general, values and variables of dotprecision data types occur during the calculation of scalar product expressions. Their values
are exactly representable in the format of the dotprecision data types, i.e. without rounding error, independent of the size of the
vectors or matrices contained in the scalar product expressions.
*/
class dotprecision
{
   private:
      // ---- Dataelements  -----------------------------------------
      Dotprecision akku;

   public:
      // ---- Constructors  -----------------------------------------
      //! Constructor of class dotprecision
      dotprecision(void)                 throw();
      //! Constructor of class dotprecision
      dotprecision(const dotprecision &) throw();
      

      //! Implementation of standard assigning operator
      dotprecision & operator =(const dotprecision &) throw();
      //! Implementation of standard assigning operator
      dotprecision & operator =(const real &)         throw();
      //! Implementation of standard assigning operator
      dotprecision & operator =(const l_real &)       throw(); // in l_real.cpp

      // ---- Typecasts     -----------------------------------------
      //! Constructor of class dotprecision
      explicit dotprecision(const real &)         throw();
      //! Constructor of class dotprecision
      explicit dotprecision(const l_real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa dotprecision(const real &)
      */
      friend inline dotprecision _dotprecision(const real &d) throw() { return dotprecision(d); }

      // ---- Destruktor    -----------------------------------------
      ~dotprecision(void) ;

      // ---- Input/Output  -----------------------------------------
      //! Implementation of standard output method
      friend std::string  & operator <<(std::string &,const dotprecision &) throw();
      //! Implementation of standard input method
      friend std::string  & operator >>(std::string &,dotprecision &)       throw();
      //! Implementation of standard input method
      friend void           operator >>(const std::string &,dotprecision &) throw();
      //! Implementation of standard input method
      friend void           operator >>(const char *,dotprecision &)        throw();
      //! Implementation of standard output method
      friend std::ostream & operator <<(std::ostream &,const dotprecision &) throw();
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &,dotprecision &)      throw();

      // ---- Std.Operators -----------------------------------------
      //! Implementation of standard algebraic negative sign operation
      friend  dotprecision operator -(const dotprecision &) throw();
      //! Implementation of standard algebraic positive sign operation
      friend  dotprecision operator +(const dotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend  dotprecision operator +(const dotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend  dotprecision operator -(const dotprecision &,const dotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline idotprecision operator |(const dotprecision &,const dotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend  dotprecision operator +(const dotprecision &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend  dotprecision operator +(const real &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend  dotprecision operator -(const dotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend  dotprecision operator -(const real &,const dotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline idotprecision operator |(const real &,const dotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline idotprecision operator |(const dotprecision &,const real &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend dotprecision & operator +=(dotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend dotprecision & operator -=(dotprecision &,const dotprecision &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend dotprecision & operator +=(dotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend dotprecision & operator -=(dotprecision &,const real &) throw();

      // ---- Comp.Operat. ------------------------------------------
      //! Implementation of standard negation operation
      friend bool operator  !(const dotprecision &) throw();
//      operator void *() const throw() { if(sign(*this)) return (void *)1; else return 0;}

      //! Implementation of standard equality operation
      friend bool operator ==(const dotprecision &,const dotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const dotprecision &,const dotprecision &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const dotprecision &,const dotprecision &) throw();
      //! Implementation of standard more-than operation
      friend bool operator  >(const dotprecision &,const dotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const dotprecision &,const dotprecision &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend bool operator >=(const dotprecision &,const dotprecision &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const real &,const dotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const real &,const dotprecision &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const real &,const dotprecision &) throw();
      //! Implementation of standard more-than operation
      friend bool operator  >(const real &,const dotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const real &,const dotprecision &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend bool operator >=(const real &,const dotprecision &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const dotprecision &,const real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const dotprecision &,const real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const dotprecision &,const real &) throw();
      //! Implementation of standard more-than operation
      friend bool operator  >(const dotprecision &,const real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const dotprecision &,const real &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend bool operator >=(const dotprecision &,const real &) throw();

      // ---- Others ------------------------------------------------
      //! Converting the exact dotprecision value with one rounding into a real value
      friend void rnd (const dotprecision&, real&, rndtype = RND_NEXT) throw();
      //! Converting the exact dotprecision value with one rounding into the nearest lower and upper real value
      friend void rnd (const dotprecision&, real&, real&)              throw();
      //! Converting the exact dotprecision value into an interval with the nearest lower and upper bound
      friend void rnd (const dotprecision&, interval&) throw(); // Blomquist
      //! Converting the exact dotprecision value with one rounding into a real value
      friend real rnd (const dotprecision&, rndtype = RND_NEXT)        throw();

      //! The absolute value of a dotprecision value
      friend dotprecision abs(const dotprecision &) throw();
      //! The sign of a dotprecision value
      friend int         sign(const dotprecision &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend dotprecision & accumulate (dotprecision&, const real&, const real&) throw();

      // ---- internal functions ------------------------------------
      inline Dotprecision* ptr() { return &akku; }
   private:
      dotprecision & negdot(void) throw(); // Negates current dotprecision
};

      void rnd (const dotprecision&, interval&) throw(); // Blomquist

} // namespace cxsc 

#endif

