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

/* CVS $Id: interval.cpp,v 1.22 2007/05/16 14:58:58 cxsc Exp $ */

#include "interval.hpp"
#include "ioflags.hpp"
#include "idot.hpp"
#include "dot.hpp"
#include "RtsFunc.h" // b_shr1 in mid()

// Auch fi_lib hat ein eigenes interval (wie RTS sein dotprecision)
//typedef struct fi_interval { double INF, SUP;} fi_interval;
// class interval;
#undef LINT_ARGS
#define CXSC_INCLUDE
#include <fi_lib.hpp>

namespace cxsc {

interval::interval(const idotprecision & a) throw()
{
   *this=rnd(a);
}

interval & interval::operator =(const idotprecision & a) throw()
{
   return *this=rnd(a);
}

interval::interval(const dotprecision & a) throw()
{
   rnd(a,inf,sup);
}

interval::interval(const dotprecision &a,const dotprecision &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
{
   inf=rnd(a,RND_DOWN);
   sup=rnd(b,RND_UP);
   if(inf>sup)
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("interval::interval(const dotprecision &,const dotprecision &)"));   
}

interval::interval(const l_real &a,const l_real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
{
   dotakku[0]=a;
   inf=rnd(dotakku[0],RND_DOWN);
   dotakku[0]=b;
   sup=rnd(dotakku[0],RND_UP);
   if(inf>sup)
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("interval::interval(const l_real &,const l_real &)"));   
}


interval & interval::operator =(const dotprecision & a) throw()
{
   rnd(a,inf,sup);
   return *this;
}

              
// ---- Standardfunkt ---- (arithmetische Operatoren)
/*
interval operator+(const interval &a, const interval &b) throw() { return interval(adddown(a.inf,b.inf),addup(a.sup,b.sup));  }
interval operator+(const interval &a, const real &b) throw() { return interval(adddown(a.inf,b),addup(a.sup,b));  }
interval operator+(const real &a, const interval &b) throw() { return interval(adddown(a,b.inf),addup(a,b.sup));  }

interval operator-(const interval &a, const interval &b) throw() { return interval(subdown(a.inf,b.inf),subup(a.sup,b.sup));  }
interval operator-(const interval &a, const real &b) throw() { return interval(subdown(a.inf,b),subup(a.sup,b));  }
interval operator-(const real &a, const interval &b) throw() { return interval(subdown(a,b.sup),subup(a,b.inf));  }

interval operator*(const interval &a, const interval &b) throw() { return mul_ii (a, b); }
interval operator*(const interval &a, const real &b) throw() { return mul_id (a, *(double*)&b); }
interval operator*(const real &a, const interval &b) throw() { return mul_di (*(double*)&a, b); }

interval operator/(const interval &a, const interval &b) throw() {  return div_ii (a, b); }
interval operator/(const interval &a, const real &b) throw() {  return div_id (a, *(double*)&b); }
interval operator/(const real &a, const interval &b) throw() {  return div_di (*(double*)&a, b); }
*/


bool operator ==(const interval &a,const dotprecision &r)     throw() {   return(r==a.inf && r==a.sup); }
bool operator !=(const interval &a,const dotprecision &r)     throw() {   return(r!=a.inf || r!=a.sup); }
bool operator ==(const dotprecision &r,const interval &a)     throw() {   return(r==a.inf && r==a.sup); }
bool operator !=(const dotprecision &r,const interval &a)     throw() {   return(r!=a.inf || r!=a.sup); }

bool operator <=(const dotprecision &a,const interval &b) throw()
{
   return(a>=b.inf && a<=b.sup);   
}
bool operator >=(const dotprecision &a,const interval &b) throw()
{
   return(a<=b.inf && a>=b.sup);   
}
bool operator <(const dotprecision &a,const interval &b) throw()
{
   return(a>b.inf && a<b.sup);   
}

bool operator <=(const interval &a,const dotprecision &b) throw()
{
   return(a.inf>=b && a.sup<=b);   
}
bool operator >=(const interval &a,const dotprecision &b) throw()
{
   return(a.inf<=b && a.sup>=b);   
}
bool operator >(const interval &a,const dotprecision &b) throw()
{
   return(a.inf<b && a.sup>b);   
}

interval operator+ (const interval& a, const interval& b) throw()
  { return interval (adddown(a.inf,b.inf),addup(a.sup,b.sup)); }
interval operator+ (const interval& a, const real& b) throw() 
  { return interval (adddown(a.inf,b),addup(a.sup,b)); }
interval operator+ (const real& a, const interval& b) throw() 
  { return interval (adddown(a,b.inf),addup(a,b.sup)); }

interval operator-  (const interval& a, const interval& b) throw() 
  { return interval ( subdown(a.inf,b.sup), subup(a.sup,b.inf)); }
interval operator-  (const interval& a, const real& b) throw() 
  { return interval ( subdown(a.inf,b), subup(a.sup,b)); }
interval operator-  (const real& a, const interval& b) throw() 
  { return interval ( subdown(a,b.sup), subup(a,b.inf)); }

  //------------------------------------------------------------------------
  //  a * b  = ueber Entscheidungstabelle :
  //
  //                      bi,bs >= 0       bi < 0, bs >=0       bi,bs < 0
  // ----------------+------------------+------------------+----------------
  // ai,as >= 0      I   ai*bi, as*bs   I   as*bi, as*bs   I   as*bi, ai*bs
  // ----------------+------------------+------------------+----------------
  // ai < 0, as >= 0 I   ai*bs, as*bs   I      ....        I   as*bi, ai*bi
  // ----------------+------------------+------------------+----------------
  // ai,as < 0       I   ai*bs, as*bi   I   ai*bs, ai*bi   I   as*bs, ai*bi
  // ----------------+------------------+------------------+----------------
  //
  //  .... :  min(ai*bs, as*bi), max(ai*ai, as*as)
  //
interval operator *(const interval& a, const interval& b) throw()
{
   interval tmp;

   if (sign(a.inf) >= 0) 
   {                                   // 1. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0)
      {                                  //     1. Spalte: [ai*bi, as*bs]
         tmp.inf = multdown(a.inf,b.inf);
         tmp.sup = multup(a.sup,b.sup);
      } else if (sign(b.sup) >= 0) 
      {                                  //     2. Spalte: [as*bi, as*bs]
         tmp.inf = multdown(a.sup,b.inf);
         tmp.sup = multup(a.sup,b.sup);
      } else 
      {                                  //     3. Spalte: [as*bi, ai*bs]
         tmp.inf = multdown(a.sup,b.inf);
         tmp.sup = multup(a.inf,b.sup);
      }
   } else if (sign(a.sup) >= 0) 
   {                                   // 2. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0) 
      {                                  //     1. Spalte: [ai*bs, as*bs]
         tmp.inf = multdown(a.inf,b.sup);
         tmp.sup = multup(a.sup,b.sup);
      } else if (sign(b.sup) >= 0) 
      {                                  //   2. Spalte: [min(ai*bs, as*bi),
         real hlp;                       //                 max(ai*ai, as*as)]

         tmp.inf = multdown(a.inf,b.sup);
         hlp     = multdown(a.sup,b.inf);

         if (hlp < tmp.inf) 
            tmp.inf = hlp;

         tmp.sup = multup(a.inf,b.inf);
         hlp     = multup(a.sup,b.sup);
         
         if (hlp > tmp.sup) 
            tmp.sup = hlp;               
      } else 
      {                                  //     3. Spalte: [as*bi, ai*bi]
         tmp.inf = multdown(a.sup,b.inf);
         tmp.sup = multup(a.inf,b.inf);
      }
   } else 
   {                                   // 3. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0) 
      {                                  //     1. Spalte: [ai*bs, as*bi]
         tmp.inf = multdown(a.inf,b.sup);
         tmp.sup = multup(a.sup,b.inf);
      } else if (sign(b.sup) >= 0) 
      {                                  //   2. Spalte: [ai*bs, ai*bi]
         tmp.inf = multdown(a.inf,b.sup);
         tmp.sup = multup(a.inf,b.inf);
      } else 
      {                                  //     3. Spalte: [as*bs, ai*bi]
         tmp.inf = multdown(a.sup,b.sup);
         tmp.sup = multup(a.inf,b.inf);
      }
   }
   return tmp;
}

  //------------------------------------------------------------------------
        //  a / b  = ueber Entscheidungstabelle :
  //
  //                      bi,bs > 0          bi,bs < 0
  // ----------------+------------------+------------------
  // ai,as >= 0      I   ai/bs, as/bi   I   as/bs, ai/bi
  // ----------------+------------------+------------------
  // ai < 0, as >= 0 I   ai/bi, as/bi   I   as/bs, ai/bs
  // ----------------+------------------+------------------
  // ai,as < 0       I   ai/bi, as/bs   I   as/bi, ai/bs
  // ----------------+------------------+------------------
  //
interval operator/  (const interval& a, const interval& b) throw(DIV_BY_ZERO)
{
   interval tmp;

   if ((sign(b.inf) <= 0) && (sign(b.sup) >= 0))
      cxscthrow(DIV_BY_ZERO("interval::interval operator/(const interval&,const interval&)"));   

   if (sign(a.inf) >= 0) 
   {                                    // 1. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) 
      {                                 //     1. Spalte: [ai/bs, as/bi]
         tmp.inf = divdown(a.inf,b.sup);
         tmp.sup = divup(a.sup,b.inf);
      } else 
      {                                 //     2. Spalte: [as/bs, ai/bi]
         tmp.inf = divdown(a.sup,b.sup);
         tmp.sup = divup(a.inf,b.inf);
      }
   } else if (sign(a.sup) >= 0) 
   {                                    // 2. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) 
      {                                 //     1. Spalte: [ai/bi, as/bi]
         tmp.inf = divdown(a.inf,b.inf);
         tmp.sup = divup(a.sup,b.inf);
      } else 
      {                                 //     2. Spalte: [as/bs, ai/bs]
         tmp.inf = divdown(a.sup,b.sup);
         tmp.sup = divup(a.inf,b.sup);
      }
   } else 
   {                                    // 3. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) 
      {                                 //     1. Spalte: [ai/bi, as/bs]
         tmp.inf = divdown(a.inf,b.inf);
         tmp.sup = divup(a.sup,b.sup);
      } else 
      {                                 //     2. Spalte: [as/bi, ai/bs]
         tmp.inf = divdown(a.sup,b.inf);
         tmp.sup = divup(a.inf,b.sup);
      }
   }

   return tmp;
}

interval operator*  (const real& a, const interval& b) throw()
{
   interval tmp;
   if (sign(a) == 0) 
   {
      tmp.inf = 0.0;
      tmp.sup = 0.0;      
   } else if (sign(a) > 0) 
   {
      tmp.inf = multdown(a,b.inf);
      tmp.sup = multup(a,b.sup);
   } else // if (sign(a) < 0) 
   {
      tmp.inf = multdown(a,b.sup);
      tmp.sup = multup(a,b.inf);
   }
   return tmp;
}

interval operator*  (const interval& a, const real& b) throw()
{
   interval tmp;
   if (sign(b) == 0) 
   {
      tmp.inf = 0.0;
      tmp.sup = 0.0;      
   } else if (sign(b) > 0) 
   {
      tmp.inf = multdown(a.inf,b);
      tmp.sup = multup(a.sup,b);
   } else // if (sign(b) < 0) 
   {
      tmp.inf = multdown(a.sup,b);
      tmp.sup = multup(a.inf,b);
   }
   return tmp;
}

interval operator/  (const real& a, const interval& b)  throw() { return (interval(a) / b); }
interval operator/  (const interval& a, const real& b)  throw() { return (a / interval(b)); }


// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const interval& a) throw()
{
   s << '['          << SaveOpt << RndDown
     << a.inf << ',' << RndUp 
     << a.sup        << RestoreOpt 
     << ']';
   return s;
}
std::string & operator << (std::string &s, const interval &a) throw()
{
   s+='[';
   s << SaveOpt << RndDown
     << a.inf;
   s+=',';
   s << RndUp 
     << a.sup   << RestoreOpt; 
   s+=']';
   return s;
}

std::istream & operator >> (std::istream &s, interval &a) throw()
{
   char c;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, '[');
   if (inpdotflag) 
      s.putback(c);

   s >> SaveOpt >> RndDown >> a.inf;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, ',');
   if (inpdotflag) s.putback(c);

   s >> RndUp >> a.sup >> RestoreOpt;

   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ']') 
         s.putback(c);
   }

   /*if (a.inf > a.sup) {
      errmon (ERR_INTERVAL(EMPTY));
   } */         
   return s;
}

std::string & operator >> (std::string &s, interval &a) throw()
{
   s = skipwhitespacessinglechar (s, '[');
   s >> SaveOpt >> RndDown >> a.inf;
   s = skipwhitespacessinglechar (s, ',');
   s >> RndUp >> a.sup >> RestoreOpt;
   s = skipwhitespaces (s);

   if (s[0] == ']') 
      s.erase(0,1);

    /*if (a.inf > a.sup) {
      errmon (ERR_INTERVAL(EMPTY));
    } */                         
   return s;
}

void operator >>(const std::string &s,interval &a) throw()
{
   std::string r(s);
   r>>a;
}
void operator >>(const char *s,interval &a) throw()
{
   std::string r(s);
   r>>a;
}

real mid(const interval & a) throw()
{
   dotakku[4] =  a.inf;
   dotakku[4] += a.sup;
    
   if (dotakku[4] != 0.0) 
   {
      // Division nur bei ungleich 0.0
      Dotprecision ptr = *dotakku[4].ptr();
                   
      // --------------------------------------------------------
      //  Dividiere dotakku[0] durch 2, mittels 1 Rechtsshift
                               
      ptr[(a_intg)++ptr[A_END]] = ZERO;
      
      b_shr1 (&ptr[(a_intg)ptr[A_BEGIN]], (a_intg)(ptr[A_END]-ptr[A_BEGIN]+1));
                                                   
      if (ptr[(a_intg)ptr[A_END]]   == ZERO) 
         --ptr[A_END];
                                                         
      if (ptr[(a_intg)ptr[A_BEGIN]] == ZERO) 
         ++ptr[A_BEGIN];   

      // --------------------------------------------------------
   }
          
   return rnd(dotakku[4]);
}

// for compatibility with CToolbox library - from former i_util.hpp
int in (const real& x, const interval& y )                  // Contained-in relation
  { return ( (Inf(y) <= x) && (x <= Sup(y)) ); }      //----------------------

int in ( const interval& x, const interval& y )       // Contained-in relation
{                                                     //----------------------
  return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) );
}

// function in dot.hpp and dot.cpp
//void rnd ( const dotprecision& d, interval& x )    // Round dotprecision to interval
//{                                            //-------------------------------
//  real Lower, Upper;
//
//  rnd(d,Lower,Upper);         // Rounding to interval bounds
//  x = interval(Lower,Upper);
//}

/*!
\param x The interval for which the epsilon inflation should be computed
\param eps The real value of epsilon
\return The inflated interval

The epsilon inflation or \f$ \epsilon \f$-inflation of a real floating-point interval \f$ \left[ x \right] \in R \f$ is defined by

\f[
[x] \bowtie \varepsilon := \left\{
\matrix{
	{[x] + [-\varepsilon , +\varepsilon ] \cdot d \left( [x] \right)} & {\mbox{if } d \left( [x] \right) \not= 0} \hfill \cr
	{[x] + [-x_{min} , +x_{min} ]} & {\mbox{otherwise}} \hfill
}
\right.
\f]

The \f$ \epsilon \f$-inflation is defined componentwise for complex intervals and interval vectors and matrices.
*/
interval Blow ( const interval& x, const real& eps )    // Epsilon inflation
{                                                       //------------------
  interval tmp;
  tmp =  (1.0+eps)*x - eps*x;
  return interval(pred(Inf(tmp)),succ(Sup(tmp)));
}

int Disjoint ( const interval& a, const interval& b )   // Test for disjointedness
{                                                       //------------------------
  return (Inf(a) > Sup(b)) || (Inf(b) > Sup(a));
}

real AbsMin ( const interval& x )                       // Absolute minimum of
{                                                       // an interval
  if ( in(0.0,x) )                                      //--------------------
    return 0.0;
  else if (Inf(x) > 0.0)
    return Inf(x);
  else
    return -Sup(x);
}

real AbsMax ( const interval& x )                       // Absolute maximum of
{                                                       // an interval
  real a, b;                                            //--------------------

  a = abs(Inf(x));
  b = abs(Sup(x));

  if (a > b)
    return a;
  else
    return b;
}

real RelDiam ( const interval& x )                               // Relative diameter
{                                                         // of an interval
  if ( in(0.0,x) )                                        //------------------
    return diam(x);
  else
    return( divu(diam(x),AbsMin(x)) );
}

//----------------------------------------------------------------------------
// Checks if the diameter of the interval 'x' is less or equal to 'n' ulps.
// Ulp is an abbreviation for units in the last place.
//----------------------------------------------------------------------------
/*!
\param x Interval to be checked
\param n Numbers of ulp's
\return Integer representing the result of the check

Ulp is an abbreviation for units in the last place.
*/
int UlpAcc ( const interval& x, int n )
{
  int  i;
  real Infimum;

  Infimum = Inf(x);
  for (i = 1; i <= n; i++) Infimum = succ(Infimum);
  return (Infimum >= Sup(x));
}

// The following interval constants are optimal inclusions:

    const real Pi_Inf = 7074237752028440.0 / 2251799813685248.0;
    const interval Pi_interval = interval(Pi_Inf,succ(Pi_Inf));     // Pi

    const real Pi2_Inf = 7074237752028440.0 / 1125899906842624.0;
    const interval Pi2_interval = interval(Pi2_Inf,succ(Pi2_Inf));  // 2*Pi

    const real Pi3_Inf = 5305678314021330.0 / 562949953421312.0;
    const interval Pi3_interval = interval(Pi3_Inf,succ(Pi3_Inf));  // 3*Pi

    const real Pid2_Inf = 7074237752028440.0 / 4503599627370496.0;
    const interval Pid2_interval = interval(Pid2_Inf,succ(Pid2_Inf));  // Pi/2

    const real Pid3_Inf = 4716158501352293.0 / 4503599627370496.0;
    const interval Pid3_interval = interval(Pid3_Inf,succ(Pid3_Inf));  // Pi/3

    const real Pid4_Inf = 7074237752028440.0 / 9007199254740992.0;
    const interval Pid4_interval = interval(Pid4_Inf,succ(Pid4_Inf));  // Pi/4

    const real Pir_Inf = 5734161139222658.0 / 18014398509481984.0;
    const interval Pir_interval = interval(Pir_Inf,succ(Pir_Inf));     // 1/Pi

    const real Pi2r_Inf = 5734161139222658.0 / 36028797018963968.0;
    const interval Pi2r_interval = interval(Pi2r_Inf,succ(Pi2r_Inf)); 
                                                                   // 1/(2*Pi)

    const real Pip2_Inf = 5556093337880030.0 / 562949953421312.0;
    const interval Pip2_interval = interval(Pip2_Inf,succ(Pip2_Inf));  // Pi^2

    const real SqrtPi_Inf = 7982422502469482.0 / 4503599627370496.0;
    const interval SqrtPi_interval = interval(SqrtPi_Inf,succ(SqrtPi_Inf));  
                                                                // sqrt(Pi)

    const real Sqrt2Pi_Inf = 5644425081792261.0 / 2251799813685248.0;
    const interval Sqrt2Pi_interval = interval(Sqrt2Pi_Inf,succ(Sqrt2Pi_Inf));
                                                               // sqrt(2Pi)

    const real SqrtPir_Inf = 5081767996463981.0 / 9007199254740992.0;
    const interval SqrtPir_interval = interval(SqrtPir_Inf,succ(SqrtPir_Inf));
                                                               // 1/sqrt(Pi)

    const real Sqrt2Pir_Inf = 7186705221432912.0 / 18014398509481984.0;
    const interval Sqrt2Pir_interval=interval(Sqrt2Pir_Inf,succ(Sqrt2Pir_Inf));
                                                              // 1/sqrt(2Pi)

    const real Sqrt2_Inf = 6369051672525772.0 / 4503599627370496.0;
    const interval Sqrt2_interval=interval(Sqrt2_Inf,succ(Sqrt2_Inf));
                                                                  // sqrt(2)

    const real Sqrt2r_Inf = 6369051672525772.0 / 9007199254740992.0;
    const interval Sqrt2r_interval=interval(Sqrt2r_Inf,succ(Sqrt2r_Inf));
                                                                  // 1/sqrt(2)

    const real Sqrt3_Inf = 7800463371553962.0 / 4503599627370496.0;
    const interval Sqrt3_interval=interval(Sqrt3_Inf,succ(Sqrt3_Inf));
                                                                  // sqrt(3)

    const real Sqrt3d2_Inf = 7800463371553962.0 / 9007199254740992.0;
    const interval Sqrt3d2_interval=interval(Sqrt3d2_Inf,succ(Sqrt3d2_Inf));
                                                                 // sqrt(3)/2

    const real Sqrt3r_Inf = 5200308914369308.0 / 9007199254740992.0;
    const interval Sqrt3r_interval=interval(Sqrt3r_Inf,succ(Sqrt3r_Inf));
                                                                  // 1/sqrt(3)

    const real Ln2_Inf = 6243314768165359.0 / 9007199254740992.0;
    const interval Ln2_interval=interval(Ln2_Inf,succ(Ln2_Inf)); // ln(2)

    const real Ln2r_Inf = 6497320848556798.0 / 4503599627370496.0;
    const interval Ln2r_interval=interval(Ln2r_Inf,succ(Ln2r_Inf)); // 1/ln(2)

    const real Ln10_Inf = 5184960683398421.0 / 2251799813685248.0;
    const interval Ln10_interval=interval(Ln10_Inf,succ(Ln10_Inf)); // ln(10)

    const real Ln10r_Inf = 7823553867474190.0 / 18014398509481984.0;
    const interval Ln10r_interval=interval(Ln10r_Inf,succ(Ln10r_Inf)); 
                                                                   // 1/ln(10)

    const real LnPi_Inf = 5155405087351229.0 / 4503599627370496.0;
    const interval LnPi_interval=interval(LnPi_Inf,succ(LnPi_Inf)); // ln(Pi)

    const real Ln2Pi_Inf = 8277062471433908.0 / 4503599627370496.0;
    const interval Ln2Pi_interval=interval(Ln2Pi_Inf,succ(Ln2Pi_Inf)); 
                                                                    // ln(2Pi)

    const real E_Inf = 6121026514868073.0 / 2251799813685248.0;
    const interval E_interval=interval(E_Inf,succ(E_Inf)); // e

    const real Er_Inf = 6627126856707895.0 / 18014398509481984.0;
    const interval Er_interval=interval(Er_Inf,succ(Er_Inf)); // 1/e

    const real Ep2_Inf = 8319337573440941.0 / 1125899906842624.0;
    const interval Ep2_interval=interval(Ep2_Inf,succ(Ep2_Inf)); // e^2

    const real Ep2r_Inf = 4875967449235915.0 / 36028797018963968.0;
    const interval Ep2r_interval=interval(Ep2r_Inf,succ(Ep2r_Inf)); // 1/e^2

    const real EpPi_Inf = 6513525919879993.0 / 281474976710656.0;
    const interval EpPi_interval=interval(EpPi_Inf,succ(EpPi_Inf)); // e^(Pi)

    const real Ep2Pi_Inf = 4710234414611993.0 / 8796093022208.0;
    const interval Ep2Pi_interval=interval(Ep2Pi_Inf,succ(Ep2Pi_Inf)); 
                                                                   // e^(2Pi)

    const real EpPid2_Inf = 5416116035097439.0 / 1125899906842624.0;
    const interval EpPid2_interval=interval(EpPid2_Inf,succ(EpPid2_Inf)); 
                                                                   // e^(Pi/2)

    const real EpPid4_Inf = 4938827609611434.0 / 2251799813685248.0;
    const interval EpPid4_interval=interval(EpPid4_Inf,succ(EpPid4_Inf)); 
                                                                   // e^(Pi/4)


} // namespace cxsc

