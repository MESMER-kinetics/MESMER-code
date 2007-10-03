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

/* CVS $Id: complex.cpp,v 1.14 2007/05/16 14:58:56 cxsc Exp $ */
/* New versions of product(...), quotient(...); Blomquist, 26.10.02; */

#include "complex.hpp"
#include "dot.hpp"
#include "cdot.hpp"
#include "real.hpp"
#include "rmath.hpp"
#include "idot.hpp"
#include "interval.hpp"

#include "cinterval.hpp"

namespace cxsc {

complex::complex(const cdotprecision & a) throw()
{
   *this=rnd(a);
}

complex & complex::operator =(const cdotprecision & a) throw()
{
   return *this=rnd(a);
}

bool operator== (const complex & a, const dotprecision & b)    throw() { return !a.im && a.re==b; }
bool operator== (const dotprecision & a, const complex & b)    throw() { return !b.im && a==b.re; }
bool operator!= (const complex & a, const dotprecision & b)    throw() { return !!a.im || a.re!=b; }
bool operator!= (const dotprecision & a, const complex & b)    throw() { return !!b.im || a!=b.re; }

complex operator *(const complex &a,const complex &b) throw()
{
   complex tmp;
    
   dotakku[4] = 0.0;
   accumulate (dotakku[4],  a.re, b.re);
   accumulate (dotakku[4], -a.im, b.im);
   rnd (dotakku[4], tmp.re, RND_NEXT);

   dotakku[4] = 0.0;
   accumulate (dotakku[4], a.re, b.im);
   accumulate (dotakku[4], a.im, b.re);
   rnd (dotakku[4], tmp.im, RND_NEXT);

   return tmp;
}

complex muld(const complex &a, const complex &b) throw()
{  // Blomquist 07.11.02;
   complex tmp;
    
   dotakku[4] = 0.0;
   accumulate (dotakku[4],  a.re, b.re);
   accumulate (dotakku[4], -a.im, b.im);
   rnd (dotakku[4], tmp.re, RND_DOWN);

   dotakku[4] = 0.0;
   accumulate (dotakku[4], a.re, b.im);
   accumulate (dotakku[4], a.im, b.re);
   rnd (dotakku[4], tmp.im, RND_DOWN);

   return tmp;
}

complex mulu(const complex &a, const complex &b) throw()
{  // Blomquist 07.11.02;
   complex tmp;
    
   dotakku[4] = 0.0;
   accumulate (dotakku[4],  a.re, b.re);
   accumulate (dotakku[4], -a.im, b.im);
   rnd (dotakku[4], tmp.re, RND_UP);

   dotakku[4] = 0.0;
   accumulate (dotakku[4], a.re, b.im);
   accumulate (dotakku[4], a.im, b.re);
   rnd (dotakku[4], tmp.im, RND_UP);

   return tmp;
}


static const int Min_Exp_ = 1074, minexpom = -914, 
                 maxexpo1 = 1022, MANT_W   = 52;

void product(real a, real b, real c, real d,
             int& overfl, real& p1, interval& p2)
// New version of function product(...) from Blomquist, 26.10.02;
// Input data: a,b,c,d;  Output data: overfl, p1, p2;
// In case of overfl=0 the interval p1+p2 is an inclusion of a*b + c*d;
// overfl=1 (overflow) signalizes that p1+p2 must be multiplied with 2^1074
// to be an inclusion of a*b + c*d;
// overfl=-1 (underflow) signalizes that p1+p2 must be multiplied with 
// 2^-1074 to be an inclusion of a*b + c*d; 

{  
   int exa, exb, exc, exd;  // Exp. von a-d
   dotprecision dot;
   int inexact;

   overfl  = 0;
   inexact = 0; // false

   dot = 0.0;
   exa = expo(a);
   exb = expo(b);
   exc = expo(c);
   exd = expo(d);

   if ( sign(a) == 0  ||   sign(b) == 0 )    //  a * b == 0.0
      if ( sign(c) == 0  ||  sign(d) == 0 )  //  a * b == c * d == 0
                                             //  dot := #(0);
         ;    // No Operation necessary
      else 
      {
         //  a * b == 0;  c * d != 0;
         if (exc+exd > maxexpo1) 
         {
            //  overflow !
	     if ( exc > exd ) c = comp( mant(c), exc-Min_Exp_ );
	     else d = comp( mant(d), exd-Min_Exp_ );
             overfl = 1;
         } else 
         if  ( exc+exd < minexpom ) 
         { // undeflow; 
           // c = comp( mant(c), exc+Min_Exp_ ); kann Overflow erzeugen!
	     if (exc < exd) c = comp( mant(c), exc+Min_Exp_ );
	     else d = comp( mant(d), exd+Min_Exp_ );
             overfl = -1;
         }
         accumulate(dot,c,d);
      } else //  a,b != 0
      if ( sign(c) == 0  ||  sign(d) == 0 ) 
      {
         //  a*b != 0, c * d == 0
         if (exa+exb > maxexpo1) 
         {
            //  overflow !
            if ( exa > exb ) a = comp( mant(a), exa-Min_Exp_ );
	    else b = comp( mant(b), exb-Min_Exp_ );
            overfl = 1;
         } else 
         if (exa+exb < minexpom) 
         { // undeflow; 
           // a = comp( mant(a), exa+Min_Exp_ ); kann Overflow erzeugen!
	     if (exa < exb) a = comp( mant(a), exa+Min_Exp_ );
	     else b = comp( mant(b), exb+Min_Exp_ );
            overfl = -1;
         }
         accumulate(dot,a,b);
      } else 
      {
         // a,b,c,d != 0
         if (exa+exb > maxexpo1) 
         {  //  overflow bei a*b
            if ( exa > exb ) a = comp( mant(a), exa-Min_Exp_ );
	    else b = comp( mant(b), exb-Min_Exp_ );
            if (exc > MANT_W) c = comp( mant(c), exc-Min_Exp_ );
            else if (exd > MANT_W)
               d = comp( mant(d), exd-Min_Exp_ );
            else 
            {
               // underflow wegen Skalierung bei c*d
                c = 0.0;
                inexact = 1; // true
            }
	    overfl = 1;  // Hat vorher gefehlt!! Blomquist, 24.10.02;
         } else if (exc+exd > maxexpo1) 
         {
            // overflow bei c*d
            if ( exc > exd ) c = comp( mant(c), exc-Min_Exp_ );
	    else d = comp( mant(d), exd-Min_Exp_ );
            if (exa > MANT_W) a = comp( mant(a), exa-Min_Exp_ );
            else if (exb > MANT_W)
               b = comp( mant(b), exb-Min_Exp_ );
            else 
            {
               // underflow wegen Skalierung bei a*b
               a = 0.0;
               inexact = 1; // true
            }
            overfl = 1;
         } else 
         if ( exa+exb < minexpom  &&  exc+exd < minexpom ) 
         {
            //  underflow bei a*b und bei c*d
	     if (exa < exb) a = comp( mant(a), exa+Min_Exp_ );
	     else b = comp( mant(b), exb+Min_Exp_ );
	     if (exc < exd) c = comp( mant(c), exc+Min_Exp_ );
	     else d = comp( mant(d), exd+Min_Exp_ );
             overfl = -1;
         }
         accumulate(dot, a, b);
         accumulate(dot, c, d);
      }

   p1 = rnd(dot);
   dot -= p1;
   rnd(dot,p2);  // Blomquist, 07.11.02;

   if (inexact)
       p2 = _interval( pred(Inf(p2)), succ(Sup(p2)) );

} // end product

real quotient (real z1, interval z2, real n1, 
               interval n2, int round, int zoverfl, int noverfl)
// z1+z2 is an inclusion of a numerator.
// n1+n2 is an inclusion of a denominator.
// quotient(...) calculates with q1 an approximation of (z1+z2)/(n1+n2)
// using staggered arithmetic.
// Rounding with round (-1,0,+1) is considered.
// zoverfl and noverfl are considered by scaling back with 2^1074 or 2^-1074
// n1+n2 > 0 is assumed.
{
   real q1=0, scale;
   interval q2, nh;
   idotprecision id;
   int vorz, anz_scale, ex = 0;

   vorz = sign(z1); // n1,n2 > 0 is assumed!!
                    // so the sign of q1 ist sign of z1.
   if ( zoverfl == -1  &&  noverfl == 1 ) 
   {
      //  result in the undeflow range:
      switch (round) 
      {
         case RND_DOWN:
            if (vorz >= 0) 
                q1 = 0.0;
            else         
		q1 = -minreal; // Blomquist: MinReal --> minreal;
            break;
         case RND_NEXT:
            q1 = 0.0;
            break;
         case RND_UP:
            if (vorz <= 0) 
		q1 = 0.0;
            else         
                q1 = minreal;  // Blomquist: MinReal --> minreal;
            break;
      } // switch
   } else if ( zoverfl==1  &&  noverfl==-1 ) 
   {  
      //  result in the overflow range:
      if (vorz >= 0) q1 = MaxReal+MaxReal;  // Overflow
      else q1 = -MaxReal-MaxReal;           // Overflow
   } else 
   {  
      q1 = divd(z1, n1);  // down, to get q2 >= 0
      nh = _interval(addd(n1, Inf(n2)),
                   addu(n1, Sup(n2)));

      // q2:= ##( z1 - q1*n1 + z2 - q1*n2 );
      id = z1;
      accumulate(id, -q1, n1);
      id += z2;
      accumulate(id, -q1, n2);
      q2 = rnd(id);

      switch (round) // Considering the rounding before scaling
      {
         case RND_DOWN:
            q1 = adddown(q1, divd(Inf(q2), Sup(nh)));
            break;
         case RND_NEXT:
            q1 = q1 + (Inf(q2)+Sup(q2))*0.5/n1;
            break;
         case RND_UP:
            q1 = addup(q1, divu(Sup(q2), Inf(nh)));
            break;
      } // switch


      //  scaling back, if numerator|denominator - over-|underflow :

      //  actuell as follows:
      //  q1:= comp( mant(q1), expo(q1) + (zoverfl-noverfl)*1074 );
      //  The scaling with 2^1074 must be done in two steps with 2^ex
      //  and with the factor scale:

      anz_scale = zoverfl - noverfl; // |anz_scale| <= 1; 
      if (anz_scale > 0)
      {
	  scale = comp(0.5, +1024);
	  ex = MANT_W-1;
      }
      else if (anz_scale < 0)
      {
	  scale = comp(0.5, -1022);
	  ex = -MANT_W+1;
      }
      else scale = 1.0;  // ex = 0 is already initialized for this case 


      if (ex) times2pown(q1,ex); // EXACT part scaling, if ex != 0.
      switch (round) // correct rounding with the second factor scale:
      {
	  case RND_DOWN:
              q1 = multdown(q1, scale);
	      break;
          case RND_NEXT:
              q1 = q1 * scale;
              break;
          case RND_UP:
              q1 = multup(q1, scale);
              break;
      }  // switch
   }
   return q1;
} // end of quotient

complex _c_division(complex a, complex b, int round) throw(DIV_BY_ZERO)
{
    if (0.0 == (sqr(Re(b))+sqr(Im(b)))) {
      cxscthrow(DIV_BY_ZERO("complex operator / (const complex&,const complex&)"));
    }      
    
   int zoverflow, noverflow;
   real z1, n1;
   interval z2, n2;
   complex tmp;

   product (Re(b), Re(b), Im(b), Im(b), noverflow, n1, n2);
   product (Re(a), Re(b), Im(a), Im(b), zoverflow, z1, z2);
   SetRe (tmp, quotient (z1, z2, n1, n2, round, zoverflow, noverflow));
   product (Im(a), Re(b), -Re(a), Im(b), zoverflow, z1, z2);
   SetIm (tmp, quotient (z1, z2, n1, n2, round, zoverflow, noverflow));
   return tmp;
}

complex divn (const complex & a, const complex & b)  
{  // Blomquist: vorher c_divd(...), 07.11.02;
   return _c_division(a, b, RND_NEXT);
}

complex divd (const complex & a, const complex & b) 
{  // Blomquist: vorher c_divd(...), 07.11.02;
   return _c_division(a, b, RND_DOWN);
}

complex divu (const complex & a, const complex & b)  
{  // Blomquist: vorher c_divu(...), 07.11.02;
   return _c_division(a, b, RND_UP);
}

complex operator / (const complex &a,const complex &b) throw()
{
   return divn(a,b);
}

real abs2(const complex &a) throw()
{
   dotakku[4]=0;
   accumulate(dotakku[4],a.re,a.re);
   accumulate(dotakku[4],a.im,a.im);
   return rnd(dotakku[4]);
}

real abs (complex z) throw()
{   //  calculation of |z|; Blomquist 06.12.02;
    return sqrtx2y2(Re(z),Im(z));
}


// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const complex& a) throw()
{
   s << '('          
     << a.re << ','  
     << a.im       
     << ')';
   return s;
}
std::string & operator << (std::string &s, const complex &a) throw()
{
   s+='(';
   s << a.re;
   s+=',';
   s << a.im; 
   s+=')';
   return s;
}

std::istream & operator >> (std::istream &s, complex &a) throw()
{
   char c;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, '(');
   if (inpdotflag) 
      s.putback(c);

   s >> a.re;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, ',');
   if (inpdotflag) s.putback(c);

   s >> a.im >> RestoreOpt;

   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ')') 
         s.putback(c);
   }

   /*if (a.re > a.im) {
      errmon (ERR_INTERVAL(EMPTY));
   } */         
   return s;
}

std::string & operator >> (std::string &s, complex &a) throw()
{
   s = skipwhitespacessinglechar (s, '(');
   s >> SaveOpt >> RndDown >> a.re;
   s = skipwhitespacessinglechar (s, ',');
   s >> RndUp >> a.im >> RestoreOpt;
   s = skipwhitespaces (s);

   if (s[0] == ')') 
      s.erase(0,1);

    /*if (a.re > a.im) {
      errmon (ERR_INTERVAL(EMPTY));
    } */                         
   return s;
}

void operator >>(const std::string &s,complex &a) throw()
{
   std::string r(s);
   r>>a;
}
void operator >>(const char *s,complex &a) throw()
{
   std::string r(s);
   r>>a;
}

complex exp(const complex& x) { return mid(exp(cinterval(x))); }
complex cos(const complex& x) { return mid(cos(cinterval(x))); }
complex sin(const complex& x) { return mid(sin(cinterval(x))); }
complex cosh(const complex& x) { return mid(cosh(cinterval(x))); }
complex sinh(const complex& x) { return mid(sinh(cinterval(x))); }

complex tan(const complex& x) throw() { return mid(tan(cinterval(x))); }
complex cot(const complex& x) { return mid(cot(cinterval(x))); }
complex tanh(const complex& x) { return mid(tanh(cinterval(x))); }
complex coth(const complex& x) { return mid(coth(cinterval(x))); }

real arg(const complex& x) { return mid(arg(cinterval(x))); }
real Arg(const complex& x) throw() { return mid(Arg(cinterval(x))); }

complex ln(const complex& x) throw() { return mid(ln(cinterval(x))); }
complex Ln(const complex& x) throw() { return mid(Ln(cinterval(x))); }

complex sqr(const complex& x) { return mid(sqr(cinterval(x))); }

complex sqrt(const complex& x) throw() { return mid(sqrt(cinterval(x))); }
complex sqrt(const complex& x, int d) throw() { return mid(sqrt(cinterval(x),d)); }

complex power_fast(const complex& x,int d) throw() { return mid(power_fast(cinterval(x),d)); }
complex power(const complex& x,int d) throw() { return mid(power(cinterval(x),d)); }
complex pow(const complex& x, const real& r) throw() { return mid(pow(cinterval(x),interval(r))); }
complex pow(const complex& x, const complex& y) throw() { return mid(pow(cinterval(x),cinterval(y))); }

complex asin(const complex& x) throw() { return mid(asin(cinterval(x))); }
complex acos(const complex& x) { return mid(acos(cinterval(x))); }
complex asinh(const complex& x) { return mid(asinh(cinterval(x))); }
complex acosh(const complex& x) throw() { return mid(acosh(cinterval(x))); }
complex atan(const complex& x) throw() { return mid(atan(cinterval(x))); }
complex acot(const complex& x) throw() { return mid(acot(cinterval(x))); }
complex atanh(const complex& x) { return mid(atanh(cinterval(x))); }
complex acoth(const complex& x) { return mid(acoth(cinterval(x))); }

} // namespace cxsc

