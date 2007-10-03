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

/* CVS $Id: imath.cpp,v 1.25 2007/05/16 14:58:58 cxsc Exp $ */

#include "imath.hpp"
#include "rmath.hpp"

// Auch fi_lib hat ein eigenes interval (wie RTS sein dotprecision)
//typedef struct fi_interval { double INF, SUP;} fi_interval;
#undef LINT_ARGS
#define CXSC_INCLUDE
#include <fi_lib.hpp>

extern "C" {
#ifndef rfcth_included
#define rfcth_included
#include "r_fcth.h"
#endif
}

namespace cxsc {


interval sqr    (const interval &a) throw() 
{ 
  interval res;
  res= a*a;
  if (Inf(res)<0) Inf(res)=0;
  return res;
}

interval sqrt   (const interval &a)         { return j_sqrt(a); }

interval sqrt   (const interval &a, int n)  throw(STD_FKT_OUT_OF_DEF)
{ 
   if(Inf(a)>=0.0) 
      return pow(a,interval(1.0,1.0)/n); 
   else { 
      cxscthrow(STD_FKT_OUT_OF_DEF("interval sqrt (const interval &a, int n)"));
      return interval(-1.0); // dummy result
   }
}

interval sqrt1px2(const interval& x) throw()
// Inclusion of sqrt(1+x^2); Blomquist, 13.12.02;
{
    interval t = abs(x),y;
    if (expo(Inf(t)) > 33)
    {
	y = t;
	Sup(y) = succ(Sup(y));
    } else if (expo(Sup(t)) > 33)
    {
	y = interval(Inf(t));  // interval(Inf(t)) is a point interval!
	y = sqrt(1+y*y);       // --->  y*y == sqr(y);
	y = interval(Inf(y),succ(Sup(t)));
     } else y = sqrt(1+sqr(t));
    return y;
}

interval sqrtx2y2(const interval& x, const interval& y) throw()
// Inclusion of sqrt(x^2+y^2); Blomquist, 13.12.02;
{
    interval a=abs(x), b=abs(y), r;
    int exa=expo(Sup(a)), exb=expo(Sup(b)), ex;
    if (exb > exa)
    {  // Permutation of a,b:
	r = a;  a = b;  b = r;
	ex = exa;  exa = exb;  exb = ex;
    }
    ex = 511 - exa;
    times2pown(a,ex);
    times2pown(b,ex);
    r = sqrt(a*a + b*b);
    times2pown(r,-ex);
    return r;
} // sqrtx2y2

//********************************************************************
//  Constants for:   interval sqrtp1m1(const interval& x) throw()
//  Blomquist 05.08.03
const real Delta_f = 2*minreal;
const real q_sqrtp1m1m = 9007199254740984.0 / 9007199254740992.0;
const real q_sqrtp1m1p = 4503599627370502.0 / 4503599627370496.0;
//********************************************************************
interval sqrtp1m1(const interval& x) throw()
// interval(a,b) is an inclusion of sqrt(x+1)-1;
// Exported by imath.hpp;       Blomquist, 05.08.03;
{
    real a=0,b=0,ix=Inf(x),sx=Sup(x);
    int ex_ix,ex_sx,sgn_ix,sgn_sx;

    // Calculating the lower bound:
    ex_ix = expo(ix);  sgn_ix = sign(ix);
    if (ex_ix<=-1021) //  ==> |ix| < 2^(-1021)
        { if (sgn_ix) a = sqrtp1m1(ix) - Delta_f; }
    else // |ix| >= 2^(-1021)
	if (ex_ix<=-51) 
	{
	    times2pown(ix,-1); // exact division by 2
	    a = pred(ix);
	}
	else 
	    if (sgn_ix>0) a = (ix>0.67) ? 
	       Inf( sqrt(interval(ix)+1)-1 ) : sqrtp1m1(ix)*q_sqrtp1m1m;
	    else a = (ix<-0.25) ? 
	       Inf( sqrt(interval(ix)+1)-1 ) : sqrtp1m1(ix)*q_sqrtp1m1p;
    // Calculating the upper bound:
    if (ix == sx) { ex_sx = ex_ix;  sgn_sx = sgn_ix; }
    else { ex_sx = expo(sx);  sgn_sx = sign(sx); }
    if (ex_sx<=-1021) //  ==> |sx| < 2^(-1021)
        { if (sgn_sx) b = sqrtp1m1(sx) + Delta_f; }
    else // |sx| >= 2^(-1021)
	if (ex_sx<=-47) { b = sx; times2pown(b,-1); }
	else 
	    if (sgn_sx>0) b = (sx>0.58) ? 
	        Sup( sqrt(interval(sx)+1)-1 ) : sqrtp1m1(sx)*q_sqrtp1m1p;
	    else b = (sx<-0.32) ? 
	         Sup( sqrt(interval(sx)+1)-1 ) : sqrtp1m1(sx)*q_sqrtp1m1m;
    return interval(a,b);
} // sqrtp1m1()

//********************************************************************
// Konstanten fr die Intervallfunktion sqrt(x^2-1):
real q_sqrtx2m1p(4503599627370501.0/4503599627370496.0);  // (1+e(f))
real q_sqrtx2m1m(9007199254740986.0/9007199254740992.0);  // (1-e(f))
//********************************************************************

interval sqrtx2m1(const interval& x)
// sqrt(x^2-1);  Blomquist, 13.04.04;
{
    interval z = abs(x);
    real r1,r2;
    r1 = sqrtx2m1(Inf(z)) * q_sqrtx2m1m;
    r2 = Sup(z); 
    if (expo(r2)<26)  
        r2 = sqrtx2m1(r2) * q_sqrtx2m1p; 
    // expo(r2) >= 26 --> r2 = Sup(z) ist optimale Oberschranke!
    return interval(r1,r2);
} // sqrtx2m1 (Intervallargumente)


//********************************************************************
// Konstanten fr die Intervallfunktion sqrt(1-x2):
real q_sqrt1mx2p(4503599627370501.0/4503599627370496.0);  // (1+e(f))
real q_sqrt1mx2m(9007199254740985.0/9007199254740992.0);  // (1-e(f))
//********************************************************************

interval sqrt1mx2(const interval& x)
// sqrt(1-x2);  Blomquist, 13.04.04;
{
    interval z = abs(x);  // sqrt(1-t2) monoton fallend in t aus z;
    real r1,r2,sz,iz;
    sz = Sup(z);  iz = Inf(z);
    // Berechnung der Unterschranke r1:
    r2 = sqrt1mx2(sz);
    r1 = (sz==0)? 1 : r2 * q_sqrt1mx2m;
    // Berechnung der Oberrschranke r2:
    if (iz<4.81e-8) r2 = 1; // r2=1 ist immer korrekte Oberschranke!
    else r2 = (sz==iz)? r2*q_sqrt1mx2p : sqrt1mx2(iz)*q_sqrt1mx2p;
    return interval(r1,r2);

} // sqrtx2m1 (Intervallargumente)


// Konstanten fuer die Intervallfunktion exp(-x^2):
const real q_exp_x2p = 4503599627370502.0 / 4503599627370496.0;  // (1+e(f))
const real q_exp_x2m = 9007199254740984.0 / 9007199254740992.0;  // (1-e(f))
  // Oberschranke expmx2_UB fuer |x| > expmx2_x0:  
const real expmx2_UB = 2.225073858507447856659E-308;
const real expmx2_x0 = 7491658466053896.0 / 281474976710656.0;

interval expmx2(const interval& x)
// e^(-x^2);  Blomquist, 05.07.04;
{
    real y,r1,r2,Sz,Iz;
    interval z = abs(x); 
    // Berechnung einer Unterschranke:
    Sz = Sup(z);  Iz = Inf(z);
    y = expmx2(Sz);
    r1 = y;
    if (Sz<=expmx2_x0) r1 = r1 * q_exp_x2m; // Abrunden
    if (Sz==0) r1 = y;
    // Berechnung einer Oberschranke:
    if (Iz>expmx2_x0) r2 = expmx2_UB;
    else r2 = (Sz==Iz)? y * q_exp_x2p : expmx2(Iz) * q_exp_x2p; 
    if (r2>1) r2 = 1.0;
    return interval(r1,r2);

} // expmx2 (Intervallargumente)


interval expm1(const interval& x)
// Interval function of exp(x)-1
{
    real y,r1,r2,Sx,Ix;
    Sx = Sup(x);  Ix = Inf(x);
    y = expm1(Ix);
    // calculation of a lower bound r1:
    r1 = (y>0)? y*q_exmm : y*q_exmp;  // rounding downwards
    if (r1<-1) r1 = -1;
    // calculation of an upper bound r2:
    if (Sx!=Ix) y = expm1(Sx);
    r2 = (y>0)? y*q_exmp : y*q_exmm;

    return interval(r1,r2);
} // expm1(...) 

// Die folgenden Konstanten q_expx2_p, q_expx2_m beziehen sich 
// auf eps(expx2) = 4.618958e-16;    f(x) = e^{+x^2};

const real q_expx2_p = 4503599627370500.0 / 4503599627370496.0;  // (1+e(f))
const real q_expx2_m = 9007199254740985.0 / 9007199254740992.0;  // (1-e(f))

interval expx2(const interval& x)
// e^(+x^2);  Blomquist, 25.07.06;
{
    real y,r1,r2,Sz,Iz;
    interval z = abs(x);
    Sz = Sup(z);  Iz = Inf(z);
    // Berechnung einer Unterschranke:
    y = expx2(Iz);
    r1 = y;
    r1 = r1 * q_expx2_m; // Abrunden
    if (r1<1.0) r1 = 1.0;
    // Berechnung einer Oberschranke:
    r2 = (Sz==Iz)? y * q_expx2_p : expx2(Sz) * q_expx2_p;
    if (Sz==0) r2 = 1.0;

    return interval(r1,r2);
} // expx2 (Intervallfunktion)

// Die folgenden Konstanten q_expx2m1_p, q_expx2m1_m beziehen sich 
// auf eps(expx2m1) = 4.813220e-16;    f(x) = e^{+x^2}-1;

const real q_expx2m1_p = 4503599627370500.0 / 4503599627370496.0; // (1+e(f))
const real q_expx2m1_m = 9007199254740985.0 / 9007199254740992.0; // (1-e(f))
const real expx2m1_0   = comp(0.5,-510);  // 2^(-511)

void sqr2uv(const real&, real&, real&);

real expx2m1_intv(const real& x)
// Zur Implementierung der Intervallfunktion;
// e^(+x^2)-1; rel. Fehlerschranke: eps = 4.813220E-16 = e(f) gilt
// fuer alle x, mit:  2^(-511) <= x <= x0 = 26.64174755704632....
// x0 = 7498985273150791.0 / 281474976710656.0;
// Fuer x > x0 --> Programmabbruch wegen Overflow;
// Fuer 0 <= x < 2^(-511) werden die Funktionswerte auf Null gesetzt.
// Ausfuehrlich getestet;  Blomquist, 10.08.2006;
{
    real t(x),u,v,y,res(0);
    int ex;
    if (t<0) t = -t;  // t >= 0;
    
    if (t>=6.5) res = expx2(t);
    else 
    {
	ex = expo(t);
	sqr2uv(x,u,v);  // u := x*x und v aus S(2,53); 
	if (ex>=2) // g4(x)
	{
	    y = exp(u); 
	    res = 1 - v*y;
	    res = y - res;
	}
	else 
	    if (ex>=-8) res = expm1(u) + v*exp(u); // g3(x)
	    else 
		if (ex>=-25) { // g2(x)
		    y = u*u;
                    times2pown(y,-1);
		    res = (1+u/3)*y + u;
		}
		else 
		    if(ex>=-510) res = u;  // g1(x)
    }

    return res;
} // expx2m1_intv

interval expx2m1(const interval& x)
// e^(+x^2)-1;  Blomquist, 10.08.06;
{
    real y,r1,r2,Sz,Iz;
    interval z = abs(x);
    Sz = Sup(z);  Iz = Inf(z);
    // Berechnung einer Unterschranke:
    y = expx2m1_intv(Iz);
    r1 = y;
    r1 = r1 * q_expx2m1_m; // Abrunden
    // Berechnung einer Oberschranke:
    if (Sz < expx2m1_0) 
    { 
	r2 = MinReal;
	if (Sz==0) r2 = 0.0;
    }
    else r2 = (Sz==Iz)? y * q_expx2m1_p : expx2m1_intv(Sz) * q_expx2m1_p;

    return interval(r1,r2);
} // expx2m1 (Intervallfunktion)





// ------  1-eps and 1+eps for function lnp1, Blomquist 28,07,03;  -----------
// ------------------------  eps = 2.5082e-16  -------------------------------
static real q_lnp1m = 9007199254740986.0 / 9007199254740992.0; // 1-eps
static real q_lnp1p = 4503599627370501.0 / 4503599627370496.0; // 1+eps
// ---------------------------------------------------------------------------

interval lnp1(const interval& x) throw()
// returns an inclusion of ln(1+t), with t in x; Blomquist 28.07.03;
{ 
    real ix=Inf(x), sx=Sup(x),a,b;   // ln(1+x) <= [a,b]
    // Calculating the lower bound a:
    int sgn_ix = sign(ix), ex_ix = expo(ix), sgn_sx,ex_sx;
    if (!sgn_ix) a = 0;  // optimal lower bound!
    else if (sgn_ix>0)
	     a = (ex_ix<=-53) ? pred(ix) : lnp1(ix) * q_lnp1m;
         else 
	     a = (ex_ix<=-54) ? pred(ix) : lnp1(ix) * q_lnp1p;
    // Calculating the upper bound b:
    if (ix == sx) { sgn_sx = sgn_ix; ex_sx = ex_ix; }
    else { sgn_sx = sign(sx); ex_sx = expo(sx); }
    if (sgn_sx>=0)
	b = (ex_sx<=-49) ? sx : lnp1(sx) * q_lnp1p;
    else // sx < 0:
	b = (ex_sx<=-50) ? sx : lnp1(sx) * q_lnp1m;
    return interval(a,b); // ln(1+x) in [a,b]
} // lnp1

interval sin    (const interval &a) throw() { return j_sin(a); }
interval cos    (const interval &a) throw() { return j_cos(a); }
interval tan    (const interval &a) throw() { return j_tan(a); }
interval cot    (const interval &a) throw() { return j_cot(a); }

interval asin   (const interval &a)         { return j_asin(a); }
interval acos   (const interval &a)         { return j_acos(a); }
interval atan   (const interval &a)         { return j_atan(a); }
interval acot   (const interval &a)         { return j_acot(a); }

interval exp    (const interval &a) throw() { return j_exp(a); }
interval ln     (const interval &a)         { return j_log(a); }

interval sinh   (const interval &a) throw() { return j_sinh(a); }
interval cosh   (const interval &a) throw() { return j_cosh(a); }
interval tanh   (const interval &a) throw() { return j_tanh(a); }
interval coth   (const interval &a) throw() { return j_coth(a); }

interval asinh  (const interval &a)         { return j_asnh(a); }
interval acosh  (const interval &a)         { return j_acsh(a); }
interval atanh  (const interval &a)         { return j_atnh(a); }
interval acoth  (const interval &a)         { return j_acth(a); }

/*!
\param a The value for which to compute the value of the error function
\return The computed result of the error function

\sa erf(const real & arg)
*/
interval erf    (const interval &a)         { return j_erf(a);  }
/*!
\param a The value for which to compute the value of the complementary error function
\return The computed result of the complementary error function

\sa erf(const real & arg)
*/
interval erfc   (const interval &a)         { return j_erfc(a); }

//interval pow    (const interval &a, const interval &b) throw(ERROR_INTERVAL_STD_FKT_OUT_OF_DEF)
//{
//	if(Inf(a)>0)
//		return j_exp(b*ln(a));
//	else if(Inf(a)==0 && Inf(b)>=0)
//	{
//		if(Sup(a)>=1)
//			return interval(0,pow(Sup(a),Sup(b)));
//		else
//			return interval(0,pow(Sup(a),Inf(b)));
//	}
//	else
//	{
//		cxscthrow(ERROR_INTERVAL_STD_FKT_OUT_OF_DEF("interval pow(const interval &,const interval &)"));
//		return interval(0,0);
//	}
//}

inline a_intv _a_intv(const interval &x)
{       
       return *((const a_intv *)(&x));
}
inline interval _interval(const a_intv &x)
{       
       return *((const interval *)(&x));
}

interval pow    (const interval &a, const interval &b) throw()
       { 
         interval res; 
         if(Inf(a)==0 && Inf(b)>=0)
         {
           if(Sup(a)>=1)
             res=interval(0,pow(Sup(a),Sup(b)));
           else
             res=interval(0,pow(Sup(a),Inf(b)));
         }
	 else res = _interval(i_pow(_a_intv(a),_a_intv(b))); 

         if (Inf(res) <= Sup(res))
           return res;
         else
           return interval(Sup(res),Inf(res));
       }
       
interval power(const interval& a, int n)
// Calculating a^n; 
// Examples: [-1,4]^2 = [0,16];  [-1,4]^3 = [-16,+64];   
{
    bool neg(n<0);
    int N(n),k(-1),r; 
    interval res,h;
    if (neg) N = -N;
    if (N==0) res = 1;
    else if (N==1) res = a;
    else // N > 1:
	if (Inf(a)>=MinReal) res = exp(N*ln(a));
	else if (Sup(a)<=-MinReal) {
	    h = interval(-Sup(a),-Inf(a));
	    res = exp(N*ln(h));
	    if (N%2 == 1) res = -res;
	} 
	else 
	{
	    h = abs(a);
	    while(N>0) 
	    {
		k++;
		r = N % 2;
		if (k==0) 
		    if (r==1) res=a; else res=1;
		else {
		    h = h*h;
		    if (r==1) res *= h;
		}
		N = N / 2;
	    }
	}
    if (neg) res = 1/res;
    return res;
} // Power

/*
interval Power(const interval& a, int n)
// Calculating a^n; 
// Examples: [-1,4]^2 = [0,16];  [-1,4]^3 = [-16,+64];   
{
    bool neg(n<0);
    int N(n),k(-1),r; 
    interval res,h;
    if (neg) N = -N;
    if (N==0) res = 1;
    else if (N==1) res = a;
    else // N > 1:
    {   
	h = abs(a);
	while(N>0) 
	{
	    k++;
	    r = N % 2;
	    if (k==0) 
		if (r==1) res=a; else res=1;
	    else {
		h = h*h;
		if (r==1) res *= h;
	    }
	    N = N / 2;
	}
    }
    if (neg) res = 1/res;
    return res;
} // Power
*/

//----------------------------------------------------------------------------
// Purpose: The local function 'Power()' is used to compute a lower or an
//    upper bound for the power function with real argument and integer
//    exponent, respectively.
// Parameters:
//    In: 'x'      : real argument
//        'n'      : integer exponent
//        'RndMode': rounding mode,
//                   (-1 = downwardly directed, +1 = upwardly directed)
// Description:
//    This function is used to speed up the interval power function defined
//    below. The exponentiation is reduced to multiplications using the
//    binary shift method. Depending on 'n', this function is up to 40 times
//    as fast as the standard power function for real argument and real
//    exponent. However, its accuracy is less than one ulp (unit in the last
//    place of the mantissa) since about log2(n) multiplications are executed
//    during computation. Since directed roundings are antisymmetric, one
//    gets
//
//       down(x^n) = -up((-x)^n)   and   up(x^n) = -down((-x)^n)
//
//    for x < 0 and odd n, where 'down' and 'up' denote the downwardly and
//    upwardly directed roundings, respectively.
//----------------------------------------------------------------------------
static real Power (const real & x, int n, int RndMode )
{                         // Signals change of the rounding mode
  int  ChangeRndMode;     // for x < 0 and odd n
  real p, z;

  ChangeRndMode = ( (x < 0.0) && (n % 2 == 1) );

  if (ChangeRndMode) {
    z = -x;
    RndMode = -RndMode;
  }
  else
    z = x;

  p = 1.0;
  switch (RndMode) {                             // Separate while-loops used
    case -1 : while (n > 0) {                    // to gain speed at runtime
                if (n % 2 == 1) p = muld(p,z);   //--------------------------
                n = n / 2;
                if (n > 0) z = muld(z,z);
              }
              break;
    case +1 : while (n > 0) {
                if (n % 2 == 1) p = mulu(p,z);
                n = n / 2;
                if (n > 0) z = mulu(z,z);
              }
              break;
  }

  if (ChangeRndMode)
    return -p;
  else
    return p;
}

//----------------------------------------------------------------------------
// Purpose: This version of the function 'Power()' is used to compute an
//    enclosure for the power function with interval argument and integer
//    exponent.
// Parameters:
//    In: 'x': interval argument
//        'n': integer exponent
// Description:
//    In general, this implementation does not deliver a result of maximum
//    accuracy, but it is about 30-40 times faster than the standard power
//    function for interval arguments and interval exponents. The resulting
//    interval has a width of approximately 2*log2(n) ulps. Since x^n is
//    considered as a monomial, we define x^0 := 1. For negative exponents
//    and 0 in 'x', the division at the end of the function will cause a
//    runtime error (division by zero).
//----------------------------------------------------------------------------
interval Power (const interval & x, int n )
{
  int  m;
  real Lower, Upper;

  if (n == 0) return(_interval(1.0,1.0));

  if (n > 0)  m = n;  else  m = -n;

  if ( (0.0 < Inf(x)) || (m % 2 == 1) ) {
    Lower = Power(Inf(x),m,-1);
    Upper = Power(Sup(x),m,+1);
  }
  else if (0.0 > Sup(x)) {
    Lower = Power(Sup(x),m,-1);
    Upper = Power(Inf(x),m,+1);
  }
  else {
    Lower = 0.0;
    Upper = Power(AbsMax(x),m,+1);
  }

  if (n > 0)
    return(_interval(Lower,Upper));
  else                                     // Causes a runtime error
    return(1.0/_interval(Lower,Upper));    // if 0 in 'x'.
}

interval Pi ( )                                    // Enclosure of constant pi
  { return(4.0*atan(_interval(1.0,1.0))); }        //-------------------------

// Error bounds for the interval function ln_sqrtx2y2:
real ln_x2y2_abs(2.225076E-308); // Absolute error bond
real q_lnx2y2p(4503599627370502.0 / 4503599627370496.0); // 1+e(f)
real q_lnx2y2m(9007199254740984.0 / 9007199254740992.0); // 1-e(f)
// With the following b0 
// real b0 = 6369051672525773.0 / 30191699398572330817932436647906151127335369763331523427009650401964993299137190816689013801421270140331747000246110759198164677039398341060491474011461568349195162615808.0;
real b0 = MakeHexReal(0,1022-510,0x0016A09E,0x667F3BCD);
// it holds:
// 1. b < bo  ==> g(b) := (0.5*b)*b < MinReal with rounding downwards
// 2. b >= b0 ==> g(b) := (0.5*b)*b >= MinReal with arbitrary rounding
//                                   modus by the two multiplications.

interval ln_sqrtx2y2(const interval& x, const interval& y) throw()
// ln( sqrt(x^2+y^2) ) == 0.5*ln(x^2+y^2);   Blomquist, 22.11.03;
{
    interval ax=abs(x), ay=abs(y);
    real Ix=Inf(ax), Sx=Sup(ax), Iy=Inf(ay), Sy=Sup(ay),f,u1,u2;
    // Calculating the lower bound u1:
    f = ln_sqrtx2y2(Ix,Iy);
    if (Ix==1 && Iy<b0 || Iy==1 && Ix<b0) {
    // f in the denormalized range!
	u1 = f - ln_x2y2_abs;   // directed rounding not necessary!
	if (sign(u1)<0) u1 = 0;
    } else  u1 = (sign(f)<0) ? f*q_lnx2y2p : f*q_lnx2y2m;
    // Calculating the upper bound u2:
    if (Ix==Sx && Iy==Sy) // x and y are point-intervals
	if (Sx==1 && Sy<b0 || Sy==1 && Sx<b0) {
        // f in the denormalized range!
	    u2 = (Sy==0 || Sx==0) ? f : f+ln_x2y2_abs; 
	} else  u2 = (sign(f)<0) ? f*q_lnx2y2m : f*q_lnx2y2p;
    else // x or y is no point-interval:
    {
	f = ln_sqrtx2y2(Sx,Sy);
	if (Sx==1 && Sy<b0 || Sy==1 && Sx<b0)
        // f in the denormalized range!
	    u2 = (sign(Sy)==0 || sign(Sx)==0) ? f : f+ln_x2y2_abs; 
	else  u2 = (sign(f)<0) ? f*q_lnx2y2m : f*q_lnx2y2p;
    }
    return interval(u1,u2);
} // ln_sqrtx2y2


// Constants for the interval function acoshp1(x) = acosh(1+x):
// (1+e(f)):
static const real q_acoshp1p(4503599627370503.0/4503599627370496.0);  
// (1-e(f))
static const real q_acoshp1m(9007199254740981.0/9007199254740992.0);  

interval acoshp1(const interval& x)
// acoshp1;  Blomquist, 28.03.2005;
{
    real r1,r2,sx,ix;
    sx = Sup(x);  ix = Inf(x);
    // Calculating of the lower bound r1:
    r2 = acoshp1(ix);
    r1 = r2 * q_acoshp1m;
    // Calculating of the upper bound r2:
    r2 = (sx==ix)? r2*q_acoshp1p : acoshp1(sx)*q_acoshp1p;
    return interval(r1,r2);

} // acoshp1 (interval arguments)


} // namespace cxsc

