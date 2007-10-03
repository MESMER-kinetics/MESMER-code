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

/* CVS $Id: rmath.cpp,v 1.14 2007/05/16 14:59:03 cxsc Exp $ */

#include "rmath.hpp"

#include "rtsrmath.h"

namespace cxsc {
// Constants and values for the function ln_sqrtx2y2:
real ln2_1067(6505485212531678.0/8796093022208.0); // 1067*ln(2)
// Exponents of the interval bounds:
int B_lnx2y2_1[22] = {21, 31, 51, 101, 151, 201, 251, 301, 351, 401,
                      451, 501, 551, 601, 651, 701, 751, 801, 851, 
                      901, 951, 1025}; 
int B_lnx2y2_2[22] = {-1021,-949,-899,-849,-799,-749,-699,-649,-599,
                      -549,-499,-449,-399,-349,-299,-249,-199,-149,
                      -99,-49,-29,-19}; 
// Optimal values N for N*ln(2):
int B_lnx2y2_N1[21] ={20, 40, 61, 122, 160, 229, 259, 320, 366, 427,
                      488, 518, 549, 610, 671, 732, 763, 825, 885, 
                      945, 976};
real B_lnx2y2_c1[21] = 
{
    7804143460206699.0 / 562949953421312.0, // N*ln(2) with N = 20;
    7804143460206699.0 / 281474976710656.0, // N*ln(2) with N = 40;
    5950659388407608.0 / 140737488355328.0, // N*ln(2) with N = 61;
    5950659388407608.0 / 70368744177664.0,  // N*ln(2) with N = 122;
    7804143460206699.0 / 70368744177664.0,  // N*ln(2) with N = 160;
    5584840163710419.0 / 35184372088832.0,  // N*ln(2) with N = 229;
    6316478613104797.0 / 35184372088832.0,  // N*ln(2) with N = 259;
    7804143460206699.0 / 35184372088832.0,  // N*ln(2) with N = 320;
    8925989082611412.0 / 35184372088832.0,  // N*ln(2) with N = 366;
    5206826964856657.0 / 17592186044416.0,  // N*ln(2) with N = 427;
    5950659388407608.0 / 17592186044416.0,  // N*ln(2) with N = 488;
    6316478613104797.0 / 17592186044416.0,  // N*ln(2) with N = 518;
    6694491811958559.0 / 17592186044416.0,  // N*ln(2) with N = 549;
    7438324235509510.0 / 17592186044416.0,  // N*ln(2) with N = 610;
    8182156659060461.0 / 17592186044416.0,  // N*ln(2) with N = 671;
    8925989082611412.0 / 17592186044416.0,  // N*ln(2) with N = 732;
    4652001140732587.0 / 8796093022208.0,   // N*ln(2) with N = 763;
    5030014339586349.0 / 8796093022208.0,   // N*ln(2) with N = 825;
    5395833564283538.0 / 8796093022208.0,   // N*ln(2) with N = 885;
    5761652788980727.0 / 8796093022208.0,   // N*ln(2) with N = 945;
    5950659388407608.0 / 8796093022208.0    // N*ln(2) with N = 976;
};

int uint_trail(const unsigned int& n)
//  Liefert die Anzahl aufeinanderfolgender Null-Bits am Ende von n;
//  n=0 --> 32;   n=1 --> 0;   n=2 --> 1;   n=4 --> 2;   n=1024 --> 10;
{
    int m1,p=0;
    unsigned int m=n;
    if (m==0) p = 32;
    else { // m != 0, d.h. die folgende while-Schleife bricht ab!
       do
       {
	   m1 = m & 1;  // m1=0(1), wenn letztes Bit von m gleich 0(1) ist.
	   if (m1==0) p++;
	   m = m >> 1;  // Bit-Muster um 1 nach rechts schieben.
       } while(m1==0);
    }
    return p;
} // uint_trail

void sqr2uv(const real& x, real& u, real& v)
// Liefert u,v für: x2 = u + v; EXAKTE Darstellung, falls kein overflow 
// auftritt und v im normalisierten Bereich liegt. u > |v|
// Vorsicht: Funktioniert zunächst nur auf INTEL-Systemen!!!
{
    real a,b,t,y1,y2;
    a = Cut26(x);
    b = x-a;  // x = a+b;
    u = x*x;
    t = u - a*a; // exakte Auswertung!
    y2 = a*b;
    times2pown(y2,1); // y2 = 2*a*b, exakt!
    t -= y2;
// Jetzt fehlt noch:  t2 - b*b, aber b*b wird nicht immer korrekt berechnet,
// daher nochmalige Aufspaltung von b in y1+y2!!
    y1 = Cut25(b);
    y2 = b - y1;   // b = y1+y2, exakt;
    t -= y1*y1;
    if (sign(y2)!=0)
    {
	a = y1*y2;
	times2pown(a,1); // a = 2*y1*y2, exakt!
	t -= a;
	t -= y2*y2;
    }
    v = -t;
} // sqr2uv


//------------------------------------------------------------------

// real-Konstante  expmx2_x0  fuer die Funktion e^(-x2);
const real expmx2_x0 = 7491658466053896.0 / 281474976710656.0;
// Fuer x > expmx2_x0 werden die Funktionswerte auf Null gesetzt.
// Die relative Fehlerschranke e(f) := 4.618919E-16 gilt fuer
// alle |x| <= expmx2_x0 = 26.61571750925.... 

real expmx2(const real& x) throw()
// e^(-x^2);  rel. Fehlerschranke:  eps = 4.618919E-16 = e(f) gilt
// fuer alle |x| <= expmx2_x0 = 26.61571750925....
// Fuer |x| > expmx2_x0 --> expmx2(x) = 0;
// Blomquist, 05.07.04;
{ 
    real t=x,u,v,res=0;
    int ex;
    if (t<0) t = -t;  // t >= 0;
    ex = expo(t);
    if (ex<=-26) res = 1;  // t < 2^(-26)
    else if (ex<=-6)       // t < 2^(-6)
    {
	u = t*t;  v = u; 
	times2pown(v,-1);  // v: 0.5*x2
	res = 1-u*( 1-v*(1-u/3) );
    } else if (t<=expmx2_x0) {
	sqr2uv(x,u,v);  // u:= x*x,v aus S(2,53); x2 = u+v (exakt!)
	res = exp(-u); 
	if (v!=0) { 
	    times2pown(res,500);  // Die Skalierung verhindert, dass 
	    v *= res; // v*exp(-u) in den denormalisierten Bereich faellt
	    res -= v;
	    times2pown(res,-500); // Rueckskalierung
	}
    }
    return res;
} // expmx2

real expx2(const real& x)
// e^(+x^2);  rel. Fehlerschranke:  eps = 4.618958E-16 = e(f) gilt
// fuer alle |x| <= x0 = 26.64174755704632....
// x0 = 7498985273150791.0 / 281474976710656.0;
// Fuer |x| > x0 --> Programmabbruch
// Ausfuehrlich getestet;  Blomquist, 26.07.06;
{
    real t=x,u,v,res;
    int ex;
    if (t<0) t = -t;  // t >= 0;
    ex = expo(t);
    if (ex<=-26) res = 1;  // t < 2^(-26)
    else 
	if (ex<=-6)       // t < 2^(-6)
	{
	    u = t*t;  v = u;
	    times2pown(v,-1);  // v: 0.5*x^2
	    res = 1+u*( 1+v*(1+u/3) );
	} 
	else 
	{
	    sqr2uv(x,u,v);  // u := x*x und v aus S(2,53); 
                            // x^2 = u+v ist exakt!
	    res = exp(u);
	    v *= res;       // v*exp(+u)
	    res += v;
	}
    return res;
} // expx2

real expx2m1(const real& x)
// e^(+x^2)-1; rel. Fehlerschranke: eps = 4.813220E-16 = e(f) gilt
// fuer alle x, mit:  2^(-511) <= x <= x0 = 26.64174755704632....
// x0 = 7498985273150791.0 / 281474976710656.0;
// Fuer x > x0 --> Programmabbruch wegen Overflow;
// Fuer 0 < x < 2^(-511) --> Programmabbruch, denorm. Bereich!!
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
		    else if (ex>=-1073) 
		    {
			std::cerr << "expx2m1: denormalized range!" 
                             << std::endl; exit(1);
		    }
    }

    return res;
} // expx2m1



//------------------------------------------------------------------

real sqrt1px2(const real& x) throw()
// sqrt(1+x^2); Blomquist 13.12.02;
{
    if (expo(x) > 33) return abs(x);
    else return sqrt(1+x*x);
}

// Folgende Konstante sqrtp1m1_s wird gebraucht in 
// real sqrtp1m1(const real& x) throw();  Blomquist, 05.08.03
const real sqrtp1m1_s = 9007199254740984.0 / 1125899906842624.0;

real sqrtp1m1(const real& x) throw()
// sqrtp1m1(x) = sqrt(x+1)-1;
// Blomquist, 05.08.03; 
{
    real y = x;
    int ex = expo(x);
    if (ex<=-50) times2pown(y,-1);  // |x|<2^(-50); fast division by 2 
    else if (ex>=105) y = sqrt(x);  // x >= 2^(+104) = 2.02824...e+31
    else if (ex>=53) y = sqrt(x)-1; // x >= 2^(+52)  = 4.50359...e+15
    else if (x>-0.5234375 && x<=sqrtp1m1_s ) y = x / (sqrt(x+1) + 1);
    else y = sqrt(x+1)-1;
    return y;
}


//------------------------------------------------------------------

real sqrtx2y2(const real& x, const real& y) throw()
// calculating sqrt(x^2 + y^2) in high accuracy. Blomquist 01.12.02
{
    real a,b,r;
    dotprecision dot;
    int exa,exb,ex;
    a = x;   b = y;
    exa = expo(a);  exb = expo(b);
    if (exb > exa)
    { // Permutation of a,b:
	r = a;  a = b;  b = r;
	ex = exa;  exa = exb;  exb = ex;
    }  // |a| >= |b|
    // a = |a| >= |b| > 0:
    ex = 511 - exa;   // scaling a with 2^ex --> expo(a) == 511, --> a*a and
    times2pown(a,ex); // b*b without overflow. An underflow of b*b will not
    times2pown(b,ex); // affect the accuracy of the result!
    dot = 0;
    accumulate(dot,a,a);
    accumulate(dot,b,b);
    r = rnd(dot);
    r = sqrt(r); // sqrt(...) declared in rmath.hpp
    times2pown(r,-ex); // back-scaling
    return r;
} // sqrtx2y2

//------------------------------------------------------------------

real sqrtx2m1(const real& x) 
// sqrt(x^2-1);  rel. Fehlerschranke:  eps = 2.221305E-16 = e(f)
// Blomquist, 13.04.04;
{    const real c1 = 1.000732421875, 
	       c2 = 44000.0,
               c3 = 1024.0;  // c1,c2,c3 werden exakt gespeichert!
    real res,ep,ep2,s1,s2,x1,x2,arg=x;
    if (sign(arg)<0) arg = -arg; // arg = |x| >= 0
    if (arg <= c1) { // x = 1+ep; x^2-1 = 2*ep + ep^2;
	ep = x - 1; // Differenz rundungsfehlerfrei!
	ep2 = ep*ep;  // ep2 i.a. fehlerbehaftet!
	times2pown(ep,1);  // ep = 2*ep;
	res = sqrt(ep+ep2); // res=y0: Startwert;
        // x - y0^2 = (2*eps - s1^2) + [eps^2 - s2*(y0 + s1)]
	s1 = Cut26(res);  s2 = res - s1; // Startwert y0 = s1 + s2;
	arg = ep - s1*s1;  // arg = 2*eps - s1^2;
	arg += (ep2 - s2*(res+s1));  // arg = x - y0^2
	if (sign(arg)>0) {
	arg = arg / res;
	times2pown(arg,-1);
	res += arg;  // 1. Newton-Schritt beendet; eps = 2.221261E-16
	}
    } else if (arg<c2) { 
        // x-y0^2 = [(x1^2-1)-s1^2] + [x2*(x+x1)-s2*(y0+s1)]
	x1 = Cut26(arg);  x2 = arg - x1;  // arg = x = x1 + x2;
	ep2 = x2*(arg+x1);  // ep2 ist fehlerbehaftet
	x2 = x1*x1;  ep = x2-1;
	res = sqrt(ep+ep2); // res ist Startwert für Newton-Verfahren
	s1 = Cut26(res);  s2 = res - s1;  // Startwert = s1 + s2;
	ep2 = ep2 - s2 * (res+s1); // ep2 = [x2*(x+x1)-s2*(y0+s1)]
	if (arg<c3) ep -= s1*s1;   // ep = (x1^2-1) - s1^2;
	else {
	    x2 -= s1*s1;  // x2 = x1^2-s1^2
	    ep = x2 - 1; }         // ep = (x1^2-s1^2) - 1;
	ep += ep2;        // ep = x - y0^2;
	ep /= res;
	times2pown(ep,-1);
	res = res + ep;  // 1. Newton-Schritt in hoher Genauigkeit
                         // beendet;  eps = 2.221305E-16
    } else { // arg = |x| >= 44000;
	res = -1/arg;
	times2pown(res,-1); // Multiplikation mit 0.5;
	res += arg;  // res = x - 1/(2*x);  eps = 2.221114E-16
    }
    return res;
} // sqrtx2m1 (Punktargumente)

//------------------------------------------------------------------

real sqrt1mx2(const real& x) throw(STD_FKT_OUT_OF_DEF)
// sqrt(1-x2);  rel. Fehlerschranke:  eps = 3.700747E-16 = e(f)
// Blomquist, 19.06.04;
{ 
    real t=x,res;
    int ex;
    if (sign(t)<0) t = -t; // Argument t >=0;
    if (t>1) cxscthrow(STD_FKT_OUT_OF_DEF("real sqrt1mx2(const real&)"));
    // For argument t now it holds: 0 <= t <=1;
    ex = expo(t);
    if (ex<=-26) res = 1; // t < 2^(-26) --> res = 1
    else if (ex<=-15) {   // t < 2^(-15) --> res = 1-x2/2
	res = t*t;  
	times2pown(res,-1);
	res = 1 - res;
    } else {
	if (ex>=0) 
	{  // ex>=0 --> t>=0.5; 
	    res = 1-t;   // res: delta = 1-t;
	    t = res * res;
	    times2pown(res,1); // res: 2*delta
	    res = res - t;     // res: 1-x2 = 2*delta - delta2
	} else res = 1-t*t;    // res: Maschinenwert von 1-x2
	res = sqrt(res);       // res: Nullte Naeherung 
    }
    return res;
} // sqrt1mx2


//------------------------------------------------------------------

int Interval_Nr(int* v, const int& n, const int& ex)
// n>0 subintervals:   |...\|...\|...\ ..... |...|
// subinterval Nr.:      0    1    2   .....  n-1  
{
    int i=0,j=n,k;  // n>0:  Number of subintervals
    do {
	k = (i+j)/2;
	if (ex < v[k]) j = k-1;
	else i = k+1;
    } while(i<=j);
    return j;  // x with ex=expo(x) lies in the subinterval number j  
}

//------------------------------------------------------------------

real ln_sqrtx2y2(const real& x, const real& y) 
                                throw(STD_FKT_OUT_OF_DEF)
// ln( sqrt(x^2+y^2) ) == 0.5*ln(x^2+y^2); Blomquist, 21.11.03;
// Relative error bound: 5.160563E-016;
// Absolute error bound: 2.225075E-308; if x=1 and 0<=y<=b0;
{
    int j,N;
    real a,b,r,r1;
    dotprecision dot;
   
    a = sign(x)<0 ? -x : x;  // a = |x| >= 0;
    b = sign(y)<0 ? -y : y;  // b = |y| >= 0;
    int exa=expo(a), exb=expo(b), ex;
    if (b > a)
    {
	r = a;   a = b;   b = r;
	ex = exa;   exa = exb;   exb = ex;
    }
    // It holds now:   0 <= b <= a 
    if (sign(a)==0)
        cxscthrow(STD_FKT_OUT_OF_DEF
                    ("real ln_sqrtx2y2(const real&, const real&)"));
    if (exa>20) // to avoid overflow by calculating a^2 + b^2
    {  // a>=2^(20):
	j = Interval_Nr(B_lnx2y2_1,21,exa); // j: No. of subinterval
	N = B_lnx2y2_N1[j];    // N: Optimal int value
	if (exb-exa > -25)
	{   // For (exb-exa>-25) we use the complete term:  
            // N*ln(2) + [ln(2^(-N)*a)+0.5*ln(1+(b/a)^2)]
	    b = b/a;  // a > 0
	    b = lnp1(b*b);
	    times2pown(b,-1);  // exact division by 2
	    times2pown(a,-N);
	    r = b + ln(a); // [ ... ] calculated!
	    r += B_lnx2y2_c1[j];
	}
	else { // For (exb-exa<=-25) only two summands!: 
	    times2pown(a,-N);
	    r = ln(a) + B_lnx2y2_c1[j];
	}
    }
    else  // exa<=20 or a<2^(20):
    {     // Now calculation of a^2+b^2 without overflow:
	if (exa<=-20) // to avoid underflow by calculating a^2+b^2
	    if (exa<=-1022) // a in the denormalized range
	    {
		r = b/a;
		r = lnp1(r*r);  times2pown(r,-1); // r: 0.5*ln(1+..)
		times2pown(a,1067);
		r += ln(a);     // [ .... ] ready
		r -= ln2_1067;  // rel. error = 2.459639e-16;
	    }
	    else // MinReal=2^(-1022) <= a < 2^(-20)
	    {    // Calculating the number j of the subinterval:
		j = 20 - Interval_Nr(B_lnx2y2_2,21,exa);
		r = a;  times2pown(r,B_lnx2y2_N1[j]);
		r = ln(r);  // r: ln(2^N*a);
		if (exb-exa > -25) { // calculating the complete term
		    b = b/a;
		    a = lnp1(b*b);
		    times2pown(a,-1);
		    r += a;  // [ ... ] ready now
		}
                // We now have: exb-exa<=-25,  ==>  b/a <= 2^(-24);
		r -= B_lnx2y2_c1[j]; // 0.5*ln(1+(b/a)^2) neglected!
                // relative error = 4.524090e-16 in both cases;
	    }
	else // calculation of a^2+b^2 without overflow or underflow:
	{   // exa>-20  respective  a>=2^(-20):
	    dot = 0;
	    accumulate(dot,a,a);
	    accumulate(dot,b,b);  // dot = a^2+b^2, exact!
	    real s = rnd(dot);    // s = a^2 + b^2, rounded!
	    if (s>=0.25 && s<=1.75)
		if (s>=0.828125 && s<=1.171875)
		{ // Series:
		    if (a==1 && exb<=-28)
		    {
			r = b;  times2pown(r,-1);
			r *= b;
		    }
		    else {
			dot -= 1;
			r = rnd(dot); // r = a^2+b^2-1 rounded!
			r = lnp1(r);
			times2pown(r,-1);
		    }
		}
		else { // Reading dot = a^2+b^2 twice:
		    r = rnd(dot);
		    dot -= r;
		    r1 = rnd(dot); // a^2+b^2 = r+r1, rounded!
		    r1 = lnp1(r1/r);
		    r = ln(r) + r1;
		    times2pown(r,-1); // exact division by 2
		}
	    else { // calculating straight from: 0.5*ln(x^2+y^2)
		r = ln(s);
		times2pown(r,-1);
	    } 
	}
    }
    return r;
} // ln_sqrtx2y2

typedef union { double f; char intern[7]; } help_real;

real Cut24(const real& x){
    // y = Cut24(x) liefert ein y, das mit den ersten 24 Mantissenbits
    // von x übereinstimmt, das hidden bit ist dabei mitgezählt!
    // Die restlichen 53-24=29 Mantissenbits werden auf Null gesetzt.
  help_real y;
  y.f = _double(x);
#if INTEL
  y.intern[3] &= 224;
  y.intern[0] = y.intern[1] = y.intern[2] = 0;
#else
  y.intern[4] &= 224;
  y.intern[7] = y.intern[6] = y.intern[5] = 0;
#endif
  return real(y.f);
}

real Cut25(const real& x){
    // y = Cut25(x) liefert ein y, das mit den ersten 25 Mantissenbits
    // von x übereinstimmt, das hidden bit ist dabei mitgezählt!
    // Die restlichen 53-25=28 Mantissenbits werden auf Null gesetzt.
  help_real y;
  y.f = _double(x);
#if INTEL
  y.intern[3] &= 240;
  y.intern[0] = y.intern[1] = y.intern[2] = 0;
#else
  y.intern[4] &= 240;
  y.intern[7] = y.intern[6] = y.intern[5] = 0;
#endif
  return real(y.f);
}

real Cut26(const real& x){
    // y = Cut26(x) liefert ein y, das mit den ersten 26 Mantissenbits
    // von x übereinstimmt, das hidden bit ist dabei mitgezählt!
    // Die restlichen 53-26=27 Mantissenbits werden auf Null gesetzt.
  help_real y;
  y.f = _double(x);
#if INTEL
  y.intern[3] &= 248;
  y.intern[0] = y.intern[1] = y.intern[2] = 0;
#else
  y.intern[4] &= 248;
  y.intern[7] = y.intern[6] = y.intern[5] = 0;
#endif
  return real(y.f);
}

static real q_acoshp1[5] = // Polynomial coefficients of Q_4(x) 
           // roundet to nearest. acosh(1+x) = sqrt(2*x)*Q_4(x)
{ 1.0 / 1.0,                                    // q_acoshp1[0]
  -6004799503160661.0 / 72057594037927936.0,
  +5404319552844595.0 / 288230376151711744.0,
  -6433713753386423.0 / 1152921504606846976.0,
  +8756999275442631.0 / 4611686018427387904.0   // q_acoshp1[4]
};

static const real c_ln2_B = 6243314768165359.0 / 9007199254740992.0;
// c_ln2_B < ln(2) is the nearest machine number for ln(2) with an
// absolute error < 2.3190469E-17;

real acoshp1(const real& x) throw()
// acoshp1(x) = acosh(1+x);  rel. error: eps = 7.792706E-16 = e(f)
// Ausfuehrlich getestet;  Blomquist, 27.03.05;
{ 
    real res;
    int ex(expo(x));
    if (x<0) 
	cxscthrow(STD_FKT_OUT_OF_DEF("real acoshp1(const real&)"));
    // For argument x now it holds: 0 <= x <= MaxReal;
    if (ex<=-50) res = sqrt(2*x); // 0<=x<2^(-50): acoshp1(x)=sqrt(2x)
    else if (ex<=-9) // 2^(-50)<=x<2^{-9}: acoshp1(x)=sqrt(2x)*Q_4(x)
      res = sqrt(2*x)*((((q_acoshp1[4]*x+q_acoshp1[3])*x+q_acoshp1[2])
			    *x+q_acoshp1[1])*x + q_acoshp1[0]);
    else if (ex<=0) res = lnp1(x+sqrt(2*x+x*x));     // range A_3
    else if (ex<=50) res = lnp1(x*(1+sqrt(1+2/x)));  // range A_4
    else if (ex<=1022) res = ln(2*x);                // range A_5
    else res = ln(x) + c_ln2_B;                      // range A_6

    return res;
} // acoshp1

extern "C" {
   void r_lfsr(void) {;} // Siehe real.hpp in real_ari...?!?!
}

} // namespace cxsc

