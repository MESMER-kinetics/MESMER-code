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

/* CVS $Id: interval.inl,v 1.15 2007/05/16 14:58:58 cxsc Exp $ */

namespace cxsc {

// ---- Konstruktoren ----

inline interval::interval(const real &a,const real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
                                                       : inf(a), sup(b)
{
   if (a > b) 
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval::interval(const real &a,const real &b)"));
}

/*inline interval::interval(int a,int b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
                                      : inf(a), sup(b)
{
   if (a > b) 
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval::interval(int a,int b)"));
}

inline interval::interval(const double & a,const double & b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
                                                             : inf(a), sup(b)
{
   if (a > b) 
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval::interval(const double & a,const double & b)"));
}
*/

inline interval& interval::operator= (const real& a)
{
   inf = sup = a;
   return *this;
}


// ---- Typwandlungen ----

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::interval::interval(const real&, const real&)
*/
inline interval _unchecked_interval(const real& a, const real& b) 
{
  interval tmp;
  tmp.inf = a;
  tmp.sup = b;
  return tmp;
}


// ---- Standardfunkt ---- (arithmetische Operatoren)

inline interval operator-(const interval &a) throw() { return interval (-a.sup, -a.inf); }
inline interval operator+(const interval &a) throw() { return a; }

inline interval & operator +=(interval &a,const interval &b) throw() { return a=a+b; }
inline interval & operator +=(interval &a,const real &b)     throw() { return a=a+b; }
inline interval & operator -=(interval &a,const interval &b) throw() { return a=a-b; }
inline interval & operator -=(interval &a,const real &b)     throw() { return a=a-b; }
inline interval & operator *=(interval &a,const interval &b) throw() { return a=a*b; }
inline interval & operator *=(interval &a,const real &b)     throw() { return a=a*b; }
inline interval & operator /=(interval &a,const interval &b) throw() { return a=a/b; }
inline interval & operator /=(interval &a,const real &b)     throw() { return a=a/b; }


inline interval operator |(const interval &a,const interval &b) throw() 
{
   return interval((a.inf<b.inf)?a.inf:b.inf,(a.sup>b.sup)?a.sup:b.sup);
}
inline interval operator &(const interval &a,const interval &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   return interval((a.inf>b.inf)?a.inf:b.inf,(a.sup<b.sup)?a.sup:b.sup);
}
inline interval operator |(const real &a,const interval &b) throw() 
{
   return interval((a<b.inf)?a:b.inf,(a>b.sup)?a:b.sup);
}
inline interval operator &(const real &a,const interval &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   return interval((a>b.inf)?a:b.inf,(a<b.sup)?a:b.sup);
}
inline interval operator |(const interval &a,const real &b) throw() 
{
   return interval((a.inf<b)?a.inf:b,(a.sup>b)?a.sup:b);
}
inline interval operator |(const real &a,const real &b) throw()
{
   if(a>b) return interval(b,a);
   else    return interval(a,b);
}
inline interval operator &(const interval &a,const real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   return interval((a.inf>b)?a.inf:b,(a.sup<b)?a.sup:b);
}
inline interval & operator |=(interval &a,const interval &b) throw() 
{
   a.inf=(a.inf<b.inf)?a.inf:b.inf,a.sup=(a.sup>b.sup)?a.sup:b.sup;
   return a;
}
inline interval & operator &=(interval &a,const interval &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   a.inf=(a.inf>b.inf)?a.inf:b.inf,a.sup=(a.sup<b.sup)?a.sup:b.sup;
   if(a.inf>a.sup)
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval & operator &=(interval &a,const interval &b)"));
   return a;
}
inline interval & operator |=(interval &a,const real &b) throw() 
{
   a.inf=(a.inf<b)?a.inf:b,a.sup=(a.sup>b)?a.sup:b;
   return a;
}
inline interval & operator &=(interval &a,const real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   a.inf=(a.inf>b)?a.inf:b,a.sup=(a.sup<b)?a.sup:b;
   if(a.inf>a.sup)
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval & operator &=(interval &a,const real &b)"));
   return a;
}

// --- Vergleichsoperationen ----
inline bool operator ==(const interval &a,const interval &b) throw() {   return(a.inf==b.inf && a.sup==b.sup); }
inline bool operator !=(const interval &a,const interval &b) throw() {   return(a.inf!=b.inf || a.sup!=b.sup); }
inline bool operator ==(const real &r,const interval &a)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const real &r,const interval &a)     throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const real &r)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const real &r)     throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const int &r,const interval &a)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const int &r,const interval &a)     throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const int &r)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const int &r)     throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const long &r,const interval &a)    throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const long &r,const interval &a)    throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const long &r)    throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const long &r)    throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const double &r,const interval &a)  throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const double &r,const interval &a)  throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const double &r)  throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const double &r)  throw() {   return(r!=a.inf || r!=a.sup); }

// --- Mengenvergleiche ---
// <,>,...
inline bool operator <=(const interval &a,const interval &b) throw()
{
   return(a.inf>=b.inf && a.sup<=b.sup);   
}
inline bool operator >=(const interval &a,const interval &b) throw()
{
   return(a.inf<=b.inf && a.sup>=b.sup);   
}
inline bool operator <(const interval &a,const interval &b) throw()
{
   return(a.inf>b.inf && a.sup<b.sup);   
}
inline bool operator >(const interval &a,const interval &b) throw()
{
   return(a.inf<b.inf && a.sup>b.sup);   
}

inline bool operator <=(const real &a,const interval &b) throw()
{
   return(a>=b.inf && a<=b.sup);   
}
inline bool operator >=(const real &a,const interval &b) throw()
{
   return(a<=b.inf && a>=b.sup);   
}
inline bool operator <(const real &a,const interval &b) throw()
{
   return(a>b.inf && a<b.sup);   
}

inline bool operator <=(const interval &a,const real &b) throw()
{
   return(a.inf>=b && a.sup<=b);   
}
inline bool operator >=(const interval &a,const real &b) throw()
{
   return(a.inf<=b && a.sup>=b);   
}
inline bool operator >(const interval &a,const real &b) throw()
{
   return(a.inf<b && a.sup>b);   
}

inline bool operator !(const interval &a) throw() { return (a.inf <= 0.0 && a.sup >= 0.0); }  

inline bool IsEmpty(const interval &a) throw() { return (a.inf>a.sup); }

inline interval abs(const interval &a) throw()
{
   real h1  = abs(a.inf);
   real h2  = abs(a.sup);

   if (IsEmpty(a)) return a;
   if (!a)         
      return interval(0.0, (h1 > h2) ? h1 : h2);
   if (h1 > h2)    
      return interval(h2, h1);

   return interval(h1, h2); 
}

inline real diam(const interval & a) throw()
{
   return subup(a.sup,a.inf); 
}

inline void times2pown(interval& x, const int& n) throw()
{
    real r1,r2;
    int j;
    r1 = x.inf;  r2 = x.sup;
    j = expo(r1) + n;
    if (j >= -1021) r1 = comp(mant(r1),j);
    else 
    {
	j += 1021;
	r1 = comp(mant(r1), -1021);
	if (j<-53)
	{
	    if (sign(r1)>=0) r1 = 0;
	    else r1 = -minreal;
	} else 
	    r1 = muld(r1,comp(0.5,j+1));
    }
    j = expo(r2) + n;
    if (j >= -1021) r2 = comp(mant(r2),j);
    else 
    {
	j += 1021;
	r2 = comp(mant(r2), -1021);
	if (j<-53)
	{
	    if (sign(r2)>0) r2 = minreal;
	    else r2 = 0;
	} else r2 = mulu(r2,comp(0.5,j+1));
    }
    x = _interval(r1,r2);
} // times2pown(...)

} // namespace cxsc

