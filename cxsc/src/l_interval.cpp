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

/* CVS $Id: l_interval.cpp,v 1.5 2007/05/16 14:59:01 cxsc Exp $ */

#include <math.h>
#include "l_interval.hpp"
#include "idot.hpp"
#include "RtsFunc.h"

namespace cxsc {

#define CXSC_Zero 0.

l_interval & l_interval::operator= (const l_interval & a)    throw()
{
   real *newdata=new real[a.prec+1];
   int i;
   for(i=0;i<=a.prec;i++)
      newdata[i]=a.data[i];
   delete [] data;
   data=newdata;
   prec=a.prec;
   return *this;
}   

// ---- Typwandlungen ----

interval::interval(const l_interval & a) throw()
{
   idotakku[0]=0.0;
   a._akku_add(idotakku[0]);
   *this=rnd(idotakku[0]);
}

interval & interval::operator =(const l_interval & a) throw()
{
   idotakku[0]=0.0;
   a._akku_add(idotakku[0]);
   return *this=rnd(idotakku[0]);
}


interval _interval(const l_interval & a) throw()
{
   idotakku[0]=0.0;
   a._akku_add(idotakku[0]);
   return rnd(idotakku[0]);
}

l_interval::l_interval(const dotprecision & a)
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   _allo(stagprec);
   idotakku[0]=a;
   _akku_out();   
}

l_interval::l_interval(const dotprecision & a,const dotprecision & b)
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const dotprecision&,const dotprecision&)"));
   _allo(stagprec);
   UncheckedSetInf(idotakku[0],a);
   UncheckedSetSup(idotakku[0],b);
   _akku_out();   
}

l_interval & l_interval::operator =(const dotprecision & a)
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   if(stagprec!=prec)
   {
      delete [] data;
      _allo(stagprec);
   }
   idotakku[0]=a;
   _akku_out();   
   return *this;
}

l_interval::l_interval(const idotprecision & a)
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   _allo(stagprec);
   idotakku[0]=a;
   _akku_out();   
}

l_interval & l_interval::operator =(const idotprecision & a)
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   if(stagprec!=prec)
   {
      delete [] data;
      _allo(stagprec);
   }
   idotakku[0]=a;
   _akku_out();   
   return *this;
}

l_interval::l_interval(const l_real &a, const l_real &b) 
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(stagprec);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const l_real &a, const l_real &b)"));
   dotakku[0]=a;
   dotakku[1]=b;
   idotakku[0]=idotprecision(dotakku[0],dotakku[1]);
   _akku_out();
}

l_interval::l_interval(const real &a, const l_real &b) 
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(stagprec);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const l_real &a, const l_real &b)"));
   dotakku[0]=a;
   dotakku[1]=b;
   idotakku[0]=idotprecision(dotakku[0],dotakku[1]);
   _akku_out();
}

void l_realz2l_interval(const l_real& lr, const interval& z,
			                        l_interval& li) throw()
{  // converts lr+z to li of type l_interval; Blomquist, 12.10.06;
   // lr+z is included by li.
   // prec(lr) <= prec(li)  must be realized!
    int p = StagPrec(lr);  int q = StagPrec(li);
    if(p>q) 
    { 
	std::cerr << "l_realz2l_interval(const l_real& lr, const interval& z, l_interval& li): incorrect precisions of lr,li !" 
		  << std:: endl; 
	exit (1); 
    }

    for (int i=1; i<=q-1; i++) li[i] = 0; 
    li[q]   = Inf(z);
    li[q+1] = Sup(z);
    if(p<q) for (int i=1; i<=p; i++) li[i] = lr[i];
    else // p=q
    {
	p--;
	for (int i=1; i<=p; i++) li[i] = lr[i];
	li[q]   = addd(lr[p+1],Inf(z)); 
	li[q+1] = addu(lr[p+1],Sup(z));
    }
} // l_realz2l_interval

l_interval::l_interval(const l_real &a, const real &b) 
#if (IndCheck)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(stagprec);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const l_real &a, const l_real &b)"));
   dotakku[0]=a;
   dotakku[1]=b;
   idotakku[0]=idotprecision(dotakku[0],dotakku[1]);
   _akku_out();
}


/*
l_interval _l_interval(const l_real & a, const l_real & b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval _l_interval(const l_real & a, const l_real & b)"));

   dotakku[0]=0.0;
   dotakku[1]=0.0;
   a._akku_add(dotakku[0]);
   b._akku_add(dotakku[1]);
   idotakku[0]=_idotprecision(dotakku[0],dotakku[1]);
   l_interval tmp;
   tmp._akku_out();
   return tmp;
}   
*/
l_interval _unchecked_l_interval(const l_real & lr1, const l_real & lr2) throw()
{
   real inf, sup, tmp;          // fuer Naeherungen, entspricht Interval z
   int i=1;
   l_interval li;
   dotakku[0]=0.0;
   dotakku[1]=0.0;
   lr1._akku_add(dotakku[0]);
   lr2._akku_add(dotakku[1]);
   inf = rnd(dotakku[0],RND_DOWN);
   sup = rnd(dotakku[1],RND_UP);
   while (i<li.prec && !(inf<=0. && sup>=0.)) 
   {
      tmp = inf + (sup-inf)/2.0; // middle(interval)
      li.elem(i) = tmp;
      dotakku[0] -= tmp;
      dotakku[1] -= tmp;
      inf = rnd(dotakku[0],RND_DOWN);
      sup = rnd(dotakku[1],RND_UP);
      i++;
   }
   li.elem(li.prec)=inf;           // Intervall in die letzten Stellen
   li.elem(li.prec+1)=sup;         // schreiben

   return li;
}

idotprecision _idotprecision(const l_interval & a) throw()
{
   return idotprecision(a);
}

idotprecision::idotprecision(const l_interval & a) throw() : inf(0),
                                                             sup(0)
{ 
   a._akku_add(*this);
}

idotprecision & idotprecision::operator =(const l_interval & a) throw()
{
   *this=0;
   a._akku_add(*this);
   return *this;
}

// ---- Standardfunkt ---- (arithmetische Operatoren)

l_interval operator-(const l_interval & a) throw()
{
   int precsave=stagprec;
   stagprec=a.prec;

   l_interval tmp;

   stagprec=precsave;

   int i;
   for(i=0;i<a.prec-1;i++)
      tmp.data[i]=-a.data[i];

   tmp.data[a.prec-1]=-a.data[a.prec];
   tmp.data[a.prec]=-a.data[a.prec-1];

   return tmp;
} 

// LI-LI

l_interval operator+(const l_interval & li1, const l_interval & li2) throw()
{
   l_interval li3;
   idotakku[0]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_add(idotakku[0]);
   li3._akku_out();
   return li3;
}

l_interval operator-(const l_interval & li1, const l_interval & li2) throw()
{
   l_interval li3;
   idotakku[0]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_sub(idotakku[0]);
   li3._akku_out();
   return li3;
} 
 
/* l_interval operator*(const l_interval & li1, const l_interval & li2) throw()
{
   l_interval li3;
   interval stdmul;
   
   stdmul = _interval(li1)*_interval(li2);

   if (abs(Inf(stdmul)) < MinReal) 
      li3 = _l_interval(stdmul) + 0.0; // Was das +0.0 soll ist mir ein Raetsel SR
   else 
   {               // if eingefgt am 30.07.92 von Frederic Toussaint
      idotakku[0]=0.0;
      accumulate(idotakku[0], li1, li2);
      li3._akku_out();
      li3 = li3 & _l_interval(stdmul);  // Nachkorrektur
      // Wie bitte? (SR)
   }
   return li3;
}  */

l_interval operator * (const l_interval& li1, const l_interval& li2) throw()
{  // Blomquist, Neue Version, 21.11.02;
    l_interval li3;
    interval stdmul;
    stdmul = _interval(li1) * _interval(li2); // grobe Einschließung
    idotakku[0]=0.0;
    accumulate(idotakku[0],li1,li2);
    li3._akku_out();       // li3: meist feinere Einschließung! Aber bei zu
             // breiten li1,li2 ist dieses Ergebnis li3 breiter als stdmul!
    li3 = li3 & _l_interval(stdmul);   // Das optimale Ergebnis liegt daher 
              // im Durchschnitt der groben und der feineren Einschließung!
    return li3;
}    

l_interval operator/(const l_interval & li1, const l_interval & li2) throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ 
// geändert am 29.07.92 von Frederic Toussaint
// 13.5.93 AW: neuer Algorithmus nach der Beschreibung von W. Kraemer

   int         i, j, k, Z_sign;
   real        dzmitte;
   interval    dn = _interval(li2),
               dz = _interval(li1),
               stddiv = dz/dn;

   l_interval  li3;

   li3._clear(1);

   if(!li2) 
   {
      cxscthrow(ERROR_LINTERVAL_DIV_BY_ZERO("l_interval operator/(const l_interval & li1, const l_interval & li2)"));
   } else 
   {
      dz = dz/dn;
      if (stagprec > 1) 
      {         // bei == 1 passiert nix!
         idotakku[0] = 0.0;          // X, also Zaehler
         li1._akku_add(idotakku[0]);
         k = 1;
         Z_sign = (Inf(dz) > 0.0);
         do {
            dzmitte = (Inf(dz) + Sup(dz)) / 2.0;  // mid(interval); grob!
            if (!!dz)                           
               if (Sup(abs(dz)) > MinReal)
                  if (diam(dz) < 1e-14 * abs(dzmitte) ) 
                  {
                     li3.elem(k) = dzmitte;
                     // Bestimme Rest in dotakku[2], [3] (Inf, Sup):
                     //
                     dotakku[2] = 0.0;
                     for (i=1; i<=k; i++)                   // Rumpf
                        for (j=1; j<=li2.prec-1; j++)
                           accumulate(dotakku[2], -li2.elem(j), li3.elem(i));

                     dotakku[3]  = dotakku[2];        // Rumpf ist fuer Inf, Sup gleich
                     dotakku[2] += Inf(idotakku[0]);  // Inf += Inf(X)
                     dotakku[3] += Sup(idotakku[0]);  // Sup += Sup(X)
                     if (Z_sign)                      // INCL(Z) > 0.0
                     for (i=1; i<=k; i++) 
                     {
                        accumulate(dotakku[2], -li2.elem(li2.prec+1), li3.elem(i));
                        accumulate(dotakku[3], -li2.elem(li2.prec),   li3.elem(i));
                     } else                             // INCL(Z) < 0.0
                        for (i=1; i<=k; i++) 
                        {
                           accumulate(dotakku[2], -li2.elem(li2.prec),   li3.elem(i));
                           accumulate(dotakku[3], -li2.elem(li2.prec+1), li3.elem(i));
                        }
                     //
                     // Rausrunden in dz
                     dz = _interval(rnd(dotakku[2], RND_DOWN),
                                    rnd(dotakku[3], RND_UP));
                     dz = dz / dn;
                  }
            k++;
         } while (k < stagprec);
      } // if (stagprec > 1)
      li3.elem(stagprec) = Inf(dz);     // Intervall in die letzten Stelle
      li3.elem(stagprec+1) = Sup(dz);
      li3 = li3 & _l_interval(stddiv);  // Nachkorrektur
   } // if keine Null im Nenner
   return li3;
}

// ---- Vergleichsop. ----
bool operator!(const l_interval & li) throw()
{
   idotakku[0]=0.0;
   li._akku_add(idotakku[0]);
   return (!idotakku[0]);
}

/*l_interval::operator void *(void) throw()
{
   idotakku[0]=0.0;
   _akku_add(idotakku[0]);
   return (idotakku[0]);
}*/

bool operator==(const l_interval & li1, const l_interval & li2) throw()
{
   idotakku[0]=0.0;
   idotakku[1]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_add(idotakku[1]);
   return (idotakku[0]==idotakku[1]);
}   

// ---- Mengenvergleiche ----

bool operator<(const l_interval & li1, const l_interval & li2) throw()
{
   idotakku[0]=0.0;
   idotakku[1]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_add(idotakku[1]);
   return (idotakku[0]<idotakku[1]);
} 

bool operator>(const l_interval & li1, const l_interval & li2) throw()
{
   idotakku[0]=0.0;
   idotakku[1]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_add(idotakku[1]);
   return (idotakku[0]>idotakku[1]);
} 

bool operator<=(const l_interval & li1, const l_interval & li2) throw()
{
   idotakku[0]=0.0;
   idotakku[1]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_add(idotakku[1]);
   return (idotakku[0]<=idotakku[1]);
} 

bool operator>=(const l_interval & li1, const l_interval & li2) throw()
{ 
   idotakku[0]=0.0;
   idotakku[1]=0.0;
   li1._akku_add(idotakku[0]);
   li2._akku_add(idotakku[1]);
   return (idotakku[0]>=idotakku[1]);
}

void ConvexHull(const l_interval & li1, const l_interval & li2, l_interval & li3, l_interval & li4) throw()
{
   if(li1<=li2) 
   {                                      // Trivialfall 1
      li3=li2;
      li4=li2;
   } else if(li2<=li1) 
   {                                 // Trivialfall 2
      li3=li1;
      li4=li1;
   } else 
   {                                              // rechne
      idotakku[0]=0.0;
      idotakku[1]=0.0;
      li1._akku_add(idotakku[0]);
      li2._akku_add(idotakku[1]);

      idotakku[0] |= idotakku[1];
      // nun steht das richtige Ergebnis in idotakku[0]
      idotakku[1] = idotakku[0];                              // sichern
      li3._akku_out_inn();                              // Rundung nach innen
      idotakku[0] = idotakku[1];                              // und wiederherstellen
      li4._akku_out();                                  // ... nach aussen
   }
}

void Intersection(const l_interval & li1, const l_interval & li2, l_interval & li3, l_interval & li4) throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
   if(li1<=li2) 
   {                                      // Trivialfall 1
      li3=li1;
      li4=li1;
   } else if(li2<=li1) 
   {                                 // Trivialfall 2
      li3=li2;
      li4=li2;
   } else 
   {                                              // rechne
      idotakku[1]=0.0;
      idotakku[2]=0.0;
      li1._akku_add(idotakku[1]);
      li2._akku_add(idotakku[2]);
      try
      {
         idotakku[0]=(idotakku[1] &= idotakku[2]);
      }
      catch(const EMPTY_INTERVAL &)
      {
         cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("void Intersection(const l_interval & li1, const l_interval & li2, l_interval & li3, l_interval & li4)"));
      }
      // nun steht das richtige Ergebnis in idotakku[0]
      li3._akku_out_inn();                              // Rundung nach innen
      idotakku[0] = idotakku[1];
      li4._akku_out();                                  // ... nach aussen
   }
}

l_real     mid  (const l_interval & li) throw()
{
   l_real lr;

   // --------------------------------------------------------
   //   dotakku[0] = Inf(li) + Sup(li)

   dotakku[0]  = 0.0;
   for (int j=1; j<=li.prec-1; j++)
      dotakku[0] += li.elem(j);

   dotakku[0] += dotakku[0]; // mal 2
   dotakku[0] += li.elem(li.prec);
   dotakku[0] += li.elem(li.prec+1);

   //  ------------------------------------------------------
   //  Division nur bei ungleich 0

   if (dotakku[0] != 0.0) 
   {
      Dotprecision ptr = *dotakku[0].ptr();

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

   lr._akku_out();

   return lr;
}
 
/* void accumulate(idotprecision & d, const l_interval & li1, const l_interval & li2) throw()
{
   // Änderungen am 24.9.92 von F. Toussaint wegen Underflow-Fehlern

   int  i, j, n = 0;
   real tmp;

   for (i=1; i<=li1.prec-1; i++)
      for (j=1; j<=li2.prec-1; j++) 
      {
         tmp = abs(li1.elem(i)*li2.elem(j));
         if (tmp < MinReal && tmp != CXSC_Zero) 
            n++;
         else 
            accumulate(d, _interval(li1.elem(i)), _interval(li2.elem(j)));
      }
      for (i=1; i<=li2.prec-1; i++) 
      {
         tmp = abs(Inf(_interval(li1.elem(li1.prec), li1.elem(li1.prec+1)) * _interval(li2.elem(i))));
         if (tmp < MinReal && tmp != CXSC_Zero) 
            n++;
         else 
            accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                          _interval(li2.elem(i)));
      }
      for (i=1; i<=li1.prec-1; i++) 
      {
         tmp = abs(Inf(_interval(li2.elem(li2.prec), li2.elem(li2.prec+1)) * _interval(li1.elem(i))));
         if (tmp < MinReal && tmp != CXSC_Zero) 
            n++;
         else 
            accumulate(d, _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)),
                          _interval(li1.elem(i)));
      }
      tmp = abs(Inf(_interval(li1.elem(li1.prec), li1.elem(li1.prec+1)) *
                    _interval(li2.elem(li2.prec), li2.elem(li2.prec+1))));
      if (tmp < MinReal && tmp != CXSC_Zero) 
         n++;
      else 
         accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                       _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)));
      if (n > 0) 
         accumulate(d, _interval(_real(n)), 
                       _interval(-MinReal, MinReal));
} */

void accumulate(idotprecision & d, const l_interval & li1, 
                                   const l_interval & li2) throw()
{   // Blomquist, Neue Version vom 21.11.02;
    // Die Änderungen von Toussaint wurden rückgängig gemacht.
   int  i,j;

   for (i=1; i<=li1.prec-1; i++)
      for (j=1; j<=li2.prec-1; j++) 
      {
            accumulate(d, _interval(li1.elem(i)), _interval(li2.elem(j)));
      }
      for (i=1; i<=li2.prec-1; i++) 
      {
            accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                          _interval(li2.elem(i)));
      }
      for (i=1; i<=li1.prec-1; i++) 
      {
            accumulate(d, _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)),
                          _interval(li1.elem(i)));
      }
      accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                    _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)));
}


/* void l_interval::_akku_out() throw()
{  
   // ein im idotakku[0] liegendes Zwischenergebnis
   // wird in der entsprechenden Praezision in das aufrufende l_interval
   // gerundet. Rundung, also Einschluss von aussen!!
   // Alle MinReal wurden geändert in minreal, Blomquist 19.11.02;
   _clear(1);
   interval z;
   real tmp, r, s;
   int i = 1,
   
   n = 0;
   z = rnd(idotakku[0]);

   while (i<prec && !(CXSC_Zero <= z)) 
   {       
      // z ist Intervall, deshalb <=
      if (succ(succ(Inf(z))) < Sup(z))
         break;                            // bei zu grobem Einschluss:
                                           // Intervall direkt in
                                           // Einschlusskomponente
                                           // schreiben
      tmp = Inf(z) + (Sup(z)-Inf(z))/2.0; // middle(interval)
      if (abs(Inf(z)) < minreal) 
      {
         if (abs(Sup(z)) < minreal) 
         {
            if (tmp != 0.0) 
            {
               tmp = 0.0;
               n++;
            }
            i = prec;
         }
      }
          
      this->elem(i) = tmp;
      idotakku[0] -= tmp;
      z = rnd(idotakku[0]);
      i++;
   } // while
   r = Inf(z);
   if (abs(r) < minreal) 
   {
      if (r < 0.0) 
         r = -minreal;
      else         
         r = 0.0;
   }
   s = Sup(z);
   if (abs(s) < minreal) 
   {
      if (s <= 0.0) 
         s = 0.0;
      else          
         s = minreal;
   }
   if (n > 0) 
   {  
      this->elem(prec) = r-_real(n+1)*minreal;   // Intervall in die letzten
      this->elem(prec+1) = s+_real(n+1)*minreal; // Stellen schreiben
      // (n+1), da Rundungsrichtung nicht bestimmt werden kann!
   } else 
   {
      this->elem(prec) = r;           // Intervall in die letzten Stellen
      this->elem(prec+1) = s;         // schreiben
   }
}  */

void l_interval::_akku_out() throw()
{  
   // ein im idotakku[0] liegendes Zwischenergebnis
   // wird in der entsprechenden Praezision in das aufrufende l_interval
   // so gerundet, dass idotakku[0] eingeschlossen wird.
   // Neueste Version: Blomquist 21.07.03;
   _clear(1);
   interval z;
   real tmp;
   int i = 1;
   
   z = rnd(idotakku[0]);
   while (i<prec && !(CXSC_Zero <= z)) 
   {       
      // z ist Intervall, deshalb <=
      if (succ(succ(Inf(z))) < Sup(z))
         break;  // bei zu grobem Einschluss: Intervall direkt in
                 // Einschlusskomponente schreiben.
      tmp = mid(z);  // middle(interval)
      this->elem(i) = tmp;
      idotakku[0] -= tmp;
      z = rnd(idotakku[0]);
      i++;
   } // while
   this->elem(prec) = Inf(z);      // Intervall in die letzten Stellen
   this->elem(prec+1) = Sup(z);    // schreiben
} // _akku_out()

void l_interval::_akku_out_inn() throw()
{ 
   // ein im idotakku[0] liegendes Zwischenergebnis
   // wird in der entsprechenden Praezision in die aufrufende l_interval Zahl
   // li so gerundet, dass gilt:  li enthalten in idotakku[0].
   _clear(1);
   real inf, sup, tmp;          // fuer Naeherungen, entspricht Interval z
   int i=1;
   inf = rnd(Inf(idotakku[0]),RND_UP);
   sup = rnd(Sup(idotakku[0]),RND_DOWN);
   
   if (inf>sup) 
      inf = sup;      // Vorsichtsmassnahme
        
   while (i<prec && !(inf<=CXSC_Zero&&sup>=CXSC_Zero)) 
   {
      tmp = inf + (sup-inf)/2.0; // middle(interval)
      this->elem(i) = tmp;
      idotakku[0] -= tmp;
      inf = rnd(Inf(idotakku[0]),RND_UP);
      sup = rnd(Sup(idotakku[0]),RND_DOWN);
      if (inf>sup) 
         inf = sup;    // Vorsichtsmassnahme
      i++;
   }
   this->elem(prec)=inf;           // Intervall in die letzten Stellen
   this->elem(prec+1)=sup;         // schreiben
}

/* void l_interval::_akku_add(idotprecision& d) const throw()
{ 
   // addiert aufrufenden l_interval auf iakku d.
   // Änderung am 23.9.92 von F. Toussaint, da Fehler im Unterlaufbereich
   // Ausgeblendet von Blomquist, 20.11.02, da Fehler im Unterlaufbereich
   // nicht erkennbar. Meine neue Version direkt anschliessend!
   int        n = 0;
   real       r, s;

   for (int i=1; i<=prec-1; i++) 
   {
      r = this->elem(i);
      if (abs(r) < MinReal) 
      {
         if (r != 0.0) 
            n++;
      } else 
         d += _interval(r);
   }
   r = this->elem(prec);
   if (abs(r) < MinReal) 
   {
      if (r < 0.0) 
         r = -MinReal;
      else         
         r = 0.0;
   }
   s = this->elem(prec+1);
   if (abs(s) < MinReal) 
   {
      if (s <= 0.0) 
         s = 0.0;
      else          
         s = MinReal;
   }
   if (r != CXSC_Zero || s != CXSC_Zero) 
      d += _interval(r, s);
   if (n > 0) 
   {
      d += _interval(-_real(n+1)*MinReal, _real(n+1)*MinReal);
      // (n+1), da Rundungsrichtung nicht bestimmt werden kann!
   }
} */

void l_interval::_akku_add(idotprecision& d) const throw()
{ 
    // Addiert aufrufendes Intervall vom Typ l_interval auf d.
    // Meine neue Version; Blomquist, 20.11.02;
   real r,s;
   for (int i=1; i<=prec-1; i++)
   {   
       r = this->elem(i);
       if (sign(r) != CXSC_Zero) // Addition nur, falls nötig!
	   d += _interval(r); 
   }
   r = this->elem(prec);
   s = this->elem(prec+1);
   if (r != CXSC_Zero || s != CXSC_Zero) // Addition nur, falls nötig!
      d += _interval(r, s);
}

/* void l_interval::_akku_sub(idotprecision& d) const throw()
{ 
   // Subtrahiert aufrufendes Intervall vom Typ l_interval von d.
   // Intervallsubtraktion!!
   // Änderung am 23.9.92 von F. Toussaint, da Fehler im Unterlaufbereich
   // Blomquist: Mir sind keine solchen Fehler bekannt, daher wurde diese
   // Version von mir ausgeklammert, 20.11.02; Neue Version anschliessend!

   int        n = 0;
   real       r, s;

   for (int i=1; i<=prec-1; i++) 
   {
      r = this->elem(i);
      if (abs(r) < MinReal) 
      {
         if (r != 0.0) 
            n++;
      } else 
         d -= _interval(r);
   }
   r = this->elem(prec);
   if (abs(r) < MinReal) 
   {
      if (r < 0.0) 
         r = -MinReal;
      else         
         r = 0.0;
   }
   s = this->elem(prec+1);
   if (abs(s) < MinReal) 
   {
      if (s <= 0.0) 
         s = 0.0;
      else          
         s = MinReal;
   }
   d -= _interval(r, s);
   if (n > 0) 
   {
      d -= _interval(-_real(n+1)*MinReal, _real(n+1)*MinReal);
      // (n+1), da Rundungsrichtung nicht bestimmt werden kann!
   }
} */

void l_interval::_akku_sub(idotprecision& d) const throw()
{ 
    // Subtrahiert aufrufendes Intervall vom Typ l_interval von d.
    // Meine neue Version; Blomquist, 20.11.02;
   real r,s;

   for (int i=1; i<=prec-1; i++)
   { 
       r = this->elem(i);
       if (sign(r) != CXSC_Zero) // Subtraktion nur, falls nötig!
	   d -= _interval(r); 
   }
   r = this->elem(prec);
   s = this->elem(prec+1);
   if (r != CXSC_Zero || s != CXSC_Zero) // Subtraktion nur, falls nötig!
      d -= _interval(r, s);
}

// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const l_interval & a) throw()
{
   idotakku[0]=0.0;
   a._akku_add(idotakku[0]);
   s << idotakku[0];
   return s;
}
std::string & operator << (std::string &s, const l_interval & a) throw()
{
   idotakku[0]=0.0;
   a._akku_add(idotakku[0]);
   s << idotakku[0];
   return s;
}

std::istream & operator >> (std::istream &s, l_interval & a) throw()
{
   s >> idotakku[0];
   a._akku_out();
   return s;
}

std::string & operator >> (std::string &s, l_interval & a) throw()
{
   s >> idotakku[0];
   a._akku_out();
   return s;

}

void operator >>(const std::string &s,l_interval & a) throw()
{
   std::string r(s);
   r >> idotakku[0];
   a._akku_out();
}
void operator >>(const char *s,l_interval & a) throw()
{
   std::string r(s);
   r>>idotakku[0];
   a._akku_out();
}

int in ( const real& x, const l_interval& y )         // Contained-in relation
  { return ( (Inf(y) <= x) && (x <= Sup(y)) ); }      //----------------------

int in ( const l_real& x, const l_interval& y )       // Contained-in relation
  { return ( (Inf(y) <= x) && (x <= Sup(y)) ); }      //----------------------

int in ( const interval& x, const l_interval& y )     // Contained-in relation
{                                                     //----------------------
  return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) );
}

int in ( const l_interval& x, const l_interval& y )   // Contained-in relation
{                                                     //----------------------
  return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) );
}

l_interval Blow (const l_interval& x, const real& eps )
{  
    int n;
    l_interval y;
    l_real lr1,lr2;
    y = x + interval(-eps,eps) * diam(x);
    lr1 = Inf(y);
    n = StagPrec(lr1);
    lr1[n] = pred(lr1[n]);
    lr2 = Sup(y);
    lr2[n] = succ(lr2[n]);
    return l_interval(lr1,lr2);
}

int Disjoint (const l_interval& a, const l_interval& b ) 
                                                 // Test for disjointedness
{                                                //------------------------
    return (Inf(a) > Sup(b)) || (Inf(b) > Sup(a));
}

l_real AbsMin ( const l_interval& x )              // Absolute minimum of
{                                                  // an interval
   if ( in(0.0,x) ) return l_real(0.0);
   else 
   { 
       l_real y(Inf(x));
       if (y > 0.0) return y;
       else return -Sup(x);
   }

}

l_real AbsMax (const l_interval& x )            // Absolute maximum of
{                                               // an interval
    l_real a, b;                                //--------------------

    a = abs(Inf(x));
    b = abs(Sup(x));

    if (a > b)
	return a;
    else
	return b;
}

l_real RelDiam ( const l_interval x )                     // Relative diameter
{                                                         // of an interval
  if ( in(0.0,x) )                                        //------------------
    return diam(x);
  else
    return( Sup( l_interval(diam(x)) / l_interval(AbsMin(x)) ) );
}


} // namespace cxsc

