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

/* CVS $Id: xscclass.hpp,v 1.13 2007/05/16 14:59:05 cxsc Exp $ */

#ifndef XSCCLASS_HPP_INCLUDED
#define XSCCLASS_HPP_INCLUDED

#include "compiler.h"

// #ifndef IndCheck // 4.10.00
// #  define IndCheck 1
// #endif
#  define IndCheck 0
#ifdef OLD_CXSC

//#define index CXSC_index


namespace cxsc {

class index   // for backwards compatibility
{
private:
  int ind;
public:
  index() {}
  index(const int i) { ind=i; }
  int _int() const { return ind; }
};

} // namespace cxsc 

//typedef class _index index;

#endif


#ifdef _CXSC_INLINE
#  define INLINE inline
#  define TINLINE inline   // templates in vector.inl, matrix.inl
#  define _CXSC_INCL_INL
#  define _CXSC_FRIEND_TPL // define friends in *vector.hpp, *matrix.hpp
#else
#  define INLINE
#  define TINLINE
#endif

#ifdef _CXSC_CPP
#  undef INLINE
#  define INLINE
#  undef TINLINE
#  define TINLINE inline
#  define _CXSC_FRIEND_TPL
#  undef _CXSC_INCL_INL
#endif

// wegen gcc 2.95.1 <----------------------------------------================
// #undef TINLINE
// #define TINLINE

#include <string>

namespace cxsc {

using std::string; // this does the job for all files it is included

class real;             //     real.hpp
class interval;         // interval.hpp
class complex;          //  complex.hpp
class cinterval;        // cinterva.hpp   complex interval

class l_real;           //   l_real.hpp   staggered real
class l_interval;       // l_interv.hpp   staggered interval
class l_complex;
class l_cinterval;

class dotprecision;     //      dot.hpp   
class idotprecision;    //     idot.hpp   interval dotprecision
class cdotprecision;    //     cdot.hpp   complex dotprecision
class cidotprecision;   //    cidot.hpp   complex interval dotprecision

class intvector;
class intvector_slice;
class rvector;
class rvector_slice;
class ivector;
class ivector_slice;
class cvector;
class cvector_slice;
class civector;
class civector_slice;
class l_rvector;
class l_rvector_slice;
class l_ivector;
class l_ivector_slice;

class intmatrix;
class intmatrix_slice;
class intmatrix_subv;
class rmatrix;
class rmatrix_slice;
class rmatrix_subv;
class imatrix;
class imatrix_slice;
class imatrix_subv;
class cmatrix;
class cmatrix_slice;
class cmatrix_subv;
class cimatrix;
class cimatrix_slice;
class cimatrix_subv;
class l_rmatrix;
class l_rmatrix_slice;
class l_rmatrix_subv;
class l_imatrix;
class l_imatrix_slice;
class l_imatrix_subv;

inline string nameof(bool)    { return "bool"; }
inline string nameof(char)    { return "char"; }
inline string nameof(int)     { return "int"; }
inline string nameof(long)    { return "long"; }

inline string nameof(const float &)   { return "float"; }
inline string nameof(const double &)  { return "double"; }

inline string nameof(const real &)     { return "real"; }
inline string nameof(const interval &) { return "interval"; }
inline string nameof(const complex &)  { return "complex"; }
inline string nameof(const cinterval &){ return "cinterval"; }

inline string nameof(const l_real &)      { return "l_real"; }
inline string nameof(const l_interval &)  { return "l_interval"; }
inline string nameof(const l_complex &)   { return "l_complex"; }
inline string nameof(const l_cinterval &) { return "l_cinterval"; }

inline string nameof(const dotprecision &)   { return "dotprecision"; }
inline string nameof(const idotprecision &)  { return "idotprecision"; }
inline string nameof(const cdotprecision &)  { return "cdotprecision"; }
inline string nameof(const cidotprecision &) { return "cidotprecision"; }

inline string nameof(const intvector &) { return "intvector"; }
inline string nameof(const intvector_slice &) { return "intvector_slice"; }
inline string nameof(const rvector &) { return "rvector"; }
inline string nameof(const rvector_slice &) { return "rvector_slice"; }
inline string nameof(const ivector &) { return "ivector"; }
inline string nameof(const ivector_slice &) { return "ivector_slice"; }
inline string nameof(const cvector &) { return "cvector"; }
inline string nameof(const cvector_slice &) { return "cvector_slice"; }
inline string nameof(const civector &) { return "civector"; }
inline string nameof(const civector_slice &) { return "civector_slice"; }
inline string nameof(const l_rvector &) { return "l_rvector"; }
inline string nameof(const l_rvector_slice &) { return "l_rvector_slice"; }
inline string nameof(const l_ivector &) { return "l_ivector"; }
inline string nameof(const l_ivector_slice &) { return "l_ivector_slice"; }

inline string nameof(const intmatrix &) { return "intmatrix"; }
inline string nameof(const intmatrix_slice &) { return "intmatrix_slice"; }
inline string nameof(const rmatrix &) { return "rmatrix"; }
inline string nameof(const rmatrix_slice &) { return "rmatrix_slice"; }
inline string nameof(const rmatrix_subv &) { return "rmatrix_subv"; }
inline string nameof(const imatrix &) { return "imatrix"; }
inline string nameof(const imatrix_slice &) { return "imatrix_slice"; }
inline string nameof(const imatrix_subv &) { return "imatrix_subv"; }
inline string nameof(const cmatrix &) { return "cmatrix"; }
inline string nameof(const cmatrix_slice &) { return "cmatrix_slice"; }
inline string nameof(const cmatrix_subv &) { return "cmatrix_subv"; }
inline string nameof(const cimatrix &) { return "cimatrix"; }
inline string nameof(const cimatrix_slice &) { return "cimatrix_slice"; }
inline string nameof(const cimatrix_subv &) { return "cimatrix_subv"; }
inline string nameof(const l_rmatrix &) { return "l_rmatrix"; }
inline string nameof(const l_rmatrix_slice &) { return "l_rmatrix_slice"; }
inline string nameof(const l_rmatrix_subv &) { return "l_rmatrix_subv"; }
inline string nameof(const l_imatrix &) { return "l_imatrix"; }
inline string nameof(const l_imatrix_slice &) { return "l_imatrix_slice"; }
inline string nameof(const l_imatrix_subv &) { return "l_imatrix_subv"; }

} // namespace cxsc 

#endif

