/********************************************************************
AP Library version 1.2.1

Copyright (c) 2003-2007, Sergey Bochkanov (ALGLIB project).
See www.alglib.net or alglib.sources.ru for details.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************/

#ifndef AP_H
#define AP_H

#include <stdlib.h>
#include <string>
#include <math.h>

/********************************************************************
Array bounds check
********************************************************************/
#define AP_ASSERT

#ifndef AP_ASSERT     //
#define NO_AP_ASSERT  // This code avoids definition of the
#endif                // both AP_ASSERT and NO_AP_ASSERT symbols
#ifdef NO_AP_ASSERT   //
#ifdef AP_ASSERT      //
#undef NO_AP_ASSERT   //
#endif                //
#endif                //

/********************************************************************
This symbol is used for debugging. Do not define it and do not remove
comments.
********************************************************************/
//#define UNSAFE_MEM_COPY


/********************************************************************
Namespace of a standard library AlgoPascal.
********************************************************************/
namespace ap
{

/********************************************************************
Service routines:
    vlen - just alias for n2-n1+1
********************************************************************/
inline int vlen(int n1, int n2){
    return n2-n1+1;
}

/********************************************************************
Exception class.
********************************************************************/
class ap_error
{
public:
    ap_error(){};
    ap_error(const char *s){ msg = s; };

    std::string msg;

    static void make_assertion(bool bClause)
        { if(!bClause) throw ap_error(); };
    static void make_assertion(bool bClause, const char *msg)
        { if(!bClause) throw ap_error(msg); };
private:
};


/********************************************************************
Templates for vector operations
********************************************************************/
#include "apvt.h"

/********************************************************************
BLAS functions
********************************************************************/
template<class T> T vdotproduct(const T *v1, const T *v2, int N);
template<class T> void vmove(T *vdst, const T* vsrc, int N);
template<class T> void vmoveneg(T *vdst, const T *vsrc, int N);
template<class T> void vmove(T *vdst, const T *vsrc, int N, T alpha);
template<class T> void vadd(T *vdst, const T *vsrc, int N);
template<class T> void vadd(T *vdst, const T *vsrc, int N, T alpha);
template<class T> void vsub(T *vdst, const T *vsrc, int N);
template<class T> void vsub(T *vdst, const T *vsrc, int N, T alpha);
template<class T> void vmul(T *vdst, int N, T alpha);

/********************************************************************
Template of a dynamical one-dimensional array
********************************************************************/
template<class T>
class real_1d_array
{
public:
    real_1d_array()
    {
        m_Vec=0;
        m_iVecSize = 0;
        m_iLow = 0;
        m_iHigh = -1;
    };

    ~real_1d_array()
    {
        if(m_Vec)
        {
                delete[] m_Vec;
        }
    };

    real_1d_array(const real_1d_array &rhs)
    {
        m_Vec=0;
        m_iVecSize = 0;
        m_iLow = 0;
        m_iHigh = -1;
        if( rhs.m_iVecSize!=0 )
            setcontent(rhs.m_iLow, rhs.m_iHigh, rhs.getcontent());
    };


    const real_1d_array& operator=(const real_1d_array &rhs)
    {
        if( this==&rhs )
            return *this;

        if( rhs.m_iVecSize!=0 )
            setcontent(rhs.m_iLow, rhs.m_iHigh, rhs.getcontent());
        else
        {
            m_Vec=0;
            m_iVecSize = 0;
            m_iLow = 0;
            m_iHigh = -1;
        }
        return *this;
    };


    const T& operator()(int i) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i>=m_iLow && i<=m_iHigh);
        #endif
        return m_Vec[ i-m_iLow ];
    };


    T& operator()(int i)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i>=m_iLow && i<=m_iHigh);
        #endif
        return m_Vec[ i-m_iLow ];
    };


    void setbounds( int iLow, int iHigh )
    {
        if(m_Vec)
        {
                delete[] m_Vec;
        }
        m_iLow = iLow;
        m_iHigh = iHigh;
        m_iVecSize = iHigh-iLow+1;
        m_Vec = new T[m_iVecSize];
    };


    void setcontent( int iLow, int iHigh, const T *pContent )
    {
        setbounds(iLow, iHigh);
        for(int i=0; i<m_iVecSize; i++)
            m_Vec[i] = pContent[i];
    };


    T* getcontent()
    {
        return m_Vec;
    };

    const T* getcontent() const
    {
        return m_Vec;
    };


    int getlowbound(int iBoundNum = 0) const
    {
        return m_iLow;
    };


    int gethighbound(int iBoundNum = 0) const
    {
        return m_iHigh;
    };

    raw_vector<T> getvector(int iStart, int iEnd)
    {
        if( iStart>iEnd || wrongIdx(iStart) || wrongIdx(iEnd) )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(m_Vec+iStart-m_iLow, iEnd-iStart+1, 1);
    };


    const_raw_vector<T> getvector(int iStart, int iEnd) const
    {
        if( iStart>iEnd || wrongIdx(iStart) || wrongIdx(iEnd) )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(m_Vec+iStart-m_iLow, iEnd-iStart+1, 1);
    };
private:
    bool wrongIdx(int i) const { return i<m_iLow || i>m_iHigh; };

    T         *m_Vec;
    long      m_iVecSize;
    long      m_iLow, m_iHigh;
};



/********************************************************************
Template of a dynamical two-dimensional array
********************************************************************/
template<class T>
class real_2d_array
{
public:
    real_2d_array()
    {
        m_Vec=0;
        m_iVecSize=0;
        m_iLow1 = 0;
        m_iHigh1 = -1;
        m_iLow2 = 0;
        m_iHigh2 = -1;
    };

    ~real_2d_array()
    {
        if(m_Vec)
        {
            delete[] m_Vec;
        }
    };

    real_2d_array(const real_2d_array &rhs)
    {
        m_Vec=0;
        m_iVecSize=0;
        m_iLow1 = 0;
        m_iHigh1 = -1;
        m_iLow2 = 0;
        m_iHigh2 = -1;
        if( rhs.m_iVecSize!=0 )
        {
            setbounds(rhs.m_iLow1, rhs.m_iHigh1, rhs.m_iLow2, rhs.m_iHigh2);
            for(int i=m_iLow1; i<=m_iHigh1; i++)
                vmove(&(operator()(i,m_iLow2)), &(rhs(i,m_iLow2)), m_iHigh2-m_iLow2+1);
        }
    };
    const real_2d_array& operator=(const real_2d_array &rhs)
    {
        if( this==&rhs )
            return *this;

        if( rhs.m_iVecSize!=0 )
        {
            setbounds(rhs.m_iLow1, rhs.m_iHigh1, rhs.m_iLow2, rhs.m_iHigh2);
            for(int i=m_iLow1; i<=m_iHigh1; i++)
                vmove(&(operator()(i,m_iLow2)), &(rhs(i,m_iLow2)), m_iHigh2-m_iLow2+1);
        }
        else
        {
            m_Vec=0;
            m_iVecSize=0;
            m_iLow1 = 0;
            m_iHigh1 = -1;
            m_iLow2 = 0;
            m_iHigh2 = -1;
        }
        return *this;
    };

    const T& operator()(int i1, int i2) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i1>=m_iLow1 && i1<=m_iHigh1);
        ap_error::make_assertion(i2>=m_iLow2 && i2<=m_iHigh2);
        #endif
        return m_Vec[ m_iConstOffset + i2 +i1*m_iLinearMember];
    };

    T& operator()(int i1, int i2)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i1>=m_iLow1 && i1<=m_iHigh1);
        ap_error::make_assertion(i2>=m_iLow2 && i2<=m_iHigh2);
        #endif
        return m_Vec[ m_iConstOffset + i2 +i1*m_iLinearMember];
    };

    void setbounds( int iLow1, int iHigh1, int iLow2, int iHigh2 )
    {
        if(m_Vec)
        {
            delete[] m_Vec;
        }
        int n1 = iHigh1-iLow1+1;
        int n2 = iHigh2-iLow2+1;
        m_iVecSize = n1*n2;
        m_Vec = new T[m_iVecSize];
        m_iLow1  = iLow1;
        m_iHigh1 = iHigh1;
        m_iLow2  = iLow2;
        m_iHigh2 = iHigh2;
        m_iConstOffset = -m_iLow2-m_iLow1*n2;
        m_iLinearMember = n2;
    };

    void setcontent( int iLow1, int iHigh1, int iLow2, int iHigh2, const T *pContent )
    {
        setbounds(iLow1, iHigh1, iLow2, iHigh2);
        for(int i=m_iLow1; i<=m_iHigh1; i++, pContent += m_iHigh2-m_iLow2+1)
            vmove(&(operator()(i,m_iLow2)), pContent, m_iHigh2-m_iLow2+1);
    };

    int getlowbound(int iBoundNum) const
    {
        return iBoundNum==1 ? m_iLow1 : m_iLow2;
    };

    int gethighbound(int iBoundNum) const
    {
        return iBoundNum==1 ? m_iHigh1 : m_iHigh2;
    };

    raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd)
    {
        if( (iRowStart>iRowEnd) || wrongColumn(iColumn) || wrongRow(iRowStart) ||wrongRow(iRowEnd) )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(&((*this)(iRowStart, iColumn)), iRowEnd-iRowStart+1, m_iLinearMember);
    };

    raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd)
    {
        if( (iColumnStart>iColumnEnd) || wrongRow(iRow) || wrongColumn(iColumnStart) || wrongColumn(iColumnEnd))
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(&((*this)(iRow, iColumnStart)), iColumnEnd-iColumnStart+1, 1);
    };

    const_raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd) const
    {
        if( (iRowStart>iRowEnd) || wrongColumn(iColumn) || wrongRow(iRowStart) ||wrongRow(iRowEnd) )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(&((*this)(iRowStart, iColumn)), iRowEnd-iRowStart+1, m_iLinearMember);
    };

    const_raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd) const
    {
        if( (iColumnStart>iColumnEnd) || wrongRow(iRow) || wrongColumn(iColumnStart) || wrongColumn(iColumnEnd))
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(&((*this)(iRow, iColumnStart)), iColumnEnd-iColumnStart+1, 1);
    };
private:
    bool wrongRow(int i) const { return i<m_iLow1 || i>m_iHigh1; };
    bool wrongColumn(int j) const { return j<m_iLow2 || j>m_iHigh2; };

    T           *m_Vec;
    long        m_iVecSize;
    long        m_iLow1, m_iLow2, m_iHigh1, m_iHigh2;
    long        m_iConstOffset, m_iLinearMember;
};


/********************************************************************
Constants and functions introduced for compatibility with AlgoPascal
********************************************************************/

template<class T> T sqr(T x);
template<class T> int maxint(int m1, int m2);
template<class T> int minint(int m1, int m2);
template<class T> T maxreal(T m1, T m2);
template<class T> T minreal(T m1, T m2);

};//namespace ap

/********************************************************************
BLAS functions
********************************************************************/
template<class T>
T ap::vdotproduct(const T *v1, const T *v2, int N)
{
    return ap::_vdotproduct<T>(v1, v2, N);
}

template<class T>
void ap::vmove(T *vdst, const T* vsrc, int N)
{
    ap::_vmove<T>(vdst, vsrc, N);
}

template<class T>
void ap::vmoveneg(T *vdst, const T *vsrc, int N)
{
    ap::_vmoveneg<T>(vdst, vsrc, N);
}

template<class T>
void ap::vmove(T *vdst, const T *vsrc, int N, T alpha)
{
    ap::_vmove<T,T>(vdst, vsrc, N, alpha);
}

template<class T>
void ap::vadd(T *vdst, const T *vsrc, int N)
{
    ap::_vadd<T>(vdst, vsrc, N);
}

template<class T>
void ap::vadd(T *vdst, const T *vsrc, int N, T alpha)
{
    ap::_vadd<T,T>(vdst, vsrc, N, alpha);
}

template<class T>
void ap::vsub(T *vdst, const T *vsrc, int N)
{
    ap::_vsub<T>(vdst, vsrc, N);
}

template<class T>
void ap::vsub(T *vdst, const T *vsrc, int N, T alpha)
{
    ap::_vsub<T,T>(vdst, vsrc, N, alpha);
}

template<class T>
void ap::vmul(T *vdst, int N, T alpha)
{
    ap::_vmul<T,T>(vdst, N, alpha);
}

/********************************************************************
standard functions
********************************************************************/

template<class T>
T ap::sqr(T x)
{ return x*x; }

template<class T>
int ap::maxint(int m1, int m2)
{
    return m1>m2 ? m1 : m2;
}

template<class T>
int ap::minint(int m1, int m2)
{
    return m1>m2 ? m2 : m1;
}

template<class T>
T ap::maxreal(T m1, T m2)
{
    return m1>m2 ? m1 : m2;
}

template<class T>
T ap::minreal(T m1, T m2)
{
    return m1>m2 ? m2 : m1;
}


#endif
