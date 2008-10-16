/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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
*************************************************************************/

#ifndef _tdevd_h
#define _tdevd_h

#include "ap.h"


template<class T>
void tdevde2(const T& a,
     const T& b,
     const T& c,
     T& rt1,
     T& rt2);
template<class T>
void tdevdev2(const T& a,
     const T& b,
     const T& c,
     T& rt1,
     T& rt2,
     T& cs1,
     T& sn1);
template<class T>
T tdevdpythag(T a, T b);
template<class T>
T tdevdextsign(T a, T b);


/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [N1..N2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
template<class T>
void applyrotationsfromtheleft(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array<T>& c,
     const ap::real_1d_array<T>& s,
     ap::real_2d_array<T>& a,
     ap::real_1d_array<T>& work)
{
    int j;
    int jp1;
    T ctemp;
    T stemp;
    T temp;

    if( m1>m2||n1>n2 )
    {
        return;
    }
    
    //
    // Form  P * A
    //
    if( isforward )
    {
        if( n1!=n2 )
        {
            
            //
            // Common case: N1<>N2
            //
            for(j = m1; j <= m2-1; j++)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    jp1 = j+1;
                    ap::vmove(&work(n1), &a(jp1, n1), ap::vlen(n1,n2), ctemp);
                    ap::vsub(&work(n1), &a(j, n1), ap::vlen(n1,n2), stemp);
                    ap::vmul(&a(j, n1), ap::vlen(n1,n2), ctemp);
                    ap::vadd(&a(j, n1), &a(jp1, n1), ap::vlen(n1,n2), stemp);
                    ap::vmove(&a(jp1, n1), &work(n1), ap::vlen(n1,n2));
                }
            }
        }
        else
        {
            
            //
            // Special case: N1=N2
            //
            for(j = m1; j <= m2-1; j++)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    temp = a(j+1,n1);
                    a(j+1,n1) = ctemp*temp-stemp*a(j,n1);
                    a(j,n1) = stemp*temp+ctemp*a(j,n1);
                }
            }
        }
    }
    else
    {
        if( n1!=n2 )
        {
            
            //
            // Common case: N1<>N2
            //
            for(j = m2-1; j >= m1; j--)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    jp1 = j+1;
                    ap::vmove(&work(n1), &a(jp1, n1), ap::vlen(n1,n2), ctemp);
                    ap::vsub(&work(n1), &a(j, n1), ap::vlen(n1,n2), stemp);
                    ap::vmul(&a(j, n1), ap::vlen(n1,n2), ctemp);
                    ap::vadd(&a(j, n1), &a(jp1, n1), ap::vlen(n1,n2), stemp);
                    ap::vmove(&a(jp1, n1), &work(n1), ap::vlen(n1,n2));
                }
            }
        }
        else
        {
            
            //
            // Special case: N1=N2
            //
            for(j = m2-1; j >= m1; j--)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    temp = a(j+1,n1);
                    a(j+1,n1) = ctemp*temp-stemp*a(j,n1);
                    a(j,n1) = stemp*temp+ctemp*a(j,n1);
                }
            }
        }
    }
}


/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [M1..M2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
template<class T>
void applyrotationsfromtheright(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array<T>& c,
     const ap::real_1d_array<T>& s,
     ap::real_2d_array<T>& a,
     ap::real_1d_array<T>& work)
{
    int j;
    int jp1;
    T ctemp;
    T stemp;
    T temp;

    
    //
    // Form A * P'
    //
    if( isforward )
    {
        if( m1!=m2 )
        {
            
            //
            // Common case: M1<>M2
            //
            for(j = n1; j <= n2-1; j++)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    jp1 = j+1;
                    ap::vmove(work.getvector(m1, m2), a.getcolumn(jp1, m1, m2), ctemp);
                    ap::vsub(work.getvector(m1, m2), a.getcolumn(j, m1, m2), stemp);
                    ap::vmul(a.getcolumn(j, m1, m2), ctemp);
                    ap::vadd(a.getcolumn(j, m1, m2), a.getcolumn(jp1, m1, m2), stemp);
                    ap::vmove(a.getcolumn(jp1, m1, m2), work.getvector(m1, m2));
                }
            }
        }
        else
        {
            
            //
            // Special case: M1=M2
            //
            for(j = n1; j <= n2-1; j++)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    temp = a(m1,j+1);
                    a(m1,j+1) = ctemp*temp-stemp*a(m1,j);
                    a(m1,j) = stemp*temp+ctemp*a(m1,j);
                }
            }
        }
    }
    else
    {
        if( m1!=m2 )
        {
            
            //
            // Common case: M1<>M2
            //
            for(j = n2-1; j >= n1; j--)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    jp1 = j+1;
                    ap::vmove(work.getvector(m1, m2), a.getcolumn(jp1, m1, m2), ctemp);
                    ap::vsub(work.getvector(m1, m2), a.getcolumn(j, m1, m2), stemp);
                    ap::vmul(a.getcolumn(j, m1, m2), ctemp);
                    ap::vadd(a.getcolumn(j, m1, m2), a.getcolumn(jp1, m1, m2), stemp);
                    ap::vmove(a.getcolumn(jp1, m1, m2), work.getvector(m1, m2));
                }
            }
        }
        else
        {
            
            //
            // Special case: M1=M2
            //
            for(j = n2-1; j >= n1; j--)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    temp = a(m1,j+1);
                    a(m1,j+1) = ctemp*temp-stemp*a(m1,j);
                    a(m1,j) = stemp*temp+ctemp*a(m1,j);
                }
            }
        }
    }
}


/*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************/
template<class T>
void generaterotation(T f, T g, T& cs, T& sn, T& r)
{
    T f1;
    T g1;

    if( g==0 )
    {
        cs = 1;
        sn = T(0.0);
        r = f;
    }
    else
    {
        if( f==0 )
        {
            cs = T(0.0);
            sn = 1;
            r = g;
        }
        else
        {
            f1 = f;
            g1 = g;
            r = sqrt(ap::sqr(f1)+ap::sqr(g1));
            cs = f1/r;
            sn = g1/r;
            if( fabs(f)>fabs(g)&&cs<0 )
            {
                cs = -cs;
                sn = -sn;
                r = -r;
            }
        }
    }
}


//
// Blas routine inplacetranspose
//
template<class T>
void inplacetranspose(ap::real_2d_array<T>& a,
     int i1,
     int i2,
     int j1,
     int j2,
     ap::real_1d_array<T>& work)
{
    int i;
    int j;
    int ips;
    int jps;
    int l;

    if( i1>i2||j1>j2 )
    {
        return;
    }
    ap::ap_error::make_assertion(i1-i2==j1-j2, "InplaceTranspose error: incorrect array size!");
    for(i = i1; i <= i2-1; i++)
    {
        j = j1+i-i1;
        ips = i+1;
        jps = j1+ips-i1;
        l = i2-i;
        ap::vmove(work.getvector(1, l), a.getcolumn(j, ips, i2));
        ap::vmove(a.getcolumn(j, ips, i2), a.getrow(i, jps, j2));
        ap::vmove(&a(i, jps), &work(1), ap::vlen(jps,j2));
    }
}


/*************************************************************************
Finding the eigenvalues and eigenvectors of a tridiagonal symmetric matrix

The algorithm finds the eigen pairs of a tridiagonal symmetric matrix by
using an QL/QR algorithm with implicit shifts.

Input parameters:
    D       -   the main diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-1].
    E       -   the secondary diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-2].
    N       -   size of matrix A.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not needed;
                 * 1, the eigenvectors of a tridiagonal matrix
                   are multiplied by the square matrix Z. It is used if the
                   tridiagonal matrix is obtained by the similarity
                   transformation of a symmetric matrix;
                 * 2, the eigenvectors of a tridiagonal matrix replace the
                   square matrix Z;
                 * 3, matrix Z contains the first row of the eigenvectors
                   matrix.
    Z       -   if ZNeeded=1, Z contains the square matrix by which the
                eigenvectors are multiplied.
                Array whose indexes range within [0..N-1, 0..N-1].

Output parameters:
    D       -   eigenvalues in ascending order.
                Array whose index ranges within [0..N-1].
    Z       -   if ZNeeded is equal to:
                 * 0, Z hasn�t changed;
                 * 1, Z contains the product of a given matrix (from the left)
                   and the eigenvectors matrix (from the right);
                 * 2, Z contains the eigenvectors.
                 * 3, Z contains the first row of the eigenvectors matrix.
                If ZNeeded<3, Z is the array whose indexes range within [0..N-1, 0..N-1].
                In that case, the eigenvectors are stored in the matrix columns.
                If ZNeeded=3, Z is the array whose indexes range within [0..0, 0..N-1].

Result:
    True, if the algorithm has converged.
    False, if the algorithm hasn't converged.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
template<class T>
bool smatrixtdevd(ap::real_1d_array<T>& d,
     ap::real_1d_array<T> e,
     int n,
     int zneeded,
     ap::real_2d_array<T>& z)
{
    bool result;
    ap::real_1d_array<T> d1;
    ap::real_1d_array<T> e1;
    ap::real_2d_array<T> z1;
    int i;

    
    //
    // Prepare 1-based task
    //
    d1.setbounds(1, n);
    e1.setbounds(1, n);
    ap::vmove(&d1(1), &d(0), ap::vlen(1,n));
    if( n>1 )
    {
        ap::vmove(&e1(1), &e(0), ap::vlen(1,n-1));
    }
    if( zneeded==1 )
    {
        z1.setbounds(1, n, 1, n);
        for(i = 1; i <= n; i++)
        {
            ap::vmove(&z1(i, 1), &z(i-1, 0), ap::vlen(1,n));
        }
    }
    
    //
    // Solve 1-based task
    //
    result = tridiagonalevd(d1, e1, n, zneeded, z1);
    if( !result )
    {
        return result;
    }
    
    //
    // Convert back to 0-based result
    //
    ap::vmove(&d(0), &d1(1), ap::vlen(0,n-1));
    if( zneeded!=0 )
    {
        if( zneeded==1 )
        {
            for(i = 1; i <= n; i++)
            {
                ap::vmove(&z(i-1, 0), &z1(i, 1), ap::vlen(0,n-1));
            }
            return result;
        }
        if( zneeded==2 )
        {
            z.setbounds(0, n-1, 0, n-1);
            for(i = 1; i <= n; i++)
            {
                ap::vmove(&z(i-1, 0), &z1(i, 1), ap::vlen(0,n-1));
            }
            return result;
        }
        if( zneeded==3 )
        {
            z.setbounds(0, 0, 0, n-1);
            ap::vmove(&z(0, 0), &z1(1, 1), ap::vlen(0,n-1));
            return result;
        }
        ap::ap_error::make_assertion(false, "SMatrixTDEVD: Incorrect ZNeeded!");
    }
    return result;
}

/*************************************************************************
DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
   [  A   B  ]
   [  B   C  ].
On return, RT1 is the eigenvalue of larger absolute value, and RT2
is the eigenvalue of smaller absolute value.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************/
template<class T>
void tdevde2(const T& a,
     const T& b,
     const T& c,
     T& rt1,
     T& rt2)
{
    T ab;
    T acmn;
    T acmx;
    T adf;
    T df;
    T rt;
    T sm;
    T tb;

    sm = a+c;
    df = a-c;
    adf = fabs(df);
    tb = b+b;
    ab = fabs(tb);
    if( fabs(a)>fabs(c) )
    {
        acmx = a;
        acmn = c;
    }
    else
    {
        acmx = c;
        acmn = a;
    }
    if( adf>ab )
    {
        rt = adf*sqrt(1+ap::sqr(ab/adf));
    }
    else
    {
        if( adf<ab )
        {
            rt = ab*sqrt(1+ap::sqr(adf/ab));
        }
        else
        {
            
            //
            // Includes case AB=ADF=0
            //
            rt = ab*sqrt(T(2));
        }
    }
    if( sm<0 )
    {
        rt1 = 0.5*(sm-rt);
        
        //
        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
        //
        rt2 = acmx/rt1*acmn-b/rt1*b;
    }
    else
    {
        if( sm>0 )
        {
            rt1 = 0.5*(sm+rt);
            
            //
            // Order of execution important.
            // To get fully accurate smaller eigenvalue,
            // next line needs to be executed in higher precision.
            //
            rt2 = acmx/rt1*acmn-b/rt1*b;
        }
        else
        {
            
            //
            // Includes case RT1 = RT2 = 0
            //
            rt1 = 0.5*rt;
            rt2 = -0.5*rt;
        }
    }
}

/*************************************************************************
Obsolete 1-based subroutine.
*************************************************************************/
template<class T>
bool tridiagonalevd(ap::real_1d_array<T>& d,
     ap::real_1d_array<T> e,
     int n,
     int zneeded,
     ap::real_2d_array<T>& z)
{
    bool result;
    int maxit;
    int i;
    int ii;
    int iscale;
    int j;
    int jtot;
    int k;
    int t;
    int l;
    int l1;
    int lend;
    int lendm1;
    int lendp1;
    int lendsv;
    int lm1;
    int lsv;
    int m;
    int mm;
    int mm1;
    int nm1;
    int nmaxit;
    int tmpint;
    T anorm;
    T b;
    T c;
    T eps;
    T eps2;
    T f;
    T g;
    T p;
    T r;
    T rt1;
    T rt2;
    T s;
    T safmax;
    T safmin;
    T ssfmax;
    T ssfmin;
    T tst;
    T tmp;
    ap::real_1d_array<T> work1;
    ap::real_1d_array<T> work2;
    ap::real_1d_array<T> workc;
    ap::real_1d_array<T> works;
    ap::real_1d_array<T> wtemp;
    bool gotoflag;
    int zrows;
    bool wastranspose;

    ap::ap_error::make_assertion(zneeded>=0&&zneeded<=3, "TridiagonalEVD: Incorrent ZNeeded");
    
    //
    // Quick return if possible
    //
    if( zneeded<0||zneeded>3 )
    {
        result = false;
        return result;
    }
    result = true;
    if( n==0 )
    {
        return result;
    }
    if( n==1 )
    {
        if( zneeded==2||zneeded==3 )
        {
            z.setbounds(1, 1, 1, 1);
            z(1,1) = 1;
        }
        return result;
    }
    maxit = 30;
    
    //
    // Initialize arrays
    //
    wtemp.setbounds(1, n);
    work1.setbounds(1, n-1);
    work2.setbounds(1, n-1);
    workc.setbounds(1, n);
    works.setbounds(1, n);
    
    //
    // Determine the unit roundoff and over/underflow thresholds.
    //
    eps = numeric_limits<T>::epsilon();
    eps2 = ap::sqr(eps);
    safmin = numeric_limits<T>::min();
    safmax = numeric_limits<T>::max();
    ssfmax = sqrt(safmax)/3;
    ssfmin = sqrt(safmin)/eps2;
    
    //
    // Prepare Z
    //
    // Here we are using transposition to get rid of column operations
    //
    //
    wastranspose = false;
    if( zneeded==0 )
    {
        zrows = 0;
    }
    if( zneeded==1 )
    {
        zrows = n;
    }
    if( zneeded==2 )
    {
        zrows = n;
    }
    if( zneeded==3 )
    {
        zrows = 1;
    }
    if( zneeded==1 )
    {
        wastranspose = true;
        inplacetranspose(z, 1, n, 1, n, wtemp);
    }
    if( zneeded==2 )
    {
        wastranspose = true;
        z.setbounds(1, n, 1, n);
        for(i = 1; i <= n; i++)
        {
            for(j = 1; j <= n; j++)
            {
                if( i==j )
                {
                    z(i,j) = 1;
                }
                else
                {
                    z(i,j) = T(0.0);
                }
            }
        }
    }
    if( zneeded==3 )
    {
        wastranspose = false;
        z.setbounds(1, 1, 1, n);
        for(j = 1; j <= n; j++)
        {
            if( j==1 )
            {
                z(1,j) = 1;
            }
            else
            {
                z(1,j) = T(0.0);
            }
        }
    }
    nmaxit = n*maxit;
    jtot = 0;
    
    //
    // Determine where the matrix splits and choose QL or QR iteration
    // for each block, according to whether top or bottom diagonal
    // element is smaller.
    //
    l1 = 1;
    nm1 = n-1;
    while(true)
    {
        if( l1>n )
        {
            break;
        }
        if( l1>1 )
        {
            e(l1-1) = T(0.0);
        }
        gotoflag = false;
        if( l1<=nm1 )
        {
            for(m = l1; m <= nm1; m++)
            {
                tst = fabs(e(m));
                if( tst==0 )
                {
                    gotoflag = true;
                    break;
                }
                if( tst<=sqrt(fabs(d(m)))*sqrt(fabs(d(m+1)))*eps )
                {
                    e(m) = T(0.0);
                    gotoflag = true;
                    break;
                }
            }
        }
        if( !gotoflag )
        {
            m = n;
        }
        
        //
        // label 30:
        //
        l = l1;
        lsv = l;
        lend = m;
        lendsv = lend;
        l1 = m+1;
        if( lend==l )
        {
            continue;
        }
        
        //
        // Scale submatrix in rows and columns L to LEND
        //
        if( l==lend )
        {
            anorm = fabs(d(l));
        }
        else
        {
            anorm = ap::maxreal(fabs(d(l))+fabs(e(l)), fabs(e(lend-1))+fabs(d(lend)));
            for(i = l+1; i <= lend-1; i++)
            {
                anorm = ap::maxreal(anorm, fabs(d(i))+fabs(e(i))+fabs(e(i-1)));
            }
        }
        iscale = 0;
        if( anorm==0 )
        {
            continue;
        }
        if( anorm>ssfmax )
        {
            iscale = 1;
            tmp = ssfmax/anorm;
            tmpint = lend-1;
            ap::vmul(&d(l), ap::vlen(l,lend), tmp);
            ap::vmul(&e(l), ap::vlen(l,tmpint), tmp);
        }
        if( anorm<ssfmin )
        {
            iscale = 2;
            tmp = ssfmin/anorm;
            tmpint = lend-1;
            ap::vmul(&d(l), ap::vlen(l,lend), tmp);
            ap::vmul(&e(l), ap::vlen(l,tmpint), tmp);
        }
        
        //
        // Choose between QL and QR iteration
        //
        if( fabs(d(lend))<fabs(d(l)) )
        {
            lend = lsv;
            l = lendsv;
        }
        if( lend>l )
        {
            
            //
            // QL Iteration
            //
            // Look for small subdiagonal element.
            //
            while(true)
            {
                gotoflag = false;
                if( l!=lend )
                {
                    lendm1 = lend-1;
                    for(m = l; m <= lendm1; m++)
                    {
                        tst = ap::sqr(fabs(e(m)));
                        if( tst<=eps2*fabs(d(m))*fabs(d(m+1))+safmin )
                        {
                            gotoflag = true;
                            break;
                        }
                    }
                }
                if( !gotoflag )
                {
                    m = lend;
                }
                if( m<lend )
                {
                    e(m) = T(0.0);
                }
                p = d(l);
                if( m!=l )
                {
                    
                    //
                    // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
                    // to compute its eigensystem.
                    //
                    if( m==l+1 )
                    {
                        if( zneeded>0 )
                        {
                            tdevdev2(d(l), e(l), d(l+1), rt1, rt2, c, s);
                            work1(l) = c;
                            work2(l) = s;
                            workc(1) = work1(l);
                            works(1) = work2(l);
                            if( !wastranspose )
                            {
                                applyrotationsfromtheright(false, 1, zrows, l, l+1, workc, works, z, wtemp);
                            }
                            else
                            {
                                applyrotationsfromtheleft(false, l, l+1, 1, zrows, workc, works, z, wtemp);
                            }
                        }
                        else
                        {
                            tdevde2(d(l), e(l), d(l+1), rt1, rt2);
                        }
                        d(l) = rt1;
                        d(l+1) = rt2;
                        e(l) = T(0.0);
                        l = l+2;
                        if( l<=lend )
                        {
                            continue;
                        }
                        
                        //
                        // GOTO 140
                        //
                        break;
                    }
                    if( jtot==nmaxit )
                    {
                        
                        //
                        // GOTO 140
                        //
                        break;
                    }
                    jtot = jtot+1;
                    
                    //
                    // Form shift.
                    //
                    g = (d(l+1)-p)/(2*e(l));
                    r = tdevdpythag(g, T(1));
                    g = d(m)-p+e(l)/(g+tdevdextsign(r, g));
                    s = 1;
                    c = 1;
                    p = T(0.0);
                    
                    //
                    // Inner loop
                    //
                    mm1 = m-1;
                    for(i = mm1; i >= l; i--)
                    {
                        f = s*e(i);
                        b = c*e(i);
                        generaterotation(g, f, c, s, r);
                        if( i!=m-1 )
                        {
                            e(i+1) = r;
                        }
                        g = d(i+1)-p;
                        r = (d(i)-g)*s+2*c*b;
                        p = s*r;
                        d(i+1) = g+p;
                        g = c*r-b;
                        
                        //
                        // If eigenvectors are desired, then save rotations.
                        //
                        if( zneeded>0 )
                        {
                            work1(i) = c;
                            work2(i) = -s;
                        }
                    }
                    
                    //
                    // If eigenvectors are desired, then apply saved rotations.
                    //
                    if( zneeded>0 )
                    {
                        for(i = l; i <= m-1; i++)
                        {
                            workc(i-l+1) = work1(i);
                            works(i-l+1) = work2(i);
                        }
                        if( !wastranspose )
                        {
                            applyrotationsfromtheright(false, 1, zrows, l, m, workc, works, z, wtemp);
                        }
                        else
                        {
                            applyrotationsfromtheleft(false, l, m, 1, zrows, workc, works, z, wtemp);
                        }
                    }
                    d(l) = d(l)-p;
                    e(l) = g;
                    continue;
                }
                
                //
                // Eigenvalue found.
                //
                d(l) = p;
                l = l+1;
                if( l<=lend )
                {
                    continue;
                }
                break;
            }
        }
        else
        {
            
            //
            // QR Iteration
            //
            // Look for small superdiagonal element.
            //
            while(true)
            {
                gotoflag = false;
                if( l!=lend )
                {
                    lendp1 = lend+1;
                    for(m = l; m >= lendp1; m--)
                    {
                        tst = ap::sqr(fabs(e(m-1)));
                        if( tst<=eps2*fabs(d(m))*fabs(d(m-1))+safmin )
                        {
                            gotoflag = true;
                            break;
                        }
                    }
                }
                if( !gotoflag )
                {
                    m = lend;
                }
                if( m>lend )
                {
                    e(m-1) = T(0.0);
                }
                p = d(l);
                if( m!=l )
                {
                    
                    //
                    // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
                    // to compute its eigensystem.
                    //
                    if( m==l-1 )
                    {
                        if( zneeded>0 )
                        {
                            tdevdev2(d(l-1), e(l-1), d(l), rt1, rt2, c, s);
                            work1(m) = c;
                            work2(m) = s;
                            workc(1) = c;
                            works(1) = s;
                            if( !wastranspose )
                            {
                                applyrotationsfromtheright(true, 1, zrows, l-1, l, workc, works, z, wtemp);
                            }
                            else
                            {
                                applyrotationsfromtheleft(true, l-1, l, 1, zrows, workc, works, z, wtemp);
                            }
                        }
                        else
                        {
                            tdevde2(d(l-1), e(l-1), d(l), rt1, rt2);
                        }
                        d(l-1) = rt1;
                        d(l) = rt2;
                        e(l-1) = T(0.0);
                        l = l-2;
                        if( l>=lend )
                        {
                            continue;
                        }
                        break;
                    }
                    if( jtot==nmaxit )
                    {
                        break;
                    }
                    jtot = jtot+1;
                    
                    //
                    // Form shift.
                    //
                    g = (d(l-1)-p)/(2*e(l-1));
                    r = tdevdpythag(g, T(1));
                    g = d(m)-p+e(l-1)/(g+tdevdextsign(r, g));
                    s = 1;
                    c = 1;
                    p = T(0.0);
                    
                    //
                    // Inner loop
                    //
                    lm1 = l-1;
                    for(i = m; i <= lm1; i++)
                    {
                        f = s*e(i);
                        b = c*e(i);
                        generaterotation(g, f, c, s, r);
                        if( i!=m )
                        {
                            e(i-1) = r;
                        }
                        g = d(i)-p;
                        r = (d(i+1)-g)*s+2*c*b;
                        p = s*r;
                        d(i) = g+p;
                        g = c*r-b;
                        
                        //
                        // If eigenvectors are desired, then save rotations.
                        //
                        if( zneeded>0 )
                        {
                            work1(i) = c;
                            work2(i) = s;
                        }
                    }
                    
                    //
                    // If eigenvectors are desired, then apply saved rotations.
                    //
                    if( zneeded>0 )
                    {
                        mm = l-m+1;
                        for(i = m; i <= l-1; i++)
                        {
                            workc(i-m+1) = work1(i);
                            works(i-m+1) = work2(i);
                        }
                        if( !wastranspose )
                        {
                            applyrotationsfromtheright(true, 1, zrows, m, l, workc, works, z, wtemp);
                        }
                        else
                        {
                            applyrotationsfromtheleft(true, m, l, 1, zrows, workc, works, z, wtemp);
                        }
                    }
                    d(l) = d(l)-p;
                    e(lm1) = g;
                    continue;
                }
                
                //
                // Eigenvalue found.
                //
                d(l) = p;
                l = l-1;
                if( l>=lend )
                {
                    continue;
                }
                break;
            }
        }
        
        //
        // Undo scaling if necessary
        //
        if( iscale==1 )
        {
            tmp = anorm/ssfmax;
            tmpint = lendsv-1;
            ap::vmul(&d(lsv), ap::vlen(lsv,lendsv), tmp);
            ap::vmul(&e(lsv), ap::vlen(lsv,tmpint), tmp);
        }
        if( iscale==2 )
        {
            tmp = anorm/ssfmin;
            tmpint = lendsv-1;
            ap::vmul(&d(lsv), ap::vlen(lsv,lendsv), tmp);
            ap::vmul(&e(lsv), ap::vlen(lsv,tmpint), tmp);
        }
        
        //
        // Check for no convergence to an eigenvalue after a total
        // of N*MAXIT iterations.
        //
        if( jtot>=nmaxit )
        {
            result = false;
            if( wastranspose )
            {
                inplacetranspose(z, 1, n, 1, n, wtemp);
            }
            return result;
        }
    }
    
    //
    // Order eigenvalues and eigenvectors.
    //
    if( zneeded==0 )
    {
        
        //
        // Sort
        //
        if( n==1 )
        {
            return result;
        }
        if( n==2 )
        {
            if( d(1)>d(2) )
            {
                tmp = d(1);
                d(1) = d(2);
                d(2) = tmp;
            }
            return result;
        }
        i = 2;
        do
        {
            t = i;
            while(t!=1)
            {
                k = t/2;
                if( d(k)>=d(t) )
                {
                    t = 1;
                }
                else
                {
                    tmp = d(k);
                    d(k) = d(t);
                    d(t) = tmp;
                    t = k;
                }
            }
            i = i+1;
        }
        while(i<=n);
        i = n-1;
        do
        {
            tmp = d(i+1);
            d(i+1) = d(1);
            d(+1) = tmp;
            t = 1;
            while(t!=0)
            {
                k = 2*t;
                if( k>i )
                {
                    t = 0;
                }
                else
                {
                    if( k<i )
                    {
                        if( d(k+1)>d(k) )
                        {
                            k = k+1;
                        }
                    }
                    if( d(t)>=d(k) )
                    {
                        t = 0;
                    }
                    else
                    {
                        tmp = d(k);
                        d(k) = d(t);
                        d(t) = tmp;
                        t = k;
                    }
                }
            }
            i = i-1;
        }
        while(i>=1);
    }
    else
    {
        
        //
        // Use Selection Sort to minimize swaps of eigenvectors
        //
        for(ii = 2; ii <= n; ii++)
        {
            i = ii-1;
            k = i;
            p = d(i);
            for(j = ii; j <= n; j++)
            {
                if( d(j)<p )
                {
                    k = j;
                    p = d(j);
                }
            }
            if( k!=i )
            {
                d(k) = d(i);
                d(i) = p;
                if( wastranspose )
                {
                    ap::vmove(&wtemp(1), &z(i, 1), ap::vlen(1,n));
                    ap::vmove(&z(i, 1), &z(k, 1), ap::vlen(1,n));
                    ap::vmove(&z(k, 1), &wtemp(1), ap::vlen(1,n));
                }
                else
                {
                    ap::vmove(wtemp.getvector(1, zrows), z.getcolumn(i, 1, zrows));
                    ap::vmove(z.getcolumn(i, 1, zrows), z.getcolumn(k, 1, zrows));
                    ap::vmove(z.getcolumn(k, 1, zrows), wtemp.getvector(1, zrows));
                }
            }
        }
        if( wastranspose )
        {
            inplacetranspose(z, 1, n, 1, n, wtemp);
        }
    }
    return result;
}


/*************************************************************************
DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix

   [  A   B  ]
   [  B   C  ].

On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
eigenvector for RT1, giving the decomposition

   [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
   [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].


  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************/
template<class T>
void tdevdev2(const T& a,
     const T& b,
     const T& c,
     T& rt1,
     T& rt2,
     T& cs1,
     T& sn1)
{
    int sgn1;
    int sgn2;
    T ab;
    T acmn;
    T acmx;
    T acs;
    T adf;
    T cs;
    T ct;
    T df;
    T rt;
    T sm;
    T tb;
    T tn;

    
    //
    // Compute the eigenvalues
    //
    sm = a+c;
    df = a-c;
    adf = fabs(df);
    tb = b+b;
    ab = fabs(tb);
    if( fabs(a)>fabs(c) )
    {
        acmx = a;
        acmn = c;
    }
    else
    {
        acmx = c;
        acmn = a;
    }
    if( adf>ab )
    {
        rt = adf*sqrt(1+ap::sqr(ab/adf));
    }
    else
    {
        if( adf<ab )
        {
            rt = ab*sqrt(1+ap::sqr(adf/ab));
        }
        else
        {
            
            //
            // Includes case AB=ADF=0
            //
            rt = ab*sqrt(T(2));
        }
    }
    if( sm<0 )
    {
        rt1 = 0.5*(sm-rt);
        sgn1 = -1;
        
        //
        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
        //
        rt2 = acmx/rt1*acmn-b/rt1*b;
    }
    else
    {
        if( sm>0 )
        {
            rt1 = 0.5*(sm+rt);
            sgn1 = 1;
            
            //
            // Order of execution important.
            // To get fully accurate smaller eigenvalue,
            // next line needs to be executed in higher precision.
            //
            rt2 = acmx/rt1*acmn-b/rt1*b;
        }
        else
        {
            
            //
            // Includes case RT1 = RT2 = 0
            //
            rt1 = 0.5*rt;
            rt2 = -0.5*rt;
            sgn1 = 1;
        }
    }
    
    //
    // Compute the eigenvector
    //
    if( df>=0 )
    {
        cs = df+rt;
        sgn2 = 1;
    }
    else
    {
        cs = df-rt;
        sgn2 = -1;
    }
    acs = fabs(cs);
    if( acs>ab )
    {
        ct = -tb/cs;
        sn1 = 1/sqrt(1+ct*ct);
        cs1 = ct*sn1;
    }
    else
    {
        if( ab==0 )
        {
            cs1 = 1;
            sn1 = T(0.0);
        }
        else
        {
            tn = -cs/tb;
            cs1 = 1/sqrt(1+tn*tn);
            sn1 = tn*cs1;
        }
    }
    if( sgn1==sgn2 )
    {
        tn = cs1;
        cs1 = -sn1;
        sn1 = tn;
    }
}


/*************************************************************************
Internal routine
*************************************************************************/
template<class T>
T tdevdpythag(T a, T b)
{
    T result;

    if( fabs(a)<fabs(b) )
    {
        result = fabs(b)*sqrt(1+ap::sqr(a/b));
    }
    else
    {
        result = fabs(a)*sqrt(1+ap::sqr(b/a));
    }
    return result;
}


/*************************************************************************
Internal routine
*************************************************************************/
template<class T>
T tdevdextsign(T a, T b)
{
    T result;

    if( b>=0 )
    {
        result = fabs(a);
    }
    else
    {
        result = -fabs(a);
    }
    return result;
}


#endif