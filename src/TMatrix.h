#ifndef GUARD_TMatrix_h
#define GUARD_TMatrix_h

//-------------------------------------------------------------------------------------------
//
// TMatrix.h
//
// Author: Struan Robertson
// Date:   30/Mar/2003
//
// This header file contains the declaration of the TMatrix class.  This class inherits from
// Matrix and wraps calls to LAPACK functions.
//
//-------------------------------------------------------------------------------------------
#include "Matrix.h"
#include <string>
#include <cmath>
#include <climits>
#include <stdio.h>
#include <vector>

namespace mesmer
{
  template<class T>
	class TMatrix : public Matrix<T> {

	public:

		// Constructor
		TMatrix(int n) : Matrix<T>(n, 0.0) { } ;

		//
		// Wrapped call to EISPACK routine to diagonalise matrix.
		//
		void diagonalize(T *rr) {

			int size ;
			size = static_cast<int>(this->size()) ;

			//  Allocate memory for work array
			T *work = new T[size] ;
			T *rrProxy = new T[size] ;

			tred2(this->m_matrix, size, rrProxy, work) ;
			tqli(rrProxy, work, size, this->m_matrix) ;

			for (int i = 0; i < size; ++i){
				rr[i] = to_double(rrProxy[i]);
			}

			delete [] work ;
			delete [] rrProxy ;

		}

		// 
		// Solve a set of linear equations with a single right hand side.
		//
		void solveLinearEquationSet(T *rr) {

			int size ;
			size = static_cast<int>(this->size()) ;

			//  Allocate memory for work array
			int *indx = new int[size] ;

			T d ;
			ludcmp(this->m_matrix, size, indx, d) ;
			lubksb(this->m_matrix, size, indx, rr) ;

      delete [] indx ;

    };

    // Matrix inversion method by Gaussian elimination
    int invert();

  private:

		//
		// EISPACK methods for diagonalizing matrix.
		//
		void    tred2   (T **a, int n, T *d, T *e) ;
		void    tqli    (T *d, T *e, int n, T **z) ;
		T  pythag  (T a, T b) ;

		//
		// NR LU methods for linear equation solving.
		//
		void ludcmp(T **a,  int n, int *indx, T d) ;
		void lubksb(T **a,  int n, int *indx, T* b) ;

	} ;

//-------------------------------------------------------------------------------------------
// EISPACK tred2 function.
//
// Householder reduction of matrix a to tridiagonal form.
//   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
//   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
//   Springer-Verlag, 1976, pp. 489-494.
//   W H Press et al., Numerical Recipes in C, Cambridge U P,
//   1988, pp. 373-374.
//-------------------------------------------------------------------------------------------

    template<class T>
    void TMatrix<T>::tred2(T **a, int n, T *d, T *e)
    {
        int l, k, j, i;
        T scale, hh, h, g, f;

        for (i=n;i>=2;--i) {
            l=i-1;
            h=scale=0.0;
            if (l > 1) {

                for (k=1;k<=l;++k)
                    scale += fabs(a[i-1][k-1]);
                if (scale == 0.0)
                    e[i-1]=a[i-1][l-1];
                else {
                    for (k=1;k<=l;++k) {
                        a[i-1][k-1] /= scale;
                        h += a[i-1][k-1]*a[i-1][k-1];
                    }
                    f=a[i-1][l-1];
                    g = f > 0.0 ? -sqrt(h) : sqrt(h);
                    e[i-1]=scale*g;
                    h -= f*g;
                    a[i-1][l-1]=f-g;
                    f=0.0;
                    for (j=1;j<=l;++j) {
                        /* Next statement can be omitted if eigenvectors not wanted */
                        a[j-1][i-1]=a[i-1][j-1]/h;
                        g=0.0;

                        for (k=1;k<=j;++k)
                            g += a[j-1][k-1]*a[i-1][k-1];

                        for (k=j+1;k<=l;++k)
                            g += a[k-1][j-1]*a[i-1][k-1];
                        e[j-1]=g/h;
                        f += e[j-1]*a[i-1][j-1];
                    }
                    hh=f/(h+h);
                    for (j=1;j<=l;++j) {
                        f=a[i-1][j-1];
                        e[j-1]=g=e[j-1]-hh*f;

                        for (k=1;k<=j;++k)
                            a[j-1][k-1] -= (f*e[k-1]+g*a[i-1][k-1]);
                    }
                }
            } else
                e[i-1]=a[i-1][l-1];
            d[i-1]=h;
        }
        /* Next statement can be omitted if eigenvectors not wanted */
        d[0]=0.0;
        e[0]=0.0;
        /* Contents of this loop can be omitted if eigenvectors not
        wanted except for statement d[i-1]=a[i-1][i-1]; */
        for (i=1;i<=n;++i) {
            l=i-1;
            if (d[i-1] != .0) {
                for (j=1;j<=l;++j) {
                    g=0.0;

                    for (k=1;k<=l;++k)
                        g += a[i-1][k-1]*a[k-1][j-1];

                    for (k=1;k<=l;++k)
                        a[k-1][j-1] -= g*a[k-1][i-1];
                }
            }
            d[i-1]=a[i-1][i-1];
            a[i-1][i-1]=1.0;

            for (j=1;j<=l;++j) a[j-1][i-1]=a[i-1][j-1]=0.0;
        }
    }

    //-------------------------------------------------------------------------------------------
    // EISPACK tqli function.
    //
    // QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
    // symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2.
    //
    // On input:
    //    d[1..n] contains the diagonal elements of the tridiagonal matrix.
    //    e[1..n] contains the subdiagonal elements of the tridiagonal matrix.
    // with e[1] arbitrary.
    // On output:
    //    d[1..n] returns the eigenvalues.
    //    e[1..n] is destroyed.
    //
    // When finding only the eigenvalues, several lines may be omitted, as noted in the comments.
    //
    // If the eigenvectors of a tridiagonal matrix are desired, the matrix z[1..n][1..n] is input
    // as the identity matrix. If the eigenvectors of a matrix that has been reduced by tred2 are
    // required, then z is input as the matrix output by tred2. In either case, the kth column of
    // z returns the normalized eigenvector corresponding to d[k].
    //
    //-------------------------------------------------------------------------------------------
    template<class T>
    void TMatrix<T>::tqli(T *d, T *e, int n, T **z)
    {
        int m,l,iter,i,k;
        T s,r,p,g,f,dd,c,b;

        for (i=2;i<=n;++i) e[i-2]=e[i-1];
        e[n-1]=0.0;
        for (l=1;l<=n;++l) {
            iter=0;
            do {
                for (m=l;m<=n-1;++m) {
                    dd=fabs(d[m-1])+fabs(d[m]);
                    if (fabs(e[m-1])+dd == dd) break;
                }
                if (m != l) {
                    // if (iter++ == 30) fprintf(stderr, "Too many iterations in TQLI");
                    if (iter++ == 60) fprintf(stderr, "Too many iterations in TQLI");
                    /* CHL
                    Source: http://www.nr.com/forum/showthread.php?t=592
                    I hope that bellow words will be useful for you.
                    See thread under the title: Possible convergence problems in svdcmp, jacobi, tqli, hqr by Saul Teukolsky
                    in Forum: Official Bug Reports with known bugs. May be this is a reason of slow convergency.
                    It is good check, that matrix is symmetric and to work with T accuracy. I have known also versions
                    with increased number of iterations (200 for example). But I think that this experimental number is right
                    in any case: if you have not convergency for 30 iterations, there is no convergency at all.
                    SVD method used in book is an intrinsic iterative procedure, 30 iterations is a good number to
                    convergency up to numerical accuracy. Evgeny
                    */
                    g=(d[l]-d[l-1])/(2.0*e[l-1]);
                    r=sqrt((g*g)+1.0);
                    g=d[m-1]-d[l-1]+e[l-1]/(g + (g < 0.0 ? -fabs(r) : fabs(r)));
                    s=c=1.0;
                    p=0.0;
                    for (i=m-1;i>=l;--i) {
                        f=s*e[i-1];
                        b=c*e[i-1];
                        if (fabs(f) >= fabs(g)) {
                            c=g/f;
                            r=sqrt((c*c)+1.0);
                            e[i]=f*r;
                            c *= (s=1.0/r);
                        } else {
                            s=f/g;
                            r=sqrt((s*s)+1.0);
                            e[i]=g*r;
                            s *= (c=1.0/r);
                        }
                        g=d[i]-p;
                        r=(d[i-1]-g)*s+2.0*c*b;
                        p=s*r;
                        d[i]=g+p;
                        g=c*r-b;
                        /* Next loop can be omitted if eigenvectors not wanted */

                        for (k=1;k<=n;++k) {
                            f=z[k-1][i];
                            z[k-1][i]=s*z[k-1][i-1]+c*f;
                            z[k-1][i-1]=c*z[k-1][i-1]-s*f;
                        }
                    }
                    d[l-1]=d[l-1]-p;
                    e[l-1]=g;
                    e[m-1]=0.0;
                }
            } while (m != l);
        }

        // Order eigenvalues and eigenvectors.

        for (int ii = 1; ii < n; ++ii) {
            i = ii - 1;
            k = i;
            p = d[i];
            for (int j = ii; j < n; ++j) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }
            if (k!=i) {
                d[k] = d[i];
                d[i] = p;
                for (int j = 0; j < n; ++j) {
                    p = z[j][i];
                    z[j][i] = z[j][k];
                    z[j][k] = p;
                }
            }
        }

    }

    //-------------------------------------------------------------------------------------------
    // EISPACK pythag function.
    //
    // Finds sqrtl(a**2+b**2) without overflow or destructive underflow.
    //
    //-------------------------------------------------------------------------------------------
    template<class T>
    T TMatrix<T>::pythag(T a, T b)
    {
        T p,r,s,t,u;
        if (fabs(a) > fabs(b)) p = fabs(a);
        else p = fabs(b);
        if (p == 0.) goto label_1;
        if (fabs(a) > fabs(b)) r = fabs(b);
        else r = fabs(a);
        r = (r/p)*(r/p);

label_2: t = 4. + r;
        if (t == 4.) goto label_1;
        s = r/t;
        u = 1. + 2. *s;
        p *= u;
        r = ((s/u)*(s/u))*r;

        goto label_2;
label_1: return(p);
    }

    //
    // NR LU methods for linear equation solving.
    //
    template<class T>
    void TMatrix<T>::ludcmp(T **a,  int n, int *indx, T d) {

        const T TINY = 1.0e-20 ;
        int i, imax, j, k ;
        T big, dum, sum, temp ;

        T *work = new T[n] ;
        d = 1.0 ;
        for (i = 1; i < n ; ++i) {
            big = 0.0 ;
            for (j = 0; j < n ; ++j) {
                if ((temp = fabs(a[i][j])) > big)
                    big = temp ;
            }
            if (big == 0.0) {
                fprintf(stderr, "Singular Matrix in routine ludcmp");
                exit(0) ;
            }
            work[i] = 1.0/big ;
        }

        for (j = 0; j < n ; ++j) {
            for (i = 0; i < j ; ++i) {
                sum = a[i][j] ;
                for (k = 0; k < i ; ++k)
                    sum -= a[i][k]*a[k][j] ;

                a[i][j] = sum ;
            }
            big = 0.0 ;
            for (i = j; i < n; ++i) {
                sum = a[i][j] ;
                for (k = 0; k < j ; ++k)
                    sum -= a[i][k]*a[k][j] ;

                a[i][j] = sum ;

                if ( (dum = work[i]*fabs(sum)) >= big) {
                    big = dum ;
                    imax = i ;
                }
            }
            if (j != imax) {
                for (k = 0; k < n; ++k) {
                    dum = a[imax][k] ;
                    a[imax][k] = a[j][k] ;
                    a[j][k] = dum ;
                }
                d = -d ;
                work[imax] = work[j] ;
            }
            indx[j] = imax ;
            if (a[i][j] == 0.0)
                a[j][i] = TINY ;

            if (j != n-1) {
                dum = 1.0/(a[j][i]) ;
                for (i = j+1; i < n; ++i)
                    a[i][j] = dum ;
            }

        }

        delete [] work ;

    }

    template<class T>
    void TMatrix<T>::lubksb(T **a,  int n, int *indx, T* b) {

        int i, ii = 0, ip, j ;
        T sum ;

        for (i = 0; i < n; ++i) {
            ip = indx[i] ;
            sum = b[ip] ;
            b[ip] = b[i] ;
            if (ii != 0) {
                for (j = ii-1; j < i; ++j)
                    sum -= a[i][j]*b[j] ;
            } else {
                if (sum != 0.0)
                    ii = i+1 ;
                b[i] = sum ;
            }
        }
        for (i = n-1; i >= 0; i--) {
            sum = b[i] ;
            for (j=i+1; j < n; ++j)
                sum -= a[i][j]*b[j] ;
            b[i] = sum/a[i][j] ;
        }
    }

  template<class T>
  int TMatrix<T>::invert(){
  /*######################################################################
    Author: Chi-Hsiu Liang
    Matrix  inversion, real symmetric a of order n.
    Gaussian elimination with interchanges.
    Sometimes this routine can cause problem, there is no disgnostic functions
    inside this routine. Users have to find it out by themselves
    To see what's the input matrix one should uncomments the section which
    can print out the matrix.
    (Useful checklist of whether a matrix has an inverse or not is below.
    A. Matrix with two columns/rows completely the same does not have an inverse
    B. Matrix with determinant zero does not have an inverse, of course this include
    the one above.)
    USAGE: @s_inv = invert(\@s_inv, orbNo);
    where
    @s_inv :a two dimensional matrix to be inverted.
    orbNo :member of the indicated matrix.
    ######################################################################*/
	  int n = static_cast<int>(this->size()) ;
    if (n > INT_MAX) return 2; // The matrix size is too large to process.
    T divide, ratio;
    std::vector<int> zeroCount(n, 0);

    //-------------------------------
    //produce a unit vector of size n, and copy the incoming matrix
    TMatrix<T> m2(n);
    TMatrix<T> m1(n);
    for (int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        m1[i][j] = this->m_matrix[i][j];
        m2[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }
    //-------------------------------

    for(int j = 0; j < n; ++j){
      for (int i = 0; i < n; ++i){
        if (m1[i][j] == 0.0){
          zeroCount[j]++; if (zeroCount[j] == n) return 1;
          /* If there is a zero column the whole array has no inverse.*/
        }
        else{
          if (i < j){
            if (m1[j][j] == 0.0){ /* Add the former row to the j'th row if the main row is empty.*/
              for (int col = 0; col < n; ++col){
                m1[i][col] = m1[j][col];
                m2[i][col] = m2[j][col];
              }
            }
          }
          else if (i > j){ // Add the later row to the j'th row.
            if (zeroCount[j] == i){
              for (int col = j; col < 3; ++col) swap(m1[i][col], m1[j][col]);
              for (int col = 0; col < 3; ++col) swap(m2[i][col], m2[j][col]);
              i = j - 1; zeroCount[j] = 0;
            }
            else{
              for (int col = 0; col < n; ++col){
                m1[i][col] = m1[j][col];
                m2[i][col] = m2[j][col];
              }
            }
          }
          //i = j;
          else{ // in this case i = j
            if (m1[i][j] != 1.0){
              divide = m1[i][j];
              for (int col = i; col < n; ++col) m1[i][col] /= divide; // normalise i'th row.
              for (int col = 0; col < n; ++col) m2[i][col] /= divide;
            }
            for (int row = 0; row < n; ++row){
              if (row == i) continue;
              ratio = m1[row][j] / m1[i][j];
              for (int col = i; col < n; ++col) m1[row][col] -= m1[i][col] * ratio;
                /* Only alterations after the i'th indice are necessary*/
              for (int col = 0; col < n; ++col) m2[row][col] -= m2[i][col] * ratio;
            }
          }
        }
      }
    }

    for(int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        this->m_matrix[i][j] = m2[i][j];
      }
    }
    return 0;
  }
}//namespacer mesmer


#endif // GUARD_TMatrix_h
