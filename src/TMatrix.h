#ifndef GUARD_TMatrix_h
#define GUARD_TMatrix_h

//-------------------------------------------------------------------------------------------
//
// TMatrix.h
//
// Author: Struan Robertson
// Date:   30/Mar/2003
//
// This header file contains the declaration of the TMatrix class. This class inherits from
// Matrix and wraps calls to EISPACK functions.
//
//-------------------------------------------------------------------------------------------
#include "Matrix.h"
#include <string>
#include <cstring>
#include <cmath>
#include <climits>
#include <stdio.h>
#include <vector>
#include <complex>
#include <stdexcept>
#include <stdlib.h>
#include "Persistence.h"
#include "qd/dd_real.h"
#include "qd/qd_real.h"

namespace mesmer
{
  template<class T>
  class TMatrix : public MeMatrix<T> {

  public:

    // Constructors
    TMatrix( size_t n, const T& init = T(0.0)) : MeMatrix<T>(n, init) { } ;

    // Copy constructor
    TMatrix(const MeMatrix<T>& rhs ) : MeMatrix<T>(rhs) { } ;

    //
    // Wrapped calls to EISPACK routines to diagonalise matrix.
    //
    void diagonalize(T *rr) {

      size_t size = this->size() ;

      //  Allocate memory for work array
      T *work = new T[size] ;
      T *rrProxy = new T[size] ;

      // vector<size_t> index(size, 0) ;
      // permuteMatrix(this->m_matrix, index) ;
      tred2(this->m_matrix, size, rrProxy, work) ;
      tqli(rrProxy, work, size, this->m_matrix) ;
      // unPermuteEigenEigenvectors(this->m_matrix, index) ;

      for (size_t i = 0; i < size; ++i){
        rr[i] = rrProxy[i];
      }

      bool diagnostic(false) ;
      if (diagnostic) {
        ctest << endl ;
        ctest << "Largest eigenvalue:  " << rr[0] << endl ;
        ctest << "Smallest eigenvalue: " << rr[size-2] << endl ;
        ctest << "Ratio:               " << rr[0]/rr[size-2] << endl << endl ;
      }

      delete [] work ;
      delete [] rrProxy ;

    }

    //
    // Solve a set of linear equations with a single right hand side.
    //
    void solveLinearEquationSet(T *rr) {

      const size_t size = this->size() ;

      //  Allocate memory for work array
      vector<size_t> indx(size, 0) ;

      if (ludcmp(this->m_matrix, size, indx)){
        throw(std::runtime_error("Error encountered in LU decompostion invoked from solveLinearEquationSet."));
      }

      lubksb(this->m_matrix, size, indx, rr) ;

    };

	// Determinant of Matrix.
	T Determinant() {

      const size_t size = this->size() ;

      //  Allocate memory for work array
      TMatrix<T> work = *this ;
	  
	  //  Allocate memory for work array
      vector<size_t> indx(size, 0) ;

      if (ludcmp(work.m_matrix, size, indx)) {
        throw(std::runtime_error("Error encountered in LU decompostion invoked from Determinant."));
      }
	  T det(T(1.0)) ;
      for (size_t i = 0; i < size ; i++) {
		det *= work[i][i] ;
	  }

	  return det ;
	}

    // Code adapted from C acording to the algorithm given at rosettacode.org/wiki/Cholesky_decomposition.
    void cholesky() {

      size_t size = this->size() ;

      //  Allocate memory for work array
      TMatrix<T> work(size,0.0) ;

      for (size_t i = 0; i < size; i++){
        for (size_t j = 0; j < (i+1); j++){
          T sum = T(0.0) ;
          for (size_t k = 0; k < j; k++){
            sum += work[i][k]*work[j][k];
          }
          if (i==j){
            work[i][j] = sqrt((*this)[i][j]-sum);
          }else{
            work[i][j] = (1.0/work[j][j])*((*this)[i][j]-sum);
          }
        }
      }
      // Copy work matrix back to orignal matrix.
      *this = work;
    };

    // Matrix inversion method by LU decomposition
    int invertLUdecomposition();

    void normalizeProbabilityMatrix();

    // Print out the contents of the matrix.
    void print(std::string& title, std::ostream& output_stream, int nr = -1, int nc = -1, int fr = -1, int fc = -1) const ;

    // Apply the Gram-Schmidt procedure to orthogonalize the current matrix.
    void GramSchimdt(size_t root_vector ) ;

    // Transpose matrix.
    void Transpose() ;

    // Write matrix to an XML stream.
    void WriteToXML(PersistPtr pp) ;

    void showFinalBits(const size_t n, bool isTabbed = false);

    //
    // EISPACK methods for diagonalizing matrix.
    //
    static void tred2 (T **a, size_t n, T *d, T *e) ;
    static void tqli  (T *d, T *e, size_t n, T **z) ;
    static void tqlev (T *d, T *e, size_t n) ;
    static T    pythag(T a, T b) ;

  private:

    //
    // The following two methods permute a matrix, in order to reduce the effects
    // of numerical rounding in the Househlider method (see numerical recipes) and
    // unpermute the associated eignvector matrix. SHR 7/Mar/2011: this experiment
    // appeared to have only a minor effect on numerical values, hence the statements
    // in diagonalize() have been commented out as the benefits do not appear to out
    // way the cost in CPU time at present.
    //
    void permuteMatrix(T **a, vector<size_t>& index) ;
    void unPermuteEigenEigenvectors(T **a, vector<size_t>& index) ;

    //
    // NR LU methods for linear equation solving.
    //
    int  ludcmp(T **a, const size_t &n, vector<size_t> &indx) ;
    void lubksb(T **a, const size_t &n, vector<size_t> &indx, T* b) ;

  } ;

  // Matrix mutiplication operator.
  template<class T>
  TMatrix<T> operator*(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {

    size_t msize = lhs.size() ;
    if (rhs.size() != msize) {
      // Throw error.
    }

    TMatrix<T> result(msize) ;

    for (size_t i(0) ; i < msize ; i++) {
      for (size_t j(0) ; j < msize ; j++) {
        T sm(0.0) ;
        for (size_t k(0) ; k < msize ; k++) {
          sm += lhs[i][k]*rhs[k][j] ;
        }
        result[i][j] = sm ;
      }
    }  

    return result ; // Note result goes via the stack!

  }

  // Matrix vector mutiplication operator.
  template<class T>
  void operator*=(vector<T>& rhs, const TMatrix<T>& lhs) {

    size_t msize = lhs.size() ;
    if (rhs.size() != msize) {
      // Throw error.
    }

    vector<T> result(msize) ;

    for (size_t j(0) ; j < msize ; j++) {
      T sm(0.0) ;
      for (size_t k(0) ; k < msize ; k++) {
        sm += lhs[j][k]*rhs[k] ;
      }
      result[j] = sm ;
    }

    rhs = result ; 

  }


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
  void TMatrix<T>::tred2(T **a, size_t n, T *d, T *e)
  {
    size_t l, k, j, i;
    T scale, hh, h, g, f;

    for (i=n-1; i>0; i--) {
      l=i-1 ;
      h=scale=T(0.0) ;
      if (l > 0) {

        for (k=0; k<l+1; k++)
          scale += fabs(a[i][k]);
        if (scale == T(0.0))
          e[i]=a[i][l];
        else {
          for (k=0; k<l+1; k++) {
            a[i][k] /= scale;
            h += a[i][k]*a[i][k];
          }
          f=a[i][l];
          g = (f > T(0.0)) ? -sqrt(h) : sqrt(h);
          e[i]=scale*g;
          h -= f*g;
          a[i][l]=f-g;
          f=T(0.0);
          for (j=0; j<l+1; j++) {
            // Next statement can be omitted if eigenvectors not wanted.
            a[j][i]=a[i][j]/h;
            g=T(0.0);

            for (k=0; k<j+1; k++)
              g += a[j][k]*a[i][k];

            for (k=j+1; k<l+1; k++)
              g += a[k][j]*a[i][k];

            e[j]=g/h;
            f += e[j]*a[i][j];
          }
          hh=f/(h+h);
          for (j=0; j<l+1; j++) {
            f=a[i][j];
            e[j]=g=e[j]-hh*f;

            for (k=0; k<j+1; k++)
              a[j][k] -= (f*e[k]+g*a[i][k]);
          }
        }
      } else
        e[i]=a[i][l];

      d[i] = h;
    }
    // Next statement can be omitted if eigenvectors not wanted.
    d[0]=0.0;
    e[0]=0.0;
    // Contents of this loop can be omitted if eigenvectors not
    // wanted except for statement "d[i]=a[i][i];".
    for (i=0; i<n; i++) {
      l=i;
      if (d[i] != 0.0) {
        for (j=0; j<=l; j++) {
          g=0.0;

          for (k=0; k<l; k++)
            g += a[i][k]*a[k][j];

          for (k=0; k<l; k++)
            a[k][j] -= g*a[k][i];
        }
      }
      d[i] = a[i][i];
      a[i][i] = 1.0;

      for (j=0; j<l; j++) 
        a[j][i] = a[i][j] = 0.0;
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
  void TMatrix<T>::tqli(T *d, T *e, size_t n, T **z)
  {
    size_t m,l,iter,i,k;
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
          if (iter++ == 60) { 
            fprintf(stderr, "Too many iterations in TQLI");
            throw(std::runtime_error("Too many iterations in TQLI."));
          }
          /* CHL
          Source: http://www.nr.com/forum/showthread.php?t=592
          I hope that bellow words will be useful for you.
          See thread under the title: Possible convergence problems in svdcmp, jacobi, tqli, hqr by Saul Teukolsky
          in Forum: Official Bug Reports with known bugs. May be this is a reason of slow convergency.
          It is good check, that matrix is symmetric and to work with double accuracy. I have known also versions
          with increased number of iterations (200 for example). But I think that this experimental number is right
          in any case: if you have not convergency for 30 iterations, there is no convergency at all.
          SVD method used in book is an intrinsic iterative procedure, 30 iterations is a good number to
          convergency up to numerical accuracy. Evgeny
          */
          g=(d[l]-d[l-1])/(2.0*e[l-1]);
          r=sqrt((g*g)+1.0);
          //r = pythag(g, 1.0) ;
          g=d[m-1]-d[l-1]+e[l-1]/(g + (g < 0.0 ? -fabs(r) : fabs(r)));
          s=c=T(1.0) ;
          p=0.0;
          for (i=m-1;i>=l;--i) {
            f=s*e[i-1];
            b=c*e[i-1];
            if (fabs(f) >= fabs(g)) {
              c=g/f;
              r=sqrt((c*c)+1.0);
              //r = pythag(c, 1.0) ;
              e[i]=f*r;
              c *= (s=1.0/r);
            } else {
              s=f/g;
              r=sqrt((s*s)+1.0);
              //r = pythag(s, 1.0) ;
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

    for (size_t ii = 1; ii < n; ++ii) {
      i = ii - 1;
      k = i;
      p = d[i];
      for (size_t j = ii; j < n; ++j) {
        if (d[j] < p) {
          k = j;
          p = d[j];
        }
      }
      if (k!=i) {
        d[k] = d[i];
        d[i] = p;
        for (size_t j = 0; j < n; ++j) {
          p = z[j][i];
          z[j][i] = z[j][k];
          z[j][k] = p;
        }
      }
    }

  }


  //-------------------------------------------------------------------------------------------
  // Function tqlev - eigenvalue only method based on tqli.
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
  //-------------------------------------------------------------------------------------------
  template<class T>
  void TMatrix<T>::tqlev(T *d, T *e, size_t n)
  {
    size_t m,l,iter,i,k;
    T s,r,p,g,f,dd,c,b;

    if (n==0) return ;

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
          // if (iter++ == 30) fprintf(stderr, "Too many iterations in tqlev");
          if (iter++ == 60) {
            fprintf(stderr, "Too many iterations in tqlev");
            throw(std::runtime_error("Too many iterations in tqlev.")); ;
          }
          /* CHL
          Source: http://www.nr.com/forum/showthread.php?t=592
          I hope that bellow words will be useful for you.
          See thread under the title: Possible convergence problems in svdcmp, jacobi, tqli, hqr by Saul Teukolsky
          in Forum: Official Bug Reports with known bugs. May be this is a reason of slow convergency.
          It is good check, that matrix is symmetric and to work with double accuracy. I have known also versions
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
          }
          d[l-1]=d[l-1]-p;
          e[l-1]=g;
          e[m-1]=0.0;
        }
      } while (m != l);
    }

    // Order eigenvalues.

    for (size_t ii = 1; ii < n; ++ii) {
      i = ii - 1;
      k = i;
      p = d[i];
      for (size_t j = ii; j < n; ++j) {
        if (d[j] < p) {
          k = j;
          p = d[j];
        }
      }
      if (k!=i) {
        d[k] = d[i];
        d[i] = p;
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

    T absa = fabs(a) ;
    T absb = fabs(b) ;

    p = (absa > absb) ? absa : absb ;
    if (p == 0.) return p ;
    r = (absa > absb) ? absb : absa ;
    r = (r/p)*(r/p);

    while ((t = T(4.0) + r) != T(4.0)) {
      s  = r/t;
      u  = T(1.0) + T(2.0)*s;
      p *= u;
      r *= ((s/u)*(s/u)) ;
    }

    return p ;
  }

  //
  // This method permutes the matrix, in order to reduce the effects
  // of numerical rounding in the Householder method (see numerical recipes).
  //
  template<class T>
  void TMatrix<T>::permuteMatrix(T **a, vector<size_t>& index) {

    size_t size = this->size() ;
    vector<T> diagonalElements(size,0.0) ;

    for (size_t i(0) ; i < size ; i++) {
      diagonalElements[i] = fabs(a[i][i]) ;
      index[i] = i ;
    }

    // Order matrix columns based on the magnitude of the diagonal elements.

    for (size_t ii(1); ii < size; ++ii) {
      size_t i = ii - 1;
      size_t k = i;
      T p = diagonalElements[i];
      for (size_t j(ii); j < size ; ++j) {
        if (diagonalElements[j] < p) {
          k = j;
          p = diagonalElements[j];
        }
      }
      if (k != i) {
        diagonalElements[k] = diagonalElements[i];
        diagonalElements[i] = p;
        swap(index[i],index[k]) ;

        // Swap columns.
        for (size_t j(0); j < size; ++j) {
          swap(a[j][i], a[j][k]) ;
        }

        // Swap rows.
        for (size_t j(0); j < size; ++j) {
          swap(a[i][j], a[k][j]) ;
        }
      }
    }

  }

  //
  // This method unpermutes the eignvector matrix rows following diagonalization.
  //
  template<class T>
  void TMatrix<T>::unPermuteEigenEigenvectors(T **a, vector<size_t>& index) {

    size_t size = this->size() ;

    vector<size_t> invIndex(size, 0) ;
    for (size_t i(0) ; i < size ; i++) {
      invIndex[index[i]] = i ;
    }

    for (size_t i(0) ; i < size ; i++) {
      size_t k = invIndex[i] ; 
      if (k != i) {
        invIndex[index[i]] = k ;
        invIndex[i] = i ; 
        swap(index[i],index[k]) ;
        for (size_t j(0) ; j < size ; j++) {
          swap(a[i][j],a[k][j]) ;
        }
      }
    }

  }

  //
  // NR LU methods for linear equation solving.
  //
  /**************************************************************
  * Given an N x N matrix A, this routine replaces it by the LU *
  * decomposition of a rowwise permutation of itself. A and N   *
  * are input. INDX is an output vector which records the row   *
  * permutation effected by the partial pivoting; D is output   *
  * as -1 or 1, depending on whether the number of row inter-   *
  * changes was even or odd, respectively. This routine is used *
  * in combination with LUBKSB to solve linear equations or to  *
  * invert a matrix. Return code is 1, if matrix is singular.   *
  **************************************************************/
  template<class T>
  int TMatrix<T>::ludcmp(T **a,  const size_t &n, vector<size_t> &indx) {

    size_t imax(0);

    T big, dum, sum, temp ;
    T tiny = numeric_limits<T>::epsilon();

    vector<T> work(n,T(0.0)) ;

    for (size_t i(0); i < n ; ++i) {
      big = T(0.0) ;
      for (size_t j(0); j < n ; ++j) {
        if ((temp = fabs(a[i][j])) > big){
          big = temp ;
        }
      }
      if (big == T(0.0)) {
        cerr << "Singular Matrix in routine ludcmp";
        return 1;
      }
      work[i] = T(1.0)/big ;
    }

    for (size_t j(0); j < n ; ++j) {
      for (size_t i(0); i < j ; ++i) {
        sum = a[i][j] ;
        for (size_t k(0); k < i ; ++k){
          sum -= a[i][k]*a[k][j] ;
        }
        a[i][j] = sum ;
      }
      big = T(0.0) ;
      for (size_t i(j); i < n; ++i) {
        sum = a[i][j] ;
        for (size_t k(0); k < j ; ++k)
          sum -= a[i][k]*a[k][j] ;

        a[i][j] = sum ;

        if ( (dum = work[i]*fabs(sum)) >= big) {
          big = dum ;
          imax = i ;
        }
      }
      if (j != imax) {
        for (size_t k(0); k < n; ++k) {
          dum = a[imax][k] ;
          a[imax][k] = a[j][k] ;
          a[j][k] = dum ;
        }
        work[imax] = work[j] ;
      }
      indx[j] = imax ;
      if (fabs(a[j][j]) < tiny){
        a[j][j] = tiny;
      }

      if (j != n-1) {
        dum = T(1.0)/(a[j][j]) ;
        for (size_t i(j+1); i < n; ++i)
          a[i][j] *= dum ;
      }

    }

    return 0;
  }

  /*****************************************************************
  * Solves the set of N linear equations A . X = B.  Here A is     *
  * input, not as the matrix A but rather as its LU decomposition, *
  * determined by the routine LUDCMP. INDX is input as the permuta-*
  * tion vector returned by LUDCMP. B is input as the right-hand   *
  * side vector B, and returns with the solution vector X. A, N and*
  * INDX are not modified by this routine and can be used for suc- *
  * cessive calls with different right-hand sides. This routine is *
  * also efficient for plain matrix inversion.                     *
  *****************************************************************/
  template<class T>
  void TMatrix<T>::lubksb(T **a,  const size_t &n,  vector<size_t> &indx, T* b) {

    size_t ii(0), ip;
    T sum ;

    for (size_t i(0); i < n; ++i) {
      ip = indx[i] ;
      sum = b[ip] ;
      b[ip] = b[i] ;
      if (ii >= 0) {
        for (size_t j(ii); j < i; ++j)
          sum -= a[i][j]*b[j] ;
      }
      else if (sum != T(0.0)){
        ii = i ;
      }
      b[i] = sum ;
    }
    for (int i(n-1); i >= 0; --i) {
      sum = b[i] ;
      if (size_t(i) < n-1){
        for (size_t j(i+1); j < n; ++j)
          sum -= a[i][j]*b[j] ;
      }
      b[i] = sum/a[i][i] ;
    }
  }

  template<class T>
  int TMatrix<T>::invertLUdecomposition(){

    const size_t size = this->size() ;

    // Allocate memory for work array.
    vector<size_t> indx(size, 0) ;

    // Construct an identity matrix as a primer for the inverse.
    MeMatrix<T> invM(size, T(0.0)) ; 
    for (size_t i(0); i < size; ++i){
      invM[i][i] = 1.0;
    }

    int rc = ludcmp(this->m_matrix, size, indx) ;

    // Call solver, if previous return code is ok,
    // to obtain inverse of A one column at a time.
    if (rc == 0) {
      T *temp = new T[size] ;
      for (size_t j(0); j < size; ++j) {
        for (size_t i(0); i < size; ++i) temp[i] = invM[i][j];
        lubksb(this->m_matrix, size, indx, temp);
        for (size_t i(0); i < size; ++i) invM[i][j] = temp[i];
      }
      for (size_t j(0); j < size; ++j) {
        for (size_t i(0); i < size; ++i){
          this->m_matrix[i][j] = invM[i][j];
        }
      }
      delete [] temp;
      return 0;
    }
    else {
      return 1;
    }

  }

  //
  // Normalize collision operator
  //
  template<class T>
  void TMatrix<T>::normalizeProbabilityMatrix(){

    //
    // Normalization of Probability matrix.
    // Normalising coefficients are found by using the fact that column sums
    // are unity. The procedure leads to a matrix that is of upper triangular
    // form and the normalisation constants are found by back substitution.
    //

    int i, j; //int makes sure the comparison to negative numbers meaningful (i >=0)

    int optrsize(int(this->size()));
    vector<T> work(optrsize) ;// Work space.

    T scaledRemain(0.0) ;
    for ( i = optrsize - 1 ; i >= 0 ; --i ) {

      T upperSum(0.0) ;
      for ( j = 0 ; j <= i ; ++j )
        upperSum += (*this)[j][i] ;

      if (upperSum > 0.0){
        if (i < optrsize - 1){
          scaledRemain = 0.0;
          for ( j = i + 1 ; j < optrsize ; ++j ){
            T scale = work[j];
            scaledRemain += (*this)[j][i] * scale ;
          }
        }
        work[i] = (1.0 - scaledRemain) / upperSum ;
      }
    }

    //
    // Apply normalization coefficients
    //
    for ( i = 0 ; i < optrsize ; ++i ) {
      (*this)[i][i] *= work[i] ;
      //T value = (*this)[i][i];
      for ( j = i + 1 ; j < optrsize ; ++j ) {
        (*this)[j][i] *= work[j] ;
        (*this)[i][j] *= work[j] ;
      }
    }

  }

  //
  // Print out the contents of the matrix.
  //
  template<class T>
  void TMatrix<T>::print(std::string& title, std::ostream& output_stream, int nr, int nc, int fr, int fc) const {

    size_t msize = this->size() ;
    size_t nrows = (nr < 0) ? msize : min(msize,size_t(nr)) ;
    size_t nclms = (nc < 0) ? msize : min(msize,size_t(nc)) ;
    size_t frow  = (fr < 0) ? 0     : min(msize,size_t(fr)) ;
    size_t fclm  = (fc < 0) ? 0     : min(msize,size_t(fc)) ;

    output_stream << endl << title << endl << "{" << endl ;
    for (size_t i(frow) ; i < nrows ; ++i) {
      for (size_t j(fclm) ; j < nclms ; ++j) {
        formatFloat(output_stream, (*this)[i][j],  6,  15) ;
		output_stream << "," ;
      }
      output_stream << endl ;
    }
    output_stream << "}" << endl ;

  }

  //
  // Apply the Gram-Schmidt procedure to orthogonalize the current matrix.
  //
  template<class T>
  void TMatrix<T>::GramSchimdt(size_t root_vector )  {

    size_t size = this->size() ;

    for (int i(size-1) ; i > -1 ; i--) { // Need to use int here as size_t is unsigned.

      size_t j ;
      T sum(0.0) ;
      //
      // Orthogonalize vector (Remove projections).
      //
      for (j = (size-1) ; j > size_t(i) ; j--) {
        sum = 0.0 ;
        for (size_t k = 0 ; k < size ; k++) {
          sum += (*this)[k][j] * (*this)[k][i] ;
        }

        for (size_t k = 0 ; k < size ; k++) {
          (*this)[k][i] -= sum * (*this)[k][j] ;
        }
      }

      //
      // Normalize vector.
      //
      sum = 0.0 ;
      size_t l ;
      for (l = 0 ; l < size ; l++) {
        sum += (*this)[l][i] * (*this)[l][i] ;
      }
      sum = sqrt(sum) ;
      for (l = 0 ; l < size ; l++) {
        (*this)[l][i] /= sum ;
      }

    }

  }

  //
  // Transpose matrix the current matrix.
  //
  template<class T>
  void TMatrix<T>::Transpose()  {

    size_t size = this->size() ;

    for (size_t i(0); i < size ; ++i) {

      for (size_t j(i+1); j < size ; ++j) {

        T tmp = (*this)[i][j] ;
        (*this)[i][j] = (*this)[j][i] ;
        (*this)[j][i] = tmp ;

      }

    }

  }

  template<class T>
  void TMatrix<T>::WriteToXML(PersistPtr pp)
  {
    //Write a CML element <matrix> as child of pp
    stringstream ss;
    for(size_t i(0) ; i < this->size() ; i++) {
      for(size_t j(0) ; j <= i ; j++)
        ss << to_double((*this)[i][j]) << ' ';
      ss << '\n'; 
    }
    PersistPtr ppmatrix = pp->XmlWriteValueElement("matrix", ss.str(),true);
    // The "true" parameter puts the matrix values in a CDATA wrapper so that
    // the new lines are preserved. If the parameter is omitted the data
    // is all space separated. Both form are read identically.
    ppmatrix->XmlWriteAttribute("matrixType", "squareSymmetricLT");
    ppmatrix->XmlWriteAttribute("rows", toString(this->size()));
  }

  template<class T>
  void TMatrix<T>::showFinalBits(const size_t n, bool isTabbed){

    size_t size = this->size() ;

    // if n == 0, print the whole matrix
    size_t fb(n);
    if (n == 0) fb = size;

    //
    // Show the final n x n square of the current matrix
    //
    ctest << "{\n";
    if (!isTabbed){
      for (size_t i = size - fb ; i < size ; ++i ) {
        for (size_t j = size - fb ; j < size ; ++j ) {
          formatFloat(ctest, (*this)[i][j], 5,  13) ;
        }
        ctest << endl;
      }
    } else {
      for (size_t i = size - fb ; i < size ; ++i ) {
        for (size_t j = size - fb ; j < size ; ++j ) {
          ctest << (*this)[i][j] << "\t";
        }
        ctest << endl;
      }
    }
    ctest << "}\n";
  }

  // Read a square CML <matrix>:
  // The data type of the values is T and they are separated by whitespace.
  // ppmatrix points to the <matrix> element, which must have an attribute
  // rows="n" where n is a non-zero integer.
  // If it has an attribute matrixType="squareSymmetricUT" the values are the upper triangle;
  // if it is "squareSymmetricLT" they are the lower triangle; and 
  // if it is omitted or is anything else all elements are assumed present.
  // The returned matrix is fully populated.
  template<class T>
  TMatrix<T>* ReadMatrix(PersistPtr ppmatrix)  {
    size_t nrows(0);
    if(ppmatrix && (nrows = ppmatrix->XmlReadInteger("rows",false))!=0)
    {
      bool upper=false, lower=false;
      const char* ptype = ppmatrix->XmlReadValue("matrixType",false);
      if(strcmp(ptype,"squareSymmetricLT")==0) 
        lower=true;
      else if(strcmp(ptype,"squareSymmetricUT")==0)
        upper=true;

      TMatrix<T>* m = new TMatrix<T>(nrows);
      std::stringstream ss;
      ss.str(ppmatrix->XmlRead());
      for(size_t nr=0; nr<nrows; ++nr)
      {
        size_t a = upper ? nr : 0;
        size_t b = lower ? nr : nrows-1;
        for(size_t nc=a; nc<=b; ++nc)
        {
          double val;
          ss >> val;
          (*m)[nr][nc] = (*m)[nc][nr] = val;
        }
      }
      return m;
    }
    return NULL; //empty
  }

  // Read a symmetrical square matrix from a CML <property> element, given
  // a pointer to the molecule or propertyList. See ReadMatrix() for details.
  template<class T>
  TMatrix<T>* ReadPropertyMatrix(const string& name, PersistPtr ppparent) {
    PersistPtr ppmatrix = ppparent->XmlMoveToProperty(name); //to <matrix>
    return ReadMatrix<T>(ppmatrix);
  }


}//namespace mesmer


#endif // GUARD_TMatrix_h
