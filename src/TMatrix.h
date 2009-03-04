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
    TMatrix( size_t n, const T& init = T()) : Matrix<T>(n, init) { } ;

    //
    // Wrapped call to EISPACK routine to diagonalise matrix.
    //
    void diagonalize(T *rr) {

      size_t size = this->size() ;

      //  Allocate memory for work array
      T *work = new T[size] ;
      T *rrProxy = new T[size] ;

      tred2(this->m_matrix, size, rrProxy, work) ;
      tqli(rrProxy, work, size, this->m_matrix) ;

      for (size_t i = 0; i < size; ++i){
        rr[i] = rrProxy[i];
      }

      delete [] work ;
      delete [] rrProxy ;

    }

    //
    // Solve a set of linear equations with a single right hand side.
    //
    void solveLinearEquationSet(T *rr) {

      size_t size = this->size() ;

      //  Allocate memory for work array
      int *indx = new int[size] ;

      if (ludcmp(this->m_matrix, size, indx)){
        exit(1);
      }

      lubksb(this->m_matrix, size, indx, rr) ;

      delete [] indx ;

    };

    // Matrix inversion method by Gaussian elimination
    int invertGaussianJordan();

    // Matrix inversion method by LU decomposition
    int invertLUdecomposition();

    // Matrix inversion method by adjoint cofactors
    int invertAdjointCofactors();

    void normalizeProbabilityMatrix();

  private:

    //
    // EISPACK methods for diagonalizing matrix.
    //
    void    tred2   (T **a, size_t n, T *d, T *e) ;
    void    tqli    (T *d, T *e, size_t n, T **z) ;
    T  pythag  (T a, T b) ;

    //
    // NR LU methods for linear equation solving.
    //
    int ludcmp(T **a,  size_t n, int *indx) ;
    void lubksb(T **a,  size_t n, int *indx, T* b) ;

    //
    // Calculate the inverse of the matrix by finding the adjoint of the cofactors matrix
    int GetMinor(T **src, T **dest, int row, int col, int order);
    T CalcDeterminant( T **mat, int order);

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
  void TMatrix<T>::tred2(T **a, size_t n, T *d, T *e)
  {
    size_t l, k, j, i;
    T scale, hh, h, g, f;

    for (i=n;i>=2;--i) {
      l=i-1;
      h=scale=0.0;
      if (l > 1) {

        for (k=1;k<=l;++k)
          scale += abs(a[i-1][k-1]);
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
          dd=abs(d[m-1])+abs(d[m]);
          if (abs(e[m-1])+dd == dd) break;
        }
        if (m != l) {
          // if (iter++ == 30) fprintf(stderr, "Too many iterations in TQLI");
          if (iter++ == 60) fprintf(stderr, "Too many iterations in TQLI");
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
          g=d[m-1]-d[l-1]+e[l-1]/(g + (g < 0.0 ? -abs(r) : abs(r)));
          s=c=1.0;
          p=0.0;
          for (i=m-1;i>=l;--i) {
            f=s*e[i-1];
            b=c*e[i-1];
            if (abs(f) >= abs(g)) {
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
  // EISPACK pythag function.
  //
  // Finds sqrtl(a**2+b**2) without overflow or destructive underflow.
  //
  //-------------------------------------------------------------------------------------------
  template<class T>
  T TMatrix<T>::pythag(T a, T b)
  {
    T p,r,s,t,u;
    if (abs(a) > abs(b)) p = abs(a);
    else p = abs(b);
    if (p == 0.) goto label_1;
    if (abs(a) > abs(b)) r = abs(b);
    else r = abs(a);
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
  int TMatrix<T>::ludcmp(T **a,  size_t n, int *indx) {

    int imax;

    T big, dum, sum, temp ;
    T tiny = numeric_limits<T>::epsilon();

    T *work = new T[n] ;

    for (size_t i(0); i < n ; ++i) {
      big = 0.0 ;
      for (size_t j(0); j < n ; ++j) {
        if ((temp = abs(a[i][j])) > big){
          big = temp ;
        }
      }
      if (big == 0.0) {
        cerr << "Singular Matrix in routine ludcmp";
        return 1;
      }
      work[i] = 1.0/big ;
    }

    for (size_t j(0); j < n ; ++j) {
      for (size_t i(0); i < j ; ++i) {
        sum = a[i][j] ;
        for (size_t k(0); k < i ; ++k){
          sum -= a[i][k]*a[k][j] ;
        }
        a[i][j] = sum ;
      }
      big = 0.0 ;
      for (size_t i(j); i < n; ++i) {
        sum = a[i][j] ;
        for (size_t k(0); k < j ; ++k)
          sum -= a[i][k]*a[k][j] ;

        a[i][j] = sum ;

        if ( (dum = work[i]*abs(sum)) >= big) {
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
      if (abs(a[j][j]) < tiny){
        a[j][j] = tiny;
      }

      if (j != n-1) {
        dum = 1.0/(a[j][j]) ;
        for (int i(j+1); i < n; ++i)
          a[i][j] *= dum ;
      }

    }

    delete [] work ;
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
  void TMatrix<T>::lubksb(T **a,  size_t n, int *indx, T* b) {

    int ii = 0, ip;
    T sum ;

    for (size_t i(0); i < n; ++i) {
      ip = indx[i] ;
      sum = b[ip] ;
      b[ip] = b[i] ;
      if (ii >= 0) {
        for (size_t j(ii); j < i; ++j)
          sum -= a[i][j]*b[j] ;
      }
      else if (sum != 0.0){
        ii = i ;
      }
      b[i] = sum ;
    }
    for (size_t i(n-1); i >= 0; --i) {

      sum = b[i] ;
      if (i < n-1){
        for (size_t j(i+1); j < n; ++j)
          sum -= a[i][j]*b[j] ;

      }
      b[i] = sum/a[i][i] ;
    }
  }

  template<class T>
  int TMatrix<T>::invertLUdecomposition(){
    int size = static_cast<int>(this->size()) ;

    //  Allocate memory for work array
    int *indx = new int[size] ;

    Matrix<T> invM(size); // an identity matrix as a primer for the inverse
    for (int i(0); i < size; ++i){
      for (int j(0); j < size; ++j){
        invM[i][j] = 0.0;
      }
      invM[i][i] = 1.0;
    }

    int rc = ludcmp(this->m_matrix, size, indx) ;

    ctest << "After ludcmp:";
    this->showFinalBits(size, true);
    //call solver if previous return code is ok
    //to obtain inverse of A one column at a time
    if (rc == 0) {
      T *temp = new T[size] ;
      for (int j(0); j < size; ++j) {
        for (int i(0); i < size; ++i) temp[i] = invM[i][j];
        lubksb(this->m_matrix, size, indx, temp);
        for (int i(0); i < size; ++i) invM[i][j] = temp[i];
      }
      for (int j(0); j < size; ++j) {
        for (int i(0); i < size; ++i){
          this->m_matrix[i][j] = invM[i][j];
        }
      }
      delete [] temp;
      delete [] indx ;
      return 0;
    }
    else{
      delete [] indx;
      return 1;
    }

  }

  template<class T>
  int TMatrix<T>::invertGaussianJordan(){
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
    ######################################################################*/
    int n = static_cast<int>(this->size()) ;
    if (n > INT_MAX) return 2; // The matrix size is too large to process.
    T divide, ratio;

    //-------------------------------
    //produce a unit vector of size n, and copy the incoming matrix
    Matrix<T> m1(n);
    Matrix<T> m2(n);
    for (int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        m1[i][j] = this->m_matrix[i][j];
        m2[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }
    //-------------------------------
    int* zeroCount = new int[n];

    for(int j = 0; j < n; ++j){
      for (int i = 0; i < n; ++i){
        if (m1[i][j] == 0.0){
          zeroCount[j]++; if (zeroCount[j] == n) return 1;
          /* If there is a zero column, the matrix has no inverse.*/
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
    delete [] zeroCount;

    for(int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
        this->m_matrix[i][j] = m2[i][j];
      }
    }

    return 0;
  }

  template<class T>
  int TMatrix<T>::invertAdjointCofactors(){
    // get the determinant of m_matrix
    int order = static_cast<int>(this->size()) ;
    T det = 1.0/CalcDeterminant(this->m_matrix,order);
    Matrix<T> Y(order);

    // memory allocation
    T *temp = new T[(order-1)*(order-1)];
    T **minor = new T*[order-1];
    for(int i=0;i<order-1;++i)
      minor[i] = temp+(i*(order-1));

    for(int j=0;j<order;++j){
      for(int i=0;i<order;++i){
        // get the co-factor (matrix) of m_matrix(j,i)
        GetMinor(this->m_matrix,minor,j,i,order);
        Y[i][j] = det*CalcDeterminant(minor,order-1);
        if( (i+j)%2 == 1)
          Y[i][j] = -Y[i][j];
      }
    }

    for(int j=0;j<order;++j){
      for(int i=0;i<order;++i){
        this->m_matrix[i][j] = Y[i][j];
      }
    }

    // release memory
    delete [] minor[0];
    delete [] minor;
    return 0;
  }

  // calculate the cofactor of element (row,col)
  template<class T>
  int TMatrix<T>::GetMinor(T **src, T **dest, int row, int col, int order)
  {
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; ++i )
    {
      if( i != row )
      {
        colCount = 0;
        for(int j = 0; j < order; ++j )
        {
          // when j is not the element
          if( j != col )
          {
            dest[rowCount][colCount] = src[i][j];
            colCount++;
          }
        }
        rowCount++;
      }
    }

    return 1;
  }

  // Calculate the determinant recursively.
  template<class T>
  T TMatrix<T>::CalcDeterminant( T **mat, int order)
  {
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
      return mat[0][0];

    // the determinant value
    T det = 0;

    // allocate the cofactor matrix
    T **minor;
    minor = new T*[order-1];
    for(int i=0;i<order-1;++i)
      minor[i] = new T[order-1];

    for(int i = 0; i < order; ++i )
    {
      // get minor of element (0,i)
      GetMinor( mat, minor, 0, i , order);
      // the recusion is here!
      det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
    }

    // release memory
    for(int i=0;i<order-1;++i)
      delete [] minor[i];
    delete [] minor;

    return det;
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
    vector<double> work(optrsize) ;// Work space.

    double scaledRemain(0.0) ;
    for ( i = optrsize - 1 ; i >= 0 ; --i ) {

      double upperSum(0.0) ;
      for ( j = 0 ; j <= i ; ++j )
        upperSum += (*this)[j][i] ;

      if (upperSum > 0.0){
        if (i < optrsize - 1){
          scaledRemain = 0.0;
          for ( j = i + 1 ; j < optrsize ; ++j ){
            double scale = work[j];
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
      double value = (*this)[i][i];
      for ( j = i + 1 ; j < optrsize ; ++j ) {
        (*this)[j][i] *= work[j] ;
        (*this)[i][j] *= work[j] ;
      }
    }

  }

}//namespacer mesmer


#endif // GUARD_TMatrix_h
