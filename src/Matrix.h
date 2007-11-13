#ifndef GUARD_Matrix_h
#define GUARD_Matrix_h

//-------------------------------------------------------------------------------------------
//
// Matrix.h
//
// Author: Struan Robertson
// Date:   23/Jan/2003
//
// This header file contains the declaration of the Matrix class.  Taken from a combination
// of sources including Numerical Recipies, TNT and Accelerated C++, BUT restricted to square
// matrices only.
//
//-------------------------------------------------------------------------------------------

//#include <stdexcept>
#include <sstream>
#include "oberror.h"

using namespace std;

namespace mesmer
{
template<class T>
class Matrix {

public:

    // Type defs

    typedef         size_t  size_type ;
    //typedef             T*  iterator ;
    //typedef const       T*  const_iterator ;
    //typedef             T   value_type ;

    // Constructors

    explicit Matrix(size_type n, const T& a = T() ) ;
    Matrix(const Matrix&) ; //Copy constructor

    // Destructor

    virtual ~Matrix() { destroy() ; }

    // Operators

    Matrix& operator=(const Matrix& rhs) ;
    T* operator[](const size_type i) { return m_matrix[i] ; }
    const T* operator[](const size_type i) const { return m_matrix[i] ; }

    // Accessors

    size_type size() const { return m_msize ; }

    // Modifiers
    void resize(const size_type n) ;
    void normalize();

protected:

    // Size of matrix.

    size_type m_msize ;

    // Pointer to location of matrix.

    T **m_matrix ;

    // Internal function that creates matrix.

    void create(size_type n) {

        m_msize = n ;

        m_matrix = allocatematrix(m_msize) ;

    }

  // Internal function that allocates space for matrix.

    T** allocatematrix(size_type n) {

        T ** matrix = new T*[n] ;

        matrix[0] = new T[n*n] ;

        for ( size_type i = 1 ; i < n ; i++ )
            matrix[i] = matrix[i - 1] + n ;

    return matrix ;

    }

    // Internal function that destroys matrix

    void destroy() {

        if (m_matrix != NULL) {
            delete[] (m_matrix[0]) ;
            delete[] (m_matrix) ;
        }
    }

private:

    // Hide default constructor - force the size of the matrix to be passed.

    Matrix() : m_msize(0), m_matrix(0) { } ;

} ;

// Construct a matrix of size n and initialized to a.

template<class T>
Matrix<T>::Matrix(size_type n, const T& a) : m_msize(0), m_matrix(0) {

    create(n) ;

    for ( size_type i = 0 ; i < m_msize ; i++ )
        for ( size_type j = 0 ; j < m_msize ; j++ )
            m_matrix[i][j] = a ;
}

// Copy constructor.

template<class T>
Matrix<T>::Matrix(const Matrix& rhs) : m_msize(0), m_matrix(0) {

    create(rhs.m_msize) ;

    for ( size_type i = 0 ; i < m_msize ; i++ )
        for ( size_type j = 0 ; j < m_msize ; j++ )
            m_matrix[i][j] = rhs[i][j] ;
}

// Operators

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &rhs) {

    if (this != rhs) {

        // If necessary resize underlying array.

        if (m_msize != rhs.m_msize) {

            destroy() ;

            create(rhs.m_msize) ;

        }

        for ( size_type i = 0 ; i < m_msize ; i++ )
            for ( size_type j = 0 ; j < m_msize ; j++ )
                m_matrix[i][j] = rhs[i][j] ;
    }
    return *this ;
}

// Modifiers.

template<class T>
void Matrix<T>::resize(const size_type n){

    if (n < 1){
      stringstream errorMsg;
      errorMsg << "Matrix must of size one or greater";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    T **matrix = allocatematrix(n)  ;

    size_type msize = std::min(n, m_msize) ;
    for ( size_type i = 0 ; i < msize ; i++ )
        for ( size_type j = 0 ; j < msize ; j++ )
            matrix[i][j] = m_matrix[i][j] ;

    destroy() ;

    m_msize = n ;

    m_matrix = matrix;

}

template<class T>
void Matrix<T>::normalize(){
    
    T* work = new T[m_msize] ;// Work space.
    //
    // Normalization of Probability matrix.
    // Normalising coefficients are found by using the fact that column sums
    // are unity. The procedure leads to a matrix that is of upper triangular
    // form and the normalisation constants are found by back substitution.
    //
    
    int i, j; //int makes sure the comparison to negative numbers meaningful (i >=0)

    for ( i = (int)m_msize - 1 ; i >= 0 ; --i ) {

      T upperSum(0.0) ;
      for ( j = 0 ; j <= i ; ++j ) upperSum += m_matrix[j][i] ;

      T scaledRemain(0.0) ;
      for ( j = i + 1 ; j < (int)m_msize ; ++j ) scaledRemain += m_matrix[j][i] * work[j] ;

      if (upperSum <= 0.0) {
        stringstream errorMsg;
        errorMsg << "Normalization coefficients in this matrix is smaller than or equal to zero";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        exit(1) ;
      }
      work[i] = (1.0 - scaledRemain) / upperSum ;
    }

    //
    // Apply normalization coefficients
    //
    for ( i = 0 ; i < (int)m_msize ; ++i ) {
      m_matrix[i][i] *= work[i] ;
      for ( j = i + 1 ; j < (int)m_msize ; ++j ) {
        m_matrix[j][i] *= work[j] ;
        m_matrix[i][j] *= work[j] ;
      }
    }
    delete [] work;

}
}//namespacer mesmer

#endif // GUARD_Matrix_h
