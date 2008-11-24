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
#include "formatfloat.h"

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
    void showFinalBits(const size_type n, bool isTabbed = false);
    void reset(const size_type n);

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

    if (this != &rhs) {

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
      cerr << "Matrix must be of size one or greater";
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
void Matrix<T>::reset(const size_type n){

    if (n < 1){
      stringstream errorMsg;
      cerr << "Matrix must be of size one or greater";
    }

    T **matrix = allocatematrix(n)  ;

    for ( size_type i = 0 ; i < n ; i++ )
        for ( size_type j = 0 ; j < n ; j++ )
            matrix[i][j] = 0.0;

    destroy() ;

    m_msize = n ;

    m_matrix = matrix;

}

template<class T>
void Matrix<T>::showFinalBits(const size_type n, bool isTabbed){

    // if n == 0, print the whole matrix
    size_type fb(n);
    if (n == 0) fb = m_msize;

    //
    // Show the final n x n square of the current matrix
    //
    ctest << "{\n";
    if (!isTabbed){
      for (size_type i = m_msize - fb ; i < m_msize ; ++i ) {
        for (size_type j = m_msize - fb ; j < m_msize ; ++j ) {
          formatFloat(ctest, m_matrix[i][j], 5,  13) ;
        }
        ctest << endl;
      }
    }
    else{
      for (size_type i = m_msize - fb ; i < m_msize ; ++i ) {
        for (size_type j = m_msize - fb ; j < m_msize ; ++j ) {
          ctest << m_matrix[i][j] << "\t";
        }
        ctest << endl;
      }
    }
    ctest << "}\n";


}
}//namespacer mesmer

#endif // GUARD_Matrix_h
