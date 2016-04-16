#ifndef GUARD_MeMatrix_h
#define GUARD_MeMatrix_h

//-------------------------------------------------------------------------------------------
//
// MeMatrix.h
//
// Author: Struan Robertson
// Date:   23/Jan/2003
//
// This header file contains the declaration of the MeMatrix class.  Taken from a combination
// of sources including Numerical Recipies, TNT and Accelerated C++, BUT restricted to square
// matrices only.
//
//-------------------------------------------------------------------------------------------

#include <sstream>
#include "error.h"
#include "formatfloat.h"

using namespace std;

namespace mesmer
{
  template<class T>
  class MeMatrix {

  public:

    // Type defs

    typedef         T*      iterator ;
    typedef const   T*      const_iterator ;
    typedef         T       value_type ;

    // Constructor

    explicit MeMatrix(size_t n, const T& a = T() ) ;

    //Copy constructor

    MeMatrix(const MeMatrix&) ; 

    // Destructor

    virtual ~MeMatrix() { destroy() ; }

    // Operators

    MeMatrix& operator=(const MeMatrix& rhs) ;
    inline T* operator[](const size_t i) { return m_matrix[i] ; }
    inline const T* operator[](const size_t i) const { return m_matrix[i] ; }

    // Accessors

    size_t size() const { return m_msize ; }

    // Modifiers
    void resize(const size_t n) ;
    void reset(const size_t n);

  protected:

    // Size of matrix.

    size_t m_msize ;

    // Pointer to location of matrix.

    T **m_matrix ;

    // Internal function that creates matrix.

    void create(size_t n) {

      m_msize = n ;

      m_matrix = allocatematrix(m_msize) ;

    }

    // Internal function that allocates space for matrix.

    T** allocatematrix(size_t n) {

      T ** matrix = new T*[n] ;

      matrix[0] = new T[n*n] ;

      for ( size_t i = 1 ; i < n ; i++ )
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

    MeMatrix() : m_msize(0), m_matrix(0) { } ;

  } ;

  // Construct a matrix of size n and initialized to a.

  template<class T>
  MeMatrix<T>::MeMatrix(size_t n, const T& a) : m_msize(0), m_matrix(0) {

    create(n) ;

    for ( size_t i = 0 ; i < m_msize ; i++ )
      for ( size_t j = 0 ; j < m_msize ; j++ )
        m_matrix[i][j] = a ;
  }

  // Copy constructor.

  template<class T>
  MeMatrix<T>::MeMatrix(const MeMatrix& rhs) : m_msize(0), m_matrix(0) {

    create(rhs.m_msize) ;

    for ( size_t i = 0 ; i < m_msize ; i++ )
      for ( size_t j = 0 ; j < m_msize ; j++ )
        m_matrix[i][j] = rhs[i][j] ;
  }

  // Operators

  template<class T>
  MeMatrix<T>& MeMatrix<T>::operator=(const MeMatrix<T> &rhs) {

    if (this != &rhs) {

      // If necessary resize underlying array.

      if (m_msize != rhs.m_msize) {

        destroy() ;

        create(rhs.m_msize) ;

      }

      for ( size_t i = 0 ; i < m_msize ; i++ )
        for ( size_t j = 0 ; j < m_msize ; j++ )
          m_matrix[i][j] = rhs[i][j] ;
    }
    return *this ;
  }

  // Modifiers.

  template<class T>
  void MeMatrix<T>::resize(const size_t n){

    if (n < 1){
      stringstream errorMsg;
      cerr << "MeMatrix must be of size one or greater";
    }

    T **matrix = allocatematrix(n)  ;

    size_t msize = std::min(n, m_msize) ;
    for ( size_t i = 0 ; i < msize ; i++ )
      for ( size_t j = 0 ; j < msize ; j++ )
        matrix[i][j] = m_matrix[i][j] ;

    destroy() ;

    m_msize = n ;

    m_matrix = matrix;

  }

  template<class T>
  void MeMatrix<T>::reset(const size_t n){

    if (n < 1){
      stringstream errorMsg;
      cerr << "MeMatrix must be of size one or greater";
    }

    T **matrix = allocatematrix(n)  ;

    for ( size_t i = 0 ; i < n ; i++ )
      for ( size_t j = 0 ; j < n ; j++ )
        matrix[i][j] = 0.0;

    destroy() ;

    m_msize = n ;

    m_matrix = matrix;

  }

}//namespacer mesmer

#endif // GUARD_MeMatrix_h
