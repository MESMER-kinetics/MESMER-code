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

template<class T>
class Matrix {

public:

   // Type defs

   typedef  size_t  size_type ;
   typedef       T* iterator ;
   typedef const T* const_iterator ;
   typedef       T  value_type ;

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

protected:

   // Size of matrix. 

   size_type m_msize ; 

   // Pointer to location of matrix. 

   T **m_matrix ;

   // Internal function that creates matrix.
   
   void create(size_type n) {

      m_msize = n ;

      m_matrix = new T*[m_msize] ;

      m_matrix[0] = new T[m_msize*m_msize] ;

      for ( size_type i = 1 ; i < m_msize ; i++ )
          m_matrix[i] = m_matrix[i - 1] + m_msize ;

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

#endif // GUARD_Matrix_h
