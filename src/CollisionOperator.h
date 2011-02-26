#ifndef GUARD_CollisionOperator_h
#define GUARD_CollisionOperator_h

//-------------------------------------------------------------------------------------------
//
// CollisionOperator.h
//
// Author: Struan Robertson
// Date:   26/Feb/2011
//
// This header file contains the declaration of the CollisionOperator class.
// This class will implement the master equation collision operator.
//
//-------------------------------------------------------------------------------------------

#include "dMatrix.h"

namespace mesmer
{

  class CollisionOperator
  {
  public:

    // Constructor
    CollisionOperator() :
        m_reactionOperator(0),
          m_eigenvectors(0),
          m_eigenvalues() {} ;

        // Destructor.
        virtual ~CollisionOperator() {} ;

  private:

    // The system transition matrix and associated eigenvalues and eigenvectors.
    qdMatrix               *m_reactionOperator ;
    qdMatrix               *m_eigenvectors;
    std::vector<qd_real>    m_eigenvalues;

  } ;

}
#endif // GUARD_CollisionOperator_h
