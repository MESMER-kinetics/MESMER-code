#ifndef GUARD_RWMatrix_h
#define GUARD_RWMatrix_h

#include "Persistence.h"
#include "dMatrix.h"

namespace mesmer
{

class RWMatrix //temporary class; somewhere to put code under development
{
public:

  // Read a symmetrical sqauare matrix from a CML <property> element
  static dMatrix ReadPropertyMatrix(const string& name, PersistPtr ppparent);

  // Read a symmetrical sqauare matrix from a CML <matrix> element
  static dMatrix ReadMatrix(const string& name, PersistPtr ppmatrix);

  // Write a CML element <matrix> as child of pp
  static void WriteHessianToXML(dMatrix& hessian, PersistPtr pp);

};

}

#endif //GUARD_FittingUtils_h

