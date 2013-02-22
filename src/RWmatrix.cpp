#include "RWMatrix.h"
#include "TimeCounter.h"

namespace mesmer
{

  
void RWMatrix::WriteHessianToXML(dMatrix& hessian, PersistPtr pp)
{
  //Write a CML element <matrix> as child of pp
  stringstream ss;
  for(size_t i(0) ; i < hessian.size() ; i++) {
    for(size_t j(0) ; j <= i ; j++)
      ss << hessian[i][j] << ' ';
    ss << '\n'; 
  }
  PersistPtr ppmatrix = pp->XmlWriteValueElement("matrix", ss.str(),true);
  // The true parameter puts the matrix values in a CDATA wrapper so that
  // the new lines are preserved. If the parameter is omitted the data
  // ia all space separated. Both form are read identically.
  ppmatrix->XmlWriteAttribute("matrixType", "squareSymmetricLT");
  ppmatrix->XmlWriteAttribute("rows", toString(hessian.size()));
}

// Reads a matrix in CML property, given a pointer to the molecule
// or propertyList. See ReadMatrix() for details.
dMatrix RWMatrix::ReadPropertyMatrix(const string& name, PersistPtr ppparent)
{
  PersistPtr ppmatrix = ppparent->XmlMoveToProperty(name); //to <matrix>
  return ReadMatrix(name, ppmatrix);
}

// With ppmatrix at the <matrix> element, reads a square, symmetric matrix
// with double values from a list of its whitespace separated lower triangle
// values. The <matrix> element must have attributes
//  matrixType="squareSymmetricLT" rows="n" where n is a non-zero integer.
dMatrix RWMatrix::ReadMatrix(const string& name, PersistPtr ppmatrix)
{
  unsigned nrows=0;
  if(ppmatrix
    && strcmp(ppmatrix->XmlReadValue("matrixType",false),"squareSymmetricLT")==0
    && (nrows = ppmatrix->XmlReadInteger("rows",false))!=0)
  {
    dMatrix m(nrows);
    std::stringstream ss;
    ss.str(ppmatrix->XmlRead());
    for(unsigned nr=0; nr<nrows; ++nr)
    {
      for(unsigned nc=0; nc<=nr; ++nc)
      {
        double val;
        ss >> val;
        m[nr][nc] = m[nc][nr] = val;
      }
    }
    return m;
  }
  cwarn << "Could not read the matrix " << name << endl;
  return dMatrix(0); //empty
}
}