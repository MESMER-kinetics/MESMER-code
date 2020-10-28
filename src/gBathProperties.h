#ifndef GUARD_gBathProperties_h
#define GUARD_gBathProperties_h

#include "MolecularComponents.h"

namespace mesmer
{
  class gBathProperties:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Bath gas related properties
    //-------------------------------------------------------------------------------------------------

  private:

    double m_Sigma ;         // Lennard-Jones sigma.
    double m_Epsilon ;       // Lennard-Jones epsilon.
    bool   m_dafaultPrmtrs ; // Indicates if defaults paramters have been read. 

  public:

    //
    // Constructor, destructor and initialization
    //
    gBathProperties(Molecule* pMol);
    ~gBathProperties() {};

    double getSigma() { return m_Sigma; }
    double getEpsilon() { return m_Epsilon; }
    bool   dafaultLJParatmeters() { return m_dafaultPrmtrs; }
  };
}
#endif
