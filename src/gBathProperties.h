#ifndef GUARD_gBathProperties_h
#define GUARD_gBathProperties_h


#include "MolecularComponents.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  class gBathProperties:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Bath gas related properties
    //-------------------------------------------------------------------------------------------------

  private:
    double         m_Sigma ;            // Lennard-Jones sigma.
    double         m_Epsilon ;          // Lennard-Jones epsilon.

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_Sigma_chk;
    int m_Epsilon_chk;
    //================================================

  public:

    //
    // Constructor, destructor and initialization
    //
    gBathProperties(Molecule* pMol);
    virtual ~gBathProperties();

    double getSigma() ;
    double getEpsilon() ;
    void   setSigma(double value);
    void   setEpsilon(double value);
  };
}
#endif
