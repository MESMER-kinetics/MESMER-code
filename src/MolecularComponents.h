//-------------------------------------------------------------------------------------------
// MolecularComponents.h
//
// Author: Chi-Hsiu Liang
//
// This file contains the base class defintion of property groups used in the Molecule class. 
//-------------------------------------------------------------------------------------------

#ifndef GUARD_MolecularComponents_h
#define GUARD_MolecularComponents_h

#include <memory>
#include "MesmerEnv.h"
#include "MesmerFlags.h"
#include "Rdouble.h"
#include "Persistence.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{

  enum RotationalTop {
    ATOMIC,
    LINEAR,
    NONLINEAR,
    SPHERICAL,
    OBLATE,
    PROLATE,
    ASYMMETRIC,
    UNDEFINED_TOP
  } ;

  // Forward class declarations.
  class Molecule;

  class MolecularComponent{
  public:
    Molecule* getHost() { return m_host; }
    const Molecule* getHost() const { return m_host; }
    static void setEnergyConvention(const std::string& convention){m_energyConvention=convention;}
    static string getEnergyConvention(){ return m_energyConvention; }

  protected:
    Molecule* m_host;
    static std::string m_energyConvention; // For all molecules.
  };

}//namespace

#endif // GUARD_MolecularComponents_h
