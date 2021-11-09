//-------------------------------------------------------------------------------------------
// MolecularComponents.h
//
// Author: Chi-Hsiu Liang
//
// This file contains property groups of class Molecule. These groups give molecules variables
// and functions to perform tasks; from the definitions of these groups a molecule can play
// roles when it is required to do so. Classes in this file do not depend on each other and
// thus they can be seperated. Any of them can be added into a molecule (with a new() to construct
// an object and then pass the pointer to the molecule) when the role of the molecule requires
// the information in that group.
//-------------------------------------------------------------------------------------------

#ifndef GUARD_MolecularComponents_h
#define GUARD_MolecularComponents_h

#include <memory>
#include "DensityOfStates.h"
#include "Distribution.h"
#include "MesmerEnv.h"
#include "MesmerFlags.h"
#include "Rdouble.h"
#include "EnergyTransferModel.h"
#include "vector3.h"
#include "dMatrix.h"
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
    //MolecularComponent():m_host(NULL){}
    static std::string m_energyConvention; //for all molecules
  };

  //-------------------------------------------------------------------------------------------------
  // Other related functions
  //-------------------------------------------------------------------------------------------------

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS(gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, std::vector<double>& rctsCellDOS);

}//namespace

#endif // GUARD_MolecularComponents_h
