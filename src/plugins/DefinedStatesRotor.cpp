
//-------------------------------------------------------------------------------------------
//
// DefinedStatesRotor.cpp
//
// This file contains the implementation of the methods for importing and testing the 
// density of states calculated or determine by other means. An example where this class
// might be deployed are the rotary-electronic states of OH which are not well approximated
// by treating these modes as decoupled, the coupling being described by one of Hund's
// coupling cases.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../gDensityOfStates.h"
#include "../Molecule.h"
#include "../Constants.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class DefinedStatesRotor : public DensityOfStatesCalculator
  {
  public:
    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol,  const MesmerEnv& env);

    // Function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    DefinedStatesRotor(const char* id) : m_id(id), 
      m_energyLevels(), 
      m_degeneracies() 
    { Register(); }

    virtual ~DefinedStatesRotor() {}
    virtual const char* getID()  { return m_id; }
    virtual bool includesRotations(){return true;}
    virtual DefinedStatesRotor* Clone() { return new DefinedStatesRotor(*this); }

  private:
    const char* m_id;
    vector<double> m_energyLevels ;   // Energies of the states and their associated degeneracies.
    vector<int>    m_degeneracies ;

  } ;

  //-------------------------------------------------------------
  //Global instance, defining its id
  DefinedStatesRotor theDefinedStatesRotor("DefinedStatesRotors");
  //-------------------------------------------------------------

  //Read data from XML and store in this instance.
  bool DefinedStatesRotor::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
  {
    // Read in defined states data.

    PersistPtr pp = ppDOSC->XmlMoveTo("me:States") ;

    if (pp) {

      const char* p = pp->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";

      while(pp = pp->XmlMoveTo("me:State"))
      {
        double energy = pp->XmlReadDouble("energy", optional);
        energy = getConvertedEnergy(units, energy);
        m_energyLevels.push_back(energy) ;
        int degeneracy = pp->XmlReadInteger("degeneracy", optional);
        m_degeneracies.push_back(degeneracy) ;
      }

    } else {

      // Default : free rotor.

      cinfo << "No rotory-electronic states defined for " << gdos->getHost()->getName() << "." << endl;

      m_energyLevels.push_back(0.0) ;
      m_degeneracies.push_back(1) ;
    }

    return true;
  }

  //
  // Calculate quantum mechanical 1D rotor densities of states of a free 
  // rotor and convolve them with the main density of states.
  //
  bool DefinedStatesRotor::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
	const size_t MaximumCell = env.MaxCell ;

    vector<double> cellDOS(MaximumCell, 0.0) ;

    // Shift eigenvalues by the zero point energy and form initial denisity of states array.

    double zeroPointEnergy(m_energyLevels[0]) ;
    for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
      size_t nr = nint((m_energyLevels[i] - zeroPointEnergy)/env.CellSize) ;
      if (nr < MaximumCell) {
        cellDOS[nr] = double(m_degeneracies[i]) ;
      }
    }

    // Initialise density of states.   

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;

  }

  //
  // Provide a function to calculate contribution to canonical partition function.
  //
  void DefinedStatesRotor::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne)
  {
    double Qrot(0.0), Erot(0.0), vErot(0.0) ;

    double zeroPointEnergy(m_energyLevels[0]) ; 
    for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
      double ene = m_energyLevels[i] - zeroPointEnergy ;
      Qrot  += double(m_degeneracies[i])*exp(-beta*ene) ;
      Erot  += ene*double(m_degeneracies[i])*exp(-beta*ene) ;
      vErot += ene*ene*double(m_degeneracies[i])*exp(-beta*ene) ;
    }
    Erot /= Qrot ;
    vErot = vErot/Qrot - Erot*Erot ;

    PrtnFn   *= Qrot ;
    IntrlEne += Erot ;
    varEne   += vErot ;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int DefinedStatesRotor::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    unsigned int nDOF(0) ;
    switch(rotorType){
      case NONLINEAR:
        nDOF = 3 ;
        break;
      case LINEAR:
        nDOF = 2 ;
        break;
      default:
        // Assume atom.
        break; 
    }

    return nDOF ;
  }

}//namespace
