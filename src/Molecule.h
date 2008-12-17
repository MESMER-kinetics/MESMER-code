#ifndef GUARD_Molecule_h
#define GUARD_Molecule_h

//-------------------------------------------------------------------------------------------
//
// Molecule.h
//
// Author: Struan Robertson
// Date:   5/Jan/2003
//
// This header file contains the declaration of the Molecule class.
//
//-------------------------------------------------------------------------------------------

#include "MolecularComponents.h"

namespace mesmer
{

  //**************************************************
  /// Basic molecule: has name and some collision parameters.
  /// Used for bath gases and unmodelled product molecules.
  class Molecule
  {

    //-------------------------------------------------------------------------------------------------
    // Basic information of a molecule
    //-------------------------------------------------------------------------------------------------

  private:

    const MesmerEnv&     m_Env;
    MesmerFlags&         m_Flags;

    PersistPtr     m_ppPersist;         // Conduit for I/O
    std::string    m_Name ;             // Molecule name.
    std::string    m_Description;       // Longer description for the structure

    double         m_Mass ;             // Mass.

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_Mass_chk;
    //================================================


  public:

    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Component pointers
    friend class MolecularComponent; // Provide access privalage of MolecularComponent to private members of Molecule.
    gBathProperties*        g_bath;  // Component pointer for bath gas properties
    gDensityOfStates*       g_dos;   // Component pointer for cell density of state properties
    gTransitionState*       g_ts;    // Component pointer for transition state properties
    gPopulation*            g_pop;   // Component pointer for population and equilibrium fraction
    gCollisionProperties*   g_coll;  // Component pointer for collision down model properties
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    //
    // Constructor
    //
    Molecule(const MesmerEnv& Env, MesmerFlags& m_Flags) ;
    virtual ~Molecule();

    // Initialize Molecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    PersistPtr  get_PersistentPointer();

    std::string getDescription() const;
    const MesmerEnv& getEnv() const;

    std::string getName() const;
    void   setName(string name) ;

    MesmerFlags& getFlags();

    double getMass() ;
    void   setMass(double value);

    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Manipulators of component pointers
    //gBathProperties        getBathGasProperties()   { return *g_bath; }
    //gDensityOfStates       getCDOSProperties()      { return *g_dos;  }
    //gTransitionState       getTSProperties()        { return *g_ts  ; }
    //gPopulation            getPopulationProperties(){ return *g_pop ; }
    //gCollisionProperties   getCollisionProperties() { return *g_coll; }

    bool activateRole(string molType);
    bool roleIsActivated(string molType);
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  };

}//namespace

// CHECK FOR INPUTFILE PARAMETERS
// for these check values, initially the values are given -1. If not provided by user the value will remain as -1.
// If provided by the user the value will be increased to 0.
// During calcualtion, every inputfile parameter if called by any function, the class will check if the parameter
// is provided. If it is provided, the chk value will increase by 1. So if a chk value is 9, it is asked for 9 times.
// However, if the user did not provide the value and the values is asked. The program will stop and report error.
// Values with obvious default values are also accounted in this check but the program will not exit; only
// warning message will be given.

#endif // GUARD_Molecule_h
