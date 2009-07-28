#ifndef GUARD_System_h
#define GUARD_System_h

//-------------------------------------------------------------------------------------------
//
// System.h
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This header file contains the declaration of the System class. This class is the route
// of all molecular and reaction data will contain information about any number of System.
// Reaction System inforamtion will be sorted in a vector of reaction maps. Molecular data
// is stored in the molecule manager.
//
//-------------------------------------------------------------------------------------------

#include "ReactionManager.h"
#include "calcmethod.h"

namespace mesmer
{
  class CalcMethod;
  class System 
  {
  public:

    System(const std::string& libraryfilename) ;
    ~System() ;

    // Read and parse a data input file.
    bool parse(PersistPtr ppIOPtr) ;

    // Begin calculation.
    void executeCalculation() ;

    // Begin fitting.
    void fitting(void);

    // Begin single calculation.
    bool calculate(double& chiSquare) ;
  
    // Print system configuration
    void configuration(void);

    // Deduce the me:type of each molecule and add it to the XML file 
    bool assignMolTypes(PersistPtr ppIOPtr) ;
    
    void WriteMetadata();

    // Mesmer control flags.
    MesmerFlags m_Flags;

  private:


    void readPTs(PersistPtr);

    bool ReadRange(const std::string&    name,
      std::vector<double>&  vals,
      PersistPtr            ppbase,
      bool                  MustBeThere=true);

    // Location of the molecule manager.
    MoleculeManager *m_pMoleculeManager;

    // Location of the reaction mananger.
    ReactionManager *m_pReactionManager ;

    // Paired concentration and pressure points.
    std::vector<CandTpair> PandTs;

    // Physical variables. Reference to this are passed to all Molecule and Reaction constructors.
    MesmerEnv m_Env;
    
    // level in XML file under <mesemer>
    PersistPtr m_ppIOPtr;

    CalcMethod* m_CalcMethod;

	// Linesearch. Performs an optimization along a specified direction.
	// This is only a temporary location for this method, it should be
	// located in a separate class that deals with optimization.

	void LineSearch(const int varID);

    const char* m_pTitle;
  } ;
}//namespace
#endif // GUARD_System_h
