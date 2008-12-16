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

namespace mesmer
{

  class System
  {
  public:

    System() ;
    ~System() ;

    // Read and parse a data input file.
    bool parse(PersistPtr ppIOPtr) ;

    // Begin calculation.
    void calculate(double& chiSquare) ;

    // Begin fitting.
    void fitting(void);

    // Begin grid search
    void gridSearch(void);
    
    // Print system configuration
    void configuration(void);
    
    /// Access a list of molecules which were not present in the data file
    /// but which were sucessfully recovered from the Library.
    vector<PersistPtr>& getLibraryMols(){return m_pMoleculeManager->getLibraryMols();};

  public:

    // Paired concentration and pressure points
    std::vector<CandTpair> PandTs;

    //Stores environmental variables
    //Reference to this are passed to the constructors of all Molecules and Reactions
    MesmerEnv m_Env;
    MesmerFlags m_Flags;

  private:

    // Location of the molecule manager.
    MoleculeManager *m_pMoleculeManager ;

    // Location of the reaction maps.
    ReactionManager *m_pReactionManager ;

    void readPTs(PersistPtr);

    bool ReadRange(const std::string&    name,
      std::vector<double>&  vals,
      PersistPtr            ppbase,
      bool                  MustBeThere=true);

    void WriteMetadata();


    // level in XML file under <mesemer>
    PersistPtr m_ppIOPtr;

    const char* m_pTitle;
  } ;
}//namespace
#endif // GUARD_System_h
