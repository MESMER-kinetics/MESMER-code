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
#include "MesmerConfig.h"
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
    void calculate() ;


    std::vector<double> Temperatures;   //Explicit emumeration of the temperatures requested in the input
    std::vector<double> Concentrations; //Bath gas concentrations: values for calculation. Use if Pressures is empty.
    std::vector<double> Pressures;      //Bath gas pressures: values for calculation

    //Stores environmental variables
    //Reference to this are passed to the constructors of all Molecules and Reactions
    MesmerEnv m_Env;

  private:

    // Location of the molecule manager.
    MoleculeManager *m_pMoleculeManager ;

    // Location of the reaction maps.
    ReactionManager *m_pReactionManager ;

    bool SetGrainParams();
    bool ReadRange(const std::string&    name,
                   std::vector<double>&  vals,
                   PersistPtr            ppbase,
                   bool                  MustBeThere=true);

    void WriteMetadata();

    // level in XML file under <mesemer>
    PersistPtr m_ppIOPtr;
  } ;
}//namespace
#endif // GUARD_System_h
