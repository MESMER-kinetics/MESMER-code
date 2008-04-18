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
  // defining the unit concentration units conversion rule
  typedef std::map<std::string, int> ConcentrationUnitConversionMapType;


  // to make sure if there is a concentration or pressure definition, there is a temperature definition.
  struct CandTpair{
    CandTpair(double _cp, double _t): concentration(_cp), temperature(_t){}
    double concentration; // particles per cubic centimeter
    double temperature; // Kelvin
  };

  class System
  {
  public:

    System() ;
    ~System() ;

    // Read and parse a data input file.
    bool parse(PersistPtr ppIOPtr) ;

    // Begin calculation.
    void calculate() ;


    std::vector<CandTpair> CPandTs;
    //    std::vector<double> Temperatures;   //Explicit emumeration of the temperatures requested in the input
    //    std::vector<double> Concentrations; //Bath gas concentrations: values for calculation. Use if Pressures is empty.
    //    std::vector<double> Pressures;      //Bath gas pressures: values for calculation

    //Stores environmental variables
    //Reference to this are passed to the constructors of all Molecules and Reactions
    MesmerEnv m_Env;

  private:

    // Location of the molecule manager.
    MoleculeManager *m_pMoleculeManager ;

    // Location of the reaction maps.
    ReactionManager *m_pReactionManager ;

    double getConvertedCP(string unitInput, double concentrationInput, double temperatureInput);
    bool SetGrainParams();
    void readCPTs(PersistPtr);
    bool ReadRange(const std::string&    name,
      std::vector<double>&  vals,
      PersistPtr            ppbase,
      bool                  MustBeThere=true);

    void WriteMetadata();

    // level in XML file under <mesemer>
    PersistPtr m_ppIOPtr;

    // mapping the conversion of concentration, pressure
    ConcentrationUnitConversionMapType concentrationUnitConversionMap;
  } ;
}//namespace
#endif // GUARD_System_h
