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
// of all molecular and reaction data and contains information about any number of species.
// Reaction System information will be stored in a vector of reaction maps. Molecular data
// is stored in the molecule manager.
//
//-------------------------------------------------------------------------------------------
#include <set>
#include "plugin.h"
#include "ReactionManager.h"
#include "calcmethod.h"
#include "CollisionOperator.h"
#include "ParallelManager.h"

#define MESMER_VERSION "6.0"

namespace mesmer
{
  ///Search the library of molecules. If PPMolList is valid copy to the main XML file,
  // replacing an existing molecule of same name. Return a pointer to the copy (or the librymol).
  extern PersistPtr GetFromLibrary(const std::string molName, PersistPtr ppMolList);
  //global variable
  extern std::string libfile;

	class ConditionsManager;

  class System 
  {
  public:

    System(const std::string& libraryfilename, ParallelManager *pParallelManager) ;
    ~System() ;

    // Initialize the System object.
    bool initialize(void) ;

    // Read and parse a data input file.
    bool parse(PersistPtr ppIOPtr) ;

    // Begin calculation.
    void executeCalculation() ;

    // Begin single calculation.
    bool calculate(double& chiSquare, vector<double> &residuals, bool writeReport = false) ;

    // Begin single calculation - wrapper function.
    bool calculate(double& chiSquare, bool writeReport = false) {
      vector<double> residuals ;
      return calculate(chiSquare, residuals, writeReport) ;
    } ;

    // Calculate rate coefficients etc. for a specific condition.
    bool calculate(size_t nCond, vector<double> &Quantities, bool writeReport = false) ;

    // Begin single calculation.
    bool calculate(double &Temperature, double &Concentration, Precision precision,
      map<string, double> &phenRates, double &MaxT, const std::string& bathGas);

    // Wrapper for single calculation.
    bool calculate(size_t nCond, map<string, double> &phenRates)  ;

    // Print system configuration
    void configuration(void);

    // Deduce the role of each molecule and add it to the XML file 
    bool assignMolTypes(PersistPtr ppIOPtr) ;

    void WriteMetadata(const std::string& infilename);

    PersistPtr getPersistPtr() { return m_ppIOPtr; }

    MoleculeManager*   getMoleculeManager()  { return m_pMoleculeManager; } ;
    ReactionManager*   getReactionManager()  { return m_pReactionManager; } ;
    ConditionsManager* getConditionsManager(){ return m_pConditionsManager; } ;
		ParallelManager*   getParallelManager()  { return m_pParallelManager; } ;

    // Mesmer control flags.
    MesmerFlags m_Flags;

    // for each <control> block
    std::vector<MesmerFlags> m_FlagsForEachControl;
    std::vector<CalcMethod*> m_CalcMethodsForEachControl;
    std::vector<ConditionsManager*> m_ConditionsForEachControl;

    MesmerEnv& getEnv() { return m_Env; } ;

    //Visit each set of Conditions to collect bath gas names
    void getAllBathGases(std::set<std::string>& bathGases);

  private:

    double calcChiSqRateCoefficients(const qdMatrix& mesmerRates, const unsigned calPoint, vector<double> &residuals) ;
    double calcChiSqYields(const unsigned calPoint, vector<double> &residuals);
    double calcChiSqEigenvalues(const unsigned calPoint, vector<double> &residuals);
    double calcChiSqRawData(const unsigned calPoint, vector<double> &residuals);

    // Location of the molecule manager.
    MoleculeManager *m_pMoleculeManager;

    // Location of the reaction mananger.
    ReactionManager *m_pReactionManager ;

		// Location of the conditions manager.
    ConditionsManager *m_pConditionsManager;

		// Location of the parallel mananger.
		ParallelManager *m_pParallelManager;

    // Physical variables. Reference to this are passed to all Molecule and Reaction constructors.
    MesmerEnv m_Env;

    // level in XML file under <mesemer>
    PersistPtr m_ppIOPtr;

    // The method used for the main calculation
    CalcMethod* m_CalcMethod;

    const char* m_pTitle;
    const char* m_pDescription;

		// Collision operator instance.

    CollisionOperator m_collisionOperator ;

  } ;

}//namespace
#endif // GUARD_System_h
