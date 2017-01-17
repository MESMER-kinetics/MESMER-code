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
#include <set>
#include "plugin.h"
#include "ReactionManager.h"
#include "calcmethod.h"
#include "CollisionOperator.h"
#include "ConditionsManager.h"
#include "ParallelManager.h"

#define MESMER_VERSION "5.0"

namespace mesmer
{
  ///Search the library of molecules. If PPMolList is valid copy to the main XML file,
  // replacing an existing molecule of same name. Return a pointer to the copy (or the librymol).
  extern PersistPtr GetFromLibrary(const std::string molName, PersistPtr ppMolList);
  //global variable
  extern std::string libfile;

  class CalcMethod;

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
    bool calculate(size_t nCond, map<string, double> &phenRates) {

      //
      // Reset microcanonical rate re-calculation flag as parameters, such
      // as reaction threshold may have been altered between invocations of
      // this method.
      //
      for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
        (*m_pReactionManager)[i]->resetCalcFlag();
      }

      double temp = m_pConditionsManager->PTPointTemp(nCond) ; 
      double conc = m_pConditionsManager->PTPointConc(nCond);
      Precision precision = m_pConditionsManager->PTPointPrecision(nCond) ;
      string bathGas(m_pConditionsManager->PTPointBathGas(nCond)) ;
      double MaxT(2.0*temp) ;
      return calculate(temp, conc, precision, phenRates, MaxT, bathGas);
    } ;

    // Print system configuration
    void configuration(void);

    // Deduce the role of each molecule and add it to the XML file 
    bool assignMolTypes(PersistPtr ppIOPtr) ;

    void WriteMetadata(const std::string& infilename);

    PersistPtr getPersistPtr() { return m_ppIOPtr; }
    PersistPtr getAnalysisPtr() { return m_ppAnalysis; }

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
    double calcChiSqRawData(const unsigned calPoint, vector<double> &residuals, bool writeReport);

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

    // current <analysis> element
    PersistPtr m_ppAnalysis;

    // The method used for the main calculation
    CalcMethod* m_CalcMethod;

    const char* m_pTitle;
    const char* m_pDescription;

		// Collision operator instance.

    CollisionOperator m_collisionOperator ;

  } ;

}//namespace
#endif // GUARD_System_h
