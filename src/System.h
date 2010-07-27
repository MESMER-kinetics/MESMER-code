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

    // Begin single calculation.
    bool calculate(double& chiSquare) ;
  
    // Print system configuration
    void configuration(void);

    // Deduce the me:type of each molecule and add it to the XML file 
    bool assignMolTypes(PersistPtr ppIOPtr) ;
    
    void WriteMetadata();

    PersistPtr getPersistPtr() { return m_ppIOPtr; }

    // Mesmer control flags.
    MesmerFlags m_Flags;

    // Build reaction operator for system.
    bool BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags) ;

    // Calculate the equilibrium fraction of each species in the system.
    bool calculateEquilibriumFractions(const double beta);

    // Diagonalize the reaction operator.
    void diagReactionOperator(const MesmerFlags &mFlags, const MesmerEnv &mEnv, const int precision, PersistPtr ppAnalysis) ;

    // Calculate the time evolution of the system
    bool timeEvolution(MesmerFlags& mFlags, PersistPtr ppPopList);

    // Calculates the Bartis-Widom macroscopic rate coefficients.
    bool BartisWidomPhenomenologicalRates(qdMatrix& rates, MesmerFlags& mFlags, PersistPtr ppBase);
    
    // Calculates the Bartis-Widom macroscopic rate coefficients, using the contracted basis set eigenvectors.
    bool BartisWidomBasisSetRates(qdMatrix& rates, MesmerFlags& mFlags);

    // Write out phenomenological rate coefficients.
    bool PrintPhenomenologicalRates(qdMatrix& Kr, qdb2D& Kp, int numberOfCemeteries, MesmerFlags& mFlags, PersistPtr ppList) ;

  private:

    typedef std::map<Reaction* , int, Reaction::ReactionPtrLess> sinkMap ;

    void readPTs(PersistPtr);

    bool ReadRange(const std::string&    name,
      std::vector<double>&  vals,
      PersistPtr            ppbase,
      bool                  MustBeThere=true);

    // sets grain parameters and determines system environment
    bool SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne);

    // Construct a transition matrix based on grains.
    void constructGrainMatrix(int msize);
    
    // Construct a transition matrix based on collision operator eigenfunctions.
    void constructBasisMatrix(void);

    void printReactionOperator(const MesmerFlags &mFlags);

    void printEigenvectors(const MesmerFlags &mFlags, std::ostream& os);

    bool produceEquilibriumVector();

    bool produceInitialPopulationVector(vector<double>& initDist);

    int getSpeciesSequenceIndex(const std::string ref);

    double calcChiSquare(const qdMatrix& mesmerRates, vector<conditionSet>& expRates, stringstream &rateCoeffTable) ;

	// This method locates all sinks and determines their location in the relevant
    // isomer or source map. 
    void locateSinks() ;

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

    // The method used for the main calculation
    CalcMethod* m_CalcMethod;

    const char* m_pTitle;

    // Maps the location of individual reactant collision operator and source terms in the reaction operator.
    Reaction::molMapType    m_isomers;
    Reaction::molMapType    m_sources;
    sinkMap                 m_sinkRxns;
    sinkMap                 m_SinkSequence;

    // Mean collision frequency.
    double                  m_meanOmega;

    // The system transition matrix and associated eigenvalues and eigenvectors.
    qdMatrix               *m_reactionOperator ;
    qdMatrix               *m_eigenvectors;
    std::vector<qd_real>    m_eigenvalues;

    // Map modelled molecules (isomers + sources) with their sequence in the transition matrix.
    Reaction::molMapType    m_SpeciesSequence;

    // Equilibrium distribution.
    std::vector<qd_real>    m_eqVector;

    bool                    m_punchSymbolGathered;
    
  } ;
}//namespace
#endif // GUARD_System_h
