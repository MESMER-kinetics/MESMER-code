#ifndef GUARD_Reaction_h
#define GUARD_Reaction_h

//-------------------------------------------------------------------------------------------
//
// Reaction.h
//
// Author: Struan Robertson
// Date:   1/Feb/2003
//
// This header file contains the declaration of the Reaction class.
//
//-------------------------------------------------------------------------------------------

#include "MoleculeManager.h"
#include "MicroRate.h"

namespace mesmer
{
  class Reaction
  {

  public:

    // Type of reaction.
    typedef enum ReactionType{ASSOCIATION,
                              DISSOCIATION,
                              ISOMERIZATION,
                              EXCHANGE,
                              ERROR_REACTION } ;

    typedef std::map<CollidingMolecule*, int> isomerMap ;
    typedef std::map<SuperMolecule*    , int> sourceMap ;

  protected:
    const MesmerEnv&   m_Env;
    std::string        m_Name ;             // Reaction name.
	ReactionType       m_reactiontype ;     // Type of reaction.
    MoleculeManager   *m_pMoleculeManager ; // Pointer to molecule manager.

    //
    // Reaction composition.
    //
    SuperMolecule     *m_srct ;              // Reactant molecules as a super-reactant
    CollidingMolecule *m_rct1 ;              // Reactant Molecule.
    ModelledMolecule  *m_rct2 ;              // Subsidiary reactant molecule.
    CollidingMolecule *m_pdt1 ;              // Product Molecule.
    ModelledMolecule  *m_pdt2 ;              // Subsidiary product molecule.
    TransitionState   *m_TransitionState;    // TransitionState

//  private:
    //
    // Reaction Rate data.
    //
    double              m_HeatOfReaction ;   // The heat of reaction corrected for zero point energies.
    double              m_kfwd ;             // Forward canonical (high pressure) rate coefficient.
    std::vector<double> m_CellKfmc ;         // Forward microcanonical rate coefficients.
    std::vector<double> m_GrainKfmc ;        // Grained averaged forward microcanonical rates.

    // Calculate reaction equilibrium constant.
    double calcEquilibriumConstant() ;

  private:
    double              m_ActEne ;           // Activation Energy
    double              m_PreExp ;           // Preexponetial factor
    double              m_NInf ;             // Modified Arrhenius parameter
    double              m_ERConc ;           // Concentration of the excess reactant

    // I/O and control
    PersistPtr          m_ppPersist;         // Conduit for I/O

    // Point to microcanoical rate coeff. calculator.
    MicroRateCalculator *m_pMicroRateCalculator ;

    // Read a molecule name from the XML file and look it up
    Molecule* GetMolRef(PersistPtr pp);

    // Grain average microcanonical rate coefficients.
    bool grnAvrgMicroRateCoeffs();

    // Wrapper function to calculate and grain average microcanoincal rate coeffcients.
    bool calcGrnAvrgMicroRateCoeffs() ;

    // Add reaction terms to collision matrix.
	virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) = 0 ;

  public:

    // Constructors.
//    Reaction();

    Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id);

    // Destructor.
    ~Reaction();

    // Copy constructor.
    //   Reaction(const Reaction& reaction) ;

    // Assignment operator.
    //   Reaction& operator=(const Reaction& reaction) ;

    // Initialize reaction.
    bool InitializeReaction(PersistPtr ppReac) ;

    std::string& getName()          { return m_Name ; } ;
    const MesmerEnv& getEnv() const { return m_Env; } ;

    // Add microcanonical terms to collision operator
    void AddMicroRates(dMatrix *CollOptr,
                       isomerMap &isomermap,
                       const double rMeanOmega);

    // Access microcanoincal rate coeffcients.
    void get_MicroRateCoeffs(std::vector<double> &kmc) ;

    double get_PreExp() const                   { return m_PreExp ; } ;
    double get_ActivationEnergy()const          { return m_ActEne ; } ;
    double get_NInf()const                      { return m_NInf   ; } ;

    TransitionState* get_TransitionState()const { return m_TransitionState ; } ;

    // Get unimolecualr species information:
    virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const = 0 ;

    // Product information:
    SuperMolecule* get_bi_molecularspecies(void) const;

    // Get the principal source reactant (i.e. reactant not in excess) if it exists.
    ModelledMolecule *get_pseudoIsomer() const ;

  } ;


}//namespace
#endif // GUARD_Reaction_h
