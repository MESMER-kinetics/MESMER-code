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

  private:

    std::string        m_Name ;             // Reaction name.
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
    ReactionType       m_reactiontype ;      // Type of reaction.

    //
    // Reaction Rate data.
    //
    double              m_HeatOfReaction ;   // The heat of reaction corrected for zero point energies.
    double              m_kfwd ;             // Forward canonical (high pressure) rate coefficient.
    std::vector<double> m_CellKfmc ;         // Forward microcanonical rate coefficients.
    std::vector<double> m_GrainKfmc ;        // Grained averaged forward microcanonical rates.

    double              m_ActEne ;           // Activation Energy
    double              m_PreExp ;           // Preexponetial factor

    // I/O and control
    PersistPtr          m_ppPersist;         // Conduit for I/O

    //
    // Point to microcanoical rate coeff. calculator.
    //
    MicroRateCalculator *m_pMicroRateCalculator ;

    // Read a molecule name from the XML file and look it up
    Molecule* GetMolRef(PersistPtr pp);

    // Calculate reaction equilibrium constant.
    double calcEquilibriumConstant(const MesmerEnv &mEnv) ;

    // Grain average microcanonical rate coefficients.
    bool grnAvrgMicroRateCoeffs(const MesmerEnv &mEnv);

    // Wrapper function to calculate and grain average microcanoincal rate coeffcients.
    bool calcGrnAvrgMicroRateCoeffs(const MesmerEnv &mEnv) ;

    // Add isomer reaction terms to collision matrix.
    void AddIsomerReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega, const MesmerEnv &mEnv) ;

    // Add (reversible) association reaction terms to collision matrix.
    void AddAssocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, sourceMap &sourcemap, const double rMeanOmega, const MesmerEnv &mEnv) ;

    // Add dissociation reaction terms to collision matrix.
    void AddDissocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

  public:

    // Constructors.
    Reaction();

    Reaction(MoleculeManager *pMoleculeManager);

    // Destructor.
    ~Reaction();

    // Copy constructor.
    //   Reaction(const Reaction& reaction) ;

    // Assignment operator.
    //   Reaction& operator=(const Reaction& reaction) ;

    // Initialize reaction.
    bool InitializeReaction(PersistPtr ppReac) ;

    std::string& getName() { return m_Name ; } ;

    // Modifier for reaction type.
    void put_Reactiontype(ReactionType reactiontype) ;

    // Accessor for reaction type.
    ReactionType get_Reactiontype() const {return m_reactiontype ; } ;

    // Add microcanonical terms to collision operator
    void AddMicroRates(dMatrix *CollOptr,
                       isomerMap &isomermap,
                       sourceMap &sourcemap,
                       const double rMeanOmega,
                       const MesmerEnv &mEnv);

    // Determine the equilibrium constant.
    void CalcEquilConst() { } ;

    // Access microcanoincal rate coeffcients.
    void get_MicroRateCoeffs(std::vector<double> &kmc, const MesmerEnv &mEnv) ;

    double get_PreExp() const                   { return m_PreExp ; } ;
    double get_ActivationEnergy()const          { return m_ActEne ; } ;

    TransitionState* get_TransitionState()const { return m_TransitionState ; } ;

    // Reactant information:
    int get_unimolecularspecies(std::vector<CollidingMolecule *> &unimolecularspecies) const ;

    // Product information:
    int get_bi_molecularspecies(SuperMolecule* bi_mol) const;

    // Get the principal source reactant (i.e. reactant not in excess) if it exists.
    CollidingMolecule *get_pseudoIsomer() const ;

  } ;

}//namespace
#endif // GUARD_Reaction_h
