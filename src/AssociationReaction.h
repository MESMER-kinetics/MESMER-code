#ifndef GUARD_AssociationReaction_h
#define GUARD_AssociationReaction_h

//-------------------------------------------------------------------------------------------
//
// AssociationReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the AssociationReaction class.
//
// This class describes a linearized association reaction in which one reactant is in such
// excess that reaction does not significantly alter its concentration. The reactant with
// the smaller concentration is deemed to be the pseudo-isomer of the reaction. Following
// regular isomerization, a number of reaction properties are delegated to the pseudo-isomer,
// e.g. the zero point energy location of the associating pair. Other quantities, such as
// the combined density of states, are properties of the reaction and are held at that level.
//
//-------------------------------------------------------------------------------------------
#include "Reaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{

  class AssociationReaction : public Reaction
  {
  public:

    // Constructors.
    AssociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id, bool isReactant)
      :Reaction(pMoleculeManager, Env, id),
      m_sourceMap(NULL),
      m_rct1(NULL),
      m_rct2(NULL),
      m_pdt1(NULL),
      deficientReactantLocation(isReactant),
      m_ERConc(0.),
      m_rctsCellEne(),
      m_rctsCellDOS(),
      m_rctsGrainEne(),
      m_rctsGrainDOS() 
    {}

    // Destructor.
    virtual ~AssociationReaction(){
      delete m_sourceMap;
    };

    void putSourceMap(sourceMap *sourcemap){m_sourceMap = sourcemap ; } ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_pdt1) ;
      return 1;
    } ;

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // Get the principal source reactant (i.e. reactant not in excess).
    virtual ModelledMolecule *get_pseudoIsomer(void) const {return m_rct1 ; } ;
    virtual ModelledMolecule *get_excessReactant(void) const {return m_rct2 ; } ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const { return m_rct1->get_zpe() + m_rct2->get_zpe() - getEnv().EMin; }
    virtual double get_relative_pdtZPE() const { return m_pdt1->get_zpe() - getEnv().EMin; }
    virtual double get_relative_TSZPE(void) const { return m_TransitionState->get_zpe() - getEnv().EMin; };

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

	// Is reaction equilibrating and therefore contributes
	// to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double &Keq, ModelledMolecule **rct, ModelledMolecule **pdt) ;

    // Get reactants cell density of states.
    virtual void getRctsCellDensityOfStates(std::vector<double> &cellDOS) ;

    // Set reactants cell  density of states.
    virtual void setRctsCellDensityOfStates(std::vector<double> &cellDOS) { m_rctsCellDOS = cellDOS ; } ;

    // Get reactants cell energies.
    virtual void getRctsCellEnergies(std::vector<double> &CellEne) ;

    // Set reactants cell energies.
    virtual void setRctsCellEnergies(std::vector<double> &CellEne) { m_rctsCellEne = CellEne ; } ;

    // Get reactants grain density of states.
    virtual void getRctsGrainDensityOfStates(std::vector<double> &grainDOS) ;

    // Get reactants grain energies.
    virtual void getRctsGrainEnergies(std::vector<double> &grainEne) ;

    virtual DensityOfStatesCalculator* get_rctsDensityOfStatesCalculator(){return get_pseudoIsomer()->get_DensityOfStatesCalculator(); }

    bool calcRctsDensityOfStates();

    double rctsRovibronicGrnCanPrtnFn();

  private:

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

    // Read parameters requires to determine reaction heats and rates.
    virtual bool ReadRateCoeffParameters(PersistPtr ppReac) ;

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    // Test k(T)
    virtual void testRateConstant();

    sourceMap *m_sourceMap ;

    // Reaction composition:

    ModelledMolecule    *m_rct1 ;   // Reactant Molecule.
    ModelledMolecule    *m_rct2 ;   // Subsidiary reactant molecule.
    CollidingMolecule   *m_pdt1 ;   // Product Molecule.

    bool deficientReactantLocation; // true if 1st rct in XML file is deficient false if 2nd reactant is deficient
    double               m_ERConc ; // Concentration of the excess reactant

    //
    // Convoluted cell and grain averages for m_rct1 and m_rct2.
    //
    std::vector<double> m_rctsCellEne ;   // Cell energies of reactants.                        
    std::vector<double> m_rctsCellDOS ;   // Convoluted cell density of states of reactants.           
    std::vector<double> m_rctsGrainEne ;  // Grain average energy array.
    std::vector<double> m_rctsGrainDOS ;  // Grain density of states array.

  } ;


}//namespace
#endif // GUARD_AssociationReaction_h
