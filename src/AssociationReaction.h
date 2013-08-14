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
    AssociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant)
      :Reaction(pMoleculeManager, Env, Flags, id),
      m_sourceMap(NULL),
      m_rct1(NULL),
      m_rct2(NULL),
      m_pdt1(NULL),
      m_deficientReactantLocation(isReactant),
      m_GrainKbmc() {} ;

    // Destructor.
    virtual ~AssociationReaction(){}

    virtual void updateSourceMap(molMapType& sourcemap) {
      if (m_rct1 && sourcemap.find(m_rct1) == sourcemap.end()){ // Reaction includes a new pseudoisomer.
        sourcemap[m_rct1] = 0 ;
      }
      m_sourceMap = &sourcemap ; 
    } ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_pdt1) ;
      return 1;
    } ;

    // Get product information:
    virtual int get_products(std::vector<Molecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_pdt1) ;
      return 1;
    } ;

    // return the colloptrsize of the reactants
    virtual int getRctColloptrsize(){return 1;}

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // Get the principal source reactant (i.e. reactant not in excess).
    virtual Molecule *get_pseudoIsomer(void) const {return m_rct1 ; } ;
    virtual Molecule *get_reactant(void) const {return m_rct1;};
    virtual Molecule *get_excessReactant(void) const {return m_rct2 ; } ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const { return m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin; }
    virtual double get_relative_pdtZPE() const { return m_pdt1->getDOS().get_zpe() - getEnv().EMin; }
    virtual double get_relative_TSZPE(void) const { return m_TransitionState->getDOS().get_zpe() - getEnv().EMin; };

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

    // Is reaction equilibrating and therefore contributes
    // to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double &Keq, Molecule **rct, Molecule **pdt) ;

    // returns the reaction type
    virtual ReactionType getReactionType(){return ASSOCIATION;};

    // Get reactants cell density of states.
    void getRctsCellDensityOfStates(std::vector<double> &cellDOS) ;

	// Get Cell Boltzmann distribution for convoluted reactant DOS
	void normalizedReactantCellBoltzmannDistribution(vector<double> &CellFrac, const int totalCellNumber);

    // Get reactants grain ZPE
    const int get_rctsGrnZPE(void);

    // calculate the effective threshold energy for utilizing in k(E) calculations, necessary for cases
    // with a negative threshold energy
    void calcEffGrnThresholds(void);

    // Get cell offset for the reactants
    int get_cellOffset(void) {
      double modulus = fmod(m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin, getEnv().GrainSize);
      return int(modulus) ;
    } ;

    bool calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne);

    // Calculate rovibronic canonical partition function in the grain level for product or reactant
    virtual double rctsRovibronicGrnCanPrtnFn();
    virtual double pdtsRovibronicGrnCanPrtnFn();

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) ;

  protected:

    // Reaction composition:

    Molecule *m_rct1 ;   // Reactant Molecule.
    Molecule *m_rct2 ;   // Subsidiary reactant molecule.
    Molecule *m_pdt1 ;   // Product Molecule.

    molMapType *m_sourceMap ;

  private:

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    // Test k(T)
    virtual void testRateConstant();

    bool m_deficientReactantLocation; // true if 1st rct in XML file is deficient false if 2nd reactant is deficient

    std::vector<double>  m_GrainKbmc ;           // Grained averaged backward microcanonical rates.

  } ;


}//namespace
#endif // GUARD_AssociationReaction_h
