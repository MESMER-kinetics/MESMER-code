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
    AssociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id)
      :Reaction(pMoleculeManager, Env, id),
      m_sourceMap(NULL), 
      m_srct(NULL), 
      m_rct1(NULL),
      m_rct2(NULL), 
      m_pdt1(NULL),
	  m_ERConc(0.) {}

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

    // Product information:
    virtual SuperMolecule* get_bi_molecularspecies(void) const {return m_srct ; } ;

    // Get the principal source reactant (i.e. reactant not in excess).
    virtual ModelledMolecule *get_pseudoIsomer(void) const {return m_rct1 ; } ;
    virtual ModelledMolecule *get_excessReactant(void) const {return m_rct2 ; } ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const { return m_rct1->get_zpe() + m_rct2->get_zpe() - getEnv().EMin; }
    virtual double get_relative_pdtZPE() const { return m_pdt1->get_zpe() - getEnv().EMin; }
    virtual double get_relative_TSZPE(void) const { return m_TransitionState->get_zpe() - getEnv().EMin; };
    
    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

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

    SuperMolecule       *m_srct ;   // Reactant molecules as a super-reactant
    ModelledMolecule    *m_rct1 ;   // Reactant Molecule.
    ModelledMolecule    *m_rct2 ;   // Subsidiary reactant molecule.
    CollidingMolecule   *m_pdt1 ;   // Product Molecule.

	double               m_ERConc ; // Concentration of the excess reactant

  } ;


}//namespace
#endif // GUARD_AssociationReaction_h
