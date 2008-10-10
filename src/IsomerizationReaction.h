#ifndef GUARD_IsomerizationReaction_h
#define GUARD_IsomerizationReaction_h

//-------------------------------------------------------------------------------------------
//
// IsomerizationReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the IsomerizationReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

namespace mesmer
{

  class IsomerizationReaction : public Reaction
  {
  public:

    // Constructors.
    IsomerizationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id)
      : Reaction(pMoleculeManager, Env, Flags, id),
      m_rct1(NULL),
      m_pdt1(NULL) {} ;

    // Destructor.
    virtual ~IsomerizationReaction() {} ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_rct1) ;
      unimolecularspecies.push_back(m_pdt1) ;
      return 2 ;
    } ;

    // return the colloptrsize of the reactants
    virtual int getRctColloptrsize(){return m_rct1->get_colloptrsize();}

    // Return products
    virtual int get_products(std::vector<ModelledMolecule *> &product) const
    {
      product.push_back(m_pdt1) ;
      return 1;
    } ;

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const {return m_rct1->get_zpe() - getEnv().EMin;}
    virtual double get_relative_pdtZPE() const {return m_pdt1->get_zpe() - getEnv().EMin;}
    virtual double get_relative_TSZPE(void) const {return m_TransitionState->get_zpe() - getEnv().EMin;};

    // Is reaction equilibrating and therefore contributes
    // to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double &Keq, ModelledMolecule **rct, ModelledMolecule **pdt) ;

    // calculate the effective threshold energy for utilizing in k(E) calculations, necessary for cases
    // with a negative threshold energy
    void calcEffGrnThresholds(void);

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual ModelledMolecule *get_reactant(void) const {return m_rct1;};

    // Calculate rovibronic canonical partition function in the grain level for product or reactant
    virtual double pdtsRovibronicGrnCanPrtnFn();
    virtual double rctsRovibronicGrnCanPrtnFn();


  private:

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    // Test k(T)
    virtual void testRateConstant();

    CollidingMolecule   *m_rct1 ;                 // Reactant Molecule.
    CollidingMolecule   *m_pdt1 ;                 // Product Molecule.

  } ;

}//namespace
#endif // GUARD_IsomerizationReaction_h
