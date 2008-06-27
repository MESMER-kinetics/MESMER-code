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
    IsomerizationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id)
      : Reaction(pMoleculeManager, Env, id),
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

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const {return m_rct1->get_zpe() - getEnv().EMin;}
    virtual double get_relative_pdtZPE() const {return m_pdt1->get_zpe() - getEnv().EMin;}
    virtual double get_relative_TSZPE(void) const {return m_TransitionState->get_zpe() - getEnv().EMin;};

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

  private:

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

    // Read parameters requires to determine reaction heats and rates.
    virtual bool ReadRateCoeffParameters(PersistPtr ppReac) ;

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    // Test k(T)
    virtual void testRateConstant();

    CollidingMolecule   *m_rct1 ;                 // Reactant Molecule.
    CollidingMolecule   *m_pdt1 ;                 // Product Molecule.

  } ;

}//namespace
#endif // GUARD_IsomerizationReaction_h
