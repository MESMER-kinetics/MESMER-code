#ifndef GUARD_DiffusiveLossReaction_h
#define GUARD_DiffusiveLossReaction_h

//-------------------------------------------------------------------------------------------
//
// DiffusiveLossReaction.h
//
// Author: Struan Robertson
// Date:   2/May/2016
//
// This header file contains the declaration of the DiffusiveLossReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"
#include "gDensityOfStates.h"

namespace mesmer
{

  class DiffusiveLossReaction : public Reaction
  {
  public:

    // Constructors.
    DiffusiveLossReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id)
      : Reaction(pMoleculeManager, Env, Flags, id),
      m_rct1(NULL)
    { m_UsesProductProperties = false; }  

    // Destructor.
    virtual ~DiffusiveLossReaction(){} ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_rct1) ;
      return unimolecularspecies.size() ;
    } ;

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // Return relative reactant, product and transition state zero-point energy (Dummy methods for diffusion).
		virtual double get_relative_rctZPE(void) const { return m_rct1->getDOS().get_zpe() - getEnv().EMin; };
		virtual double get_relative_pdtZPE(void) const { return m_rct1->getDOS().get_zpe() - getEnv().EMin; };
		virtual double get_relative_TSZPE(void)  const { return m_rct1->getDOS().get_zpe() - getEnv().EMin; };

    // Calculate high pressure rate coefficients at current T.
    virtual void HighPresRateCoeffs(vector<double> *pCoeffs) ;

	  // Calculate reaction equilibrium constant (Dummy method for diffusion).
		virtual double calcEquilibriumConstant() { 
			double Keq(0.0);  
			return Keq; 
		};

    // Return products
    virtual int get_products(std::vector<Molecule *> &product) const
    {
      product.push_back(m_rct1) ;
      return product.size() ;
    } ;

    virtual int get_reactants(std::vector<Molecule *> &reactants) const
    {
      reactants.push_back(m_rct1);
      return reactants.size() ;
    } ;

    // returns the reaction type
    virtual ReactionType getReactionType(){
      return DIFFUSION ;
    };

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual Molecule *get_reactant(void) const {return m_rct1;};

    // calculate the effective threshold energy for utilizing in k(E) calculations, necessary for cases
    // with a negative threshold energy
    void calcEffGrnThresholds(void) {
			double RxnHeat(0.0);
			set_EffGrnFwdThreshold(int(RxnHeat)) ;
		}

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) {
			throw std::runtime_error("Contracted basis Set not yet implemeneted for difffusive loss.");
		};

  private:

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    void calcFluxFirstNonZeroIdx(void);

    Molecule    *m_rct1 ; // Reactant Molecule.

  } ;


}//namespace
#endif // GUARD_DiffusiveLossReaction_h
