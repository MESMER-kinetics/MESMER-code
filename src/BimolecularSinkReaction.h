#ifndef GUARD_BimolecularSinkReaction_h
#define GUARD_BimolecularSinkReaction_h

//-------------------------------------------------------------------------------------------
//
// BimolecularSinkReaction.h
//
// Author: david glowacki
// Date:   12 Dec 2011
//
// This header file contains the declaration of the BimolecularSinkReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{

  class BimolecularSinkReaction : public Reaction
  {
  public:
		    // Constructors.
    BimolecularSinkReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant)
      : Reaction(pMoleculeManager, Env, Flags, id),
      m_rct1(NULL),
      m_rct2(NULL),
      m_pdt1(NULL),
      m_pdt2(NULL),
      m_pdt3(NULL),
      excessReactantLocation(isReactant)
    { m_UsesProductProperties = false; }  

    // Destructor.
    virtual ~BimolecularSinkReaction(){} ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_rct1) ;
      return 1;
    } ;

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const { return m_rct1->getDOS().get_zpe() - getEnv().EMin; }   // copied from Association Reaction
    virtual double get_relative_pdtZPE() const {return m_pdt1->getDOS().get_zpe() + m_pdt2->getDOS().get_zpe() - getEnv().EMin;} // copied from Irreversible Exchange Reaction
    virtual double get_relative_TSZPE(void) const { return m_TransitionState->getDOS().get_zpe() - getEnv().EMin; };  // copied from Association Reaction

    // Return products - copied from IrreversibleUnimolecular rxn
    virtual int get_products(std::vector<Molecule *> &product) const
    {
      product.push_back(m_pdt1) ;
      if(m_pdt2){
        product.push_back(m_pdt2) ;
        if(m_pdt3){
          product.push_back(m_pdt3) ;
          return 3;
        }
        return 2;
      }
      return 1;
    } ;

    virtual int get_reactants(std::vector<Molecule *> &reactants) const
    {
      reactants.push_back(m_rct1);
      reactants.push_back(m_rct2);
      return 2;
    } ;

    const int get_pdtsGrnZPE();

	// Returns header needed by reaction rate test method.
	virtual std::string TestRateCoeffHeader() const {return string("   T(K)         kf/s-1 kb/cm3mlc-1s-1   Keq/mlc cm-3") ; } ; 

	// Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

    // calculate the effective threshold energy for utilizing in k(E) calculations, necessary for cases
    // with a negative threshold energy
    void calcEffGrnThresholds(void);

    // returns the reaction type
    virtual ReactionType getReactionType(){
      return BIMOLECULAR_SINK;
    };

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual Molecule *get_reactant(void) const {return m_rct1;};

    virtual double rctsRovibronicGrnCanPrtnFn();
    virtual double pdtsRovibronicGrnCanPrtnFn() { return m_pdt1->getDOS().rovibronicGrnCanPrtnFn();}

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) ;

  private:

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    void calcFluxFirstNonZeroIdx(void);

    // Test k(T)
    virtual void testRateConstant();

    Molecule    *m_rct1 ;                 // Reactant Molecule.
    Molecule    *m_rct2 ;                 // Reactant Molecule.
    Molecule    *m_pdt1 ;                 // Product Molecule.
    Molecule    *m_pdt2 ;                 // Subsidiary product molecule.
    Molecule    *m_pdt3 ;                 // Subsidiary product molecule.

    bool excessReactantLocation;				// true if 1st rct in XML file is deficient false if 2nd reactant is deficient

	};

}//namespace
#endif // GUARD_BimolecularSinkReaction_h
