#ifndef GUARD_IrreversibleExchangeReaction_h
#define GUARD_IrreversibleExchangeReaction_h

//-------------------------------------------------------------------------------------------
//
// IrreversibleExchangeReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the IrreversibleExchangeReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

namespace mesmer
{

  class IrreversibleExchangeReaction : public Reaction
  {
  public:

    // Constructors.
    IrreversibleExchangeReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant)
      : Reaction(pMoleculeManager, Env, Flags, id),
      m_sourceMap(NULL),
      m_rct1(NULL),
      m_rct2(NULL), 
      m_pdt1(NULL), 
      m_pdt2(NULL),
      deficientReactantLocation(isReactant)
      {} 

    // Destructor.
    virtual ~IrreversibleExchangeReaction(){}

    void putSourceMap(sourceMap *sourcemap){m_sourceMap = sourcemap ; } ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const 
    {unimolecularspecies.push_back(m_rct1); return 1;} ;

    // Get the principal source reactant (i.e. reactant not in excess).
    virtual Molecule *get_pseudoIsomer(void) const {return m_rct1 ; } ;
    virtual Molecule *get_excessReactant(void) const {return m_rct2 ; } ;

    // Return products
    virtual int get_products(std::vector<Molecule *> &product) const
    {
      product.push_back(m_pdt1) ;
      if(m_pdt2){
        product.push_back(m_pdt2) ;
        return 2;
      }
      return 1;
    } ;

    // return the colloptrsize of the reactants
    virtual int getRctColloptrsize(){return 1;}

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // Get reactants grain ZPE
    const int get_rctsGrnZPE(void);

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const {return m_rct1->g_dos->get_zpe() + m_rct2->g_dos->get_zpe() - getEnv().EMin;}
    virtual double get_relative_pdtZPE() const {return m_pdt1->g_dos->get_zpe() + m_pdt2->g_dos->get_zpe() - getEnv().EMin;}
    virtual double get_relative_TSZPE(void) const {return m_TransitionState->g_dos->get_zpe() - getEnv().EMin;};
    
    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

    // is the reaction an irreversible reaction
    virtual bool isIrreversible(){return true;};

    // is reaction unimolecular
    virtual bool isUnimolecular(){return false;};

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual Molecule *get_reactant(void) const {return m_rct1;};

    // calculate the effective threshold energy for utilizing in k(E) calculations, necessary for cases
    // with a negative threshold energy
    void calcEffGrnThresholds(void);
    
    // calculate reactant/product DOS/partition function
    bool calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne);
    bool calcPdtsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne);
    void getRctsCellDensityOfStates(vector<double> &cellDOS);
    void getPdtsCellDensityOfStates(vector<double> &cellDOS);
    virtual DensityOfStatesCalculator* get_rctsDensityOfStatesCalculator(){return get_pseudoIsomer()->g_dos->get_DensityOfStatesCalculator(); }
    virtual DensityOfStatesCalculator* get_pdtsDensityOfStatesCalculator(){return m_pdt1->g_dos->get_DensityOfStatesCalculator(); }
    virtual double pdtsRovibronicGrnCanPrtnFn();
    virtual double rctsRovibronicGrnCanPrtnFn();

  private:

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();

    // Test k(T)
    virtual void testRateConstant();

    sourceMap *m_sourceMap ;

    Molecule    *m_rct1 ;                 // Reactant Molecule.
    Molecule    *m_rct2 ;                 // Subsidiary reactant molecule. 
    Molecule    *m_pdt1 ;                 // Product Molecule.
    Molecule    *m_pdt2 ;                 // Subsidiary product molecule.

    bool deficientReactantLocation; // true if 1st rct in XML file is deficient false if 2nd reactant is deficient
  } ;


}//namespace
#endif // GUARD_IrreversibleExchangeReaction_h
