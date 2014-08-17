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
      { m_UsesProductProperties = false; } 

    // Destructor.
    virtual ~IrreversibleExchangeReaction(){}

    void updateSourceMap(molMapType& sourcemap){
      if (m_rct1 && sourcemap.find(m_rct1) == sourcemap.end()){ // Reaction includes a new pseudoisomer.
        sourcemap[m_rct1] = 0 ;
      }
      m_sourceMap = &sourcemap ; 
    } ;

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const 
    {return 0;} ;

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

    virtual int get_reactants(std::vector<Molecule *> &reactants) const
    {
      reactants.push_back(m_rct1);
      reactants.push_back(m_rct2);
      return 2;
    } ;

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // Get reactants grain ZPE
    const int get_rctsGrnZPE(void);

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const {return m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin;}
    virtual double get_relative_pdtZPE() const {return m_pdt1->getDOS().get_zpe() + m_pdt2->getDOS().get_zpe() - getEnv().EMin;}
    virtual double get_relative_TSZPE(void) const {return m_TransitionState->getDOS().get_zpe() - getEnv().EMin;};
    
	// Returns header needed by reaction rate test method.
	virtual std::string TestRateCoeffHeader() const {return string("   T(K) kf/cm3mlc-1s-1 kb/cm3mlc-1s-1            Keq") ; } ; 

	// Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

    // returns the reaction type
    virtual ReactionType getReactionType(){return IRREVERSIBLE_EXCHANGE;};

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

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
		virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap){} ;

  private:

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs();
    virtual void calcFluxFirstNonZeroIdx(void);

    // Test k(T)
    virtual void testRateConstant();
    double pdtsRovibronicGrnCanPrtnFn();
    double rctsRovibronicGrnCanPrtnFn();

    molMapType *m_sourceMap ;

    Molecule    *m_rct1 ; // Reactant Molecule.
    Molecule    *m_rct2 ; // Subsidiary reactant molecule. 
    Molecule    *m_pdt1 ; // Product Molecule.
    Molecule    *m_pdt2 ; // Subsidiary product molecule.

    bool deficientReactantLocation; // true if 1st rct in XML file is deficient false if 2nd reactant is deficient
  } ;


}//namespace
#endif // GUARD_IrreversibleExchangeReaction_h
