#ifndef GUARD_ExchangeReaction_h
#define GUARD_ExchangeReaction_h

//-------------------------------------------------------------------------------------------
//
// ExchangeReaction.h
//
// Author: Struan Robertson
// Date:   Aug/2019
//
// This header file contains the declaration of the ExchangeReaction class. This reaction 
// allows unimolecular reactant to be prepared in a non-Boltzmann state as a consequence of
// a bimolecular reaction.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"
#include "gDensityOfStates.h"

namespace mesmer
{
  // Forward declarations.

  class FragDist;

  class ExchangeReaction : public Reaction
  {
  public:

    // Constructors.
    ExchangeReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant, bool isProduct)
      : Reaction(pMoleculeManager, Env, Flags, id),
      m_sourceMap(NULL),
      m_rct1(NULL),
      m_rct2(NULL),
      m_pdt1(NULL),
      m_pdt2(NULL),
      m_deficientReactantLocation(isReactant),
      m_modelledProductLocation(isProduct),
      m_EPConc(0.0),
      m_fragDist(NULL)
    {
      m_UsesProductProperties = false;
    }

    // Destructor.
    virtual ~ExchangeReaction() {}

    void updateSourceMap(molMapType& sourcemap) {
      if (m_rct1 && sourcemap.find(m_rct1) == sourcemap.end()) { // Reaction includes a new pseudoisomer.
        sourcemap[m_rct1] = 0;
      }
      m_sourceMap = &sourcemap;
    };

    // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const { 
      unimolecularspecies.push_back(m_pdt1);
      return unimolecularspecies.size();
    }

    // Get the principal source reactant (i.e. reactant not in excess).
    virtual Molecule *get_pseudoIsomer(void) const { return m_rct1; };
    virtual Molecule *get_excessReactant(void) const { return m_rct2; };

    // Return products
    virtual int get_products(std::vector<Molecule *> &products) const
    {
      products.push_back(m_pdt1);
      products.push_back(m_pdt2);
      return products.size();
    };

    virtual int get_reactants(std::vector<Molecule *> &reactants) const
    {
      reactants.push_back(m_rct1);
      reactants.push_back(m_rct2);
      return reactants.size();
    };

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac);

    // Get reactants grain ZPE
    const int get_rctsGrnZPE(void);

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const { return m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin; }
    virtual double get_relative_pdtZPE() const { return m_pdt1->getDOS().get_zpe() + m_pdt2->getDOS().get_zpe() - getEnv().EMin; }
    virtual double get_relative_TSZPE(void) const { return m_TransitionState->getDOS().get_zpe() - getEnv().EMin; };

    virtual double get_ThresholdEnergy(void) { return get_relative_TSZPE() - get_relative_rctZPE(); };

    // Calculate high pressure rate coefficients at current T.
    virtual void HighPresRateCoeffs(vector<double> *pCoeffs);

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant();

    // Is reaction equilibrating and therefore contributes
    // to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double& Keq, Molecule** rct, Molecule** pdt);

    // returns the reaction type
    virtual ReactionType getReactionType() { return BIMOLECULAR_EXCHANGE; };

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual Molecule *get_reactant(void) const { return m_rct1; };

    // calculate the effective threshold energy for utilizing in k(E) calculations.
    void calcEffGrnThresholds(void);

    // Wrapper function to calculate and grain average microcanoincal rate coeffcients.
    virtual bool calcGrnAvrgMicroRateCoeffs() ;

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega);

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) {};

    // The following method takes an effective unimolecular rate 
    // coefficient and normalizes it by concentration to obtain
    // a second order rate coefficient.
    virtual void normalizeRateCoefficient(double &rateCoefficient, std::string ref = "") const {
      // Check the sense of the rate coefficient by checking for the principal reactant.
      if (ref == m_rct1->getName())
        rateCoefficient /= get_concExcessReactant();
      else if (ref == m_pdt1->getName())
        rateCoefficient /= m_EPConc ;
      else
        throw(std::runtime_error("Reference, " + ref + ", is not a principal species in reaction " + getName() + "."));
    };

  private:

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs() {};
    virtual void calcFluxFirstNonZeroIdx(void);

    molMapType *m_sourceMap;

    Molecule    *m_rct1; // Reactant Molecule.
    Molecule    *m_rct2; // Subsidiary reactant molecule, in excess. 
    Molecule    *m_pdt1; // Product Molecule.
    Molecule    *m_pdt2; // Subsidiary product molecule, in excess.

    bool m_deficientReactantLocation; // True if 1st reactant in XML file is deficient false if 2nd reactant is deficient.
    bool m_modelledProductLocation;   // True if 1st product in XML file is modelled false if 2nd reactant is deficient.

    double    m_EPConc;  // Concentration of the excess product.

    FragDist* m_fragDist;
  };

}//namespace
#endif // GUARD_ExchangeReaction_h
