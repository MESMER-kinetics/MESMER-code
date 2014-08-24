#ifndef GUARD_SecondOrderAssocReaction_h
#define GUARD_SecondOrderAssocReaction_h

//-------------------------------------------------------------------------------------------
//
// SecondOrderAssocReaction.h
//
// Author: Struan Robertson
// Date:   24/Aug/2014
//
// This header file contains the declaration of the SecondOrderAssocReaction class.
//
// This class describes a linearized second order association reaction. The basis of the 
// linearization is a Taylor expansion about the equilibrium point, such that the reacton 
// is the linear regime and the concentration of the fragments can be treated linearly and
// so can treated as the pseudo-isomer of the reaction. As with other association reactions,
// a number of reaction properties are delegated to the pseudo-isomer, e.g. the zero point
// energy location of the associating pair. Other quantities, such as the combined density
// of states, are properties of the reaction and are held at that level.
//
//-------------------------------------------------------------------------------------------
#include "AssociationReaction.h"

using namespace Constants ;
using namespace mesmer;

namespace mesmer
{

  class SecondOrderAssocReaction : public AssociationReaction
  {
  public:

    // Constructors.
    SecondOrderAssocReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant)
      :AssociationReaction(pMoleculeManager, Env, Flags, id, isReactant) {} ;

    // Destructor.
    virtual ~SecondOrderAssocReaction(){}

    virtual void updateSourceMap(molMapType& sourcemap) {
      if (m_rct1 && sourcemap.find(m_rct1) == sourcemap.end()){ // Reaction includes a new pseudoisomer.
        sourcemap[m_rct1] = 0 ;
      }
      m_sourceMap = &sourcemap ; 
    } ;

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

    // Is reaction equilibrating and therefore contributes
    // to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double &Keq, Molecule **rct, Molecule **pdt) ;

    // returns the reaction type
    virtual ReactionType getReactionType(){return SECONDORDERASSOCIATION;};

    // Get reactants cell density of states.
    void getRctsCellDensityOfStates(std::vector<double> &cellDOS) ;

    // Get reactants grain ZPE
    const int get_rctsGrnZPE(void);

    // Calculate the effective threshold energy for utilizing in k(E)
    // calculations, necessary for cases with a negative threshold energy.
    void calcEffGrnThresholds(void);

    // Get cell offset for the reactants.
    size_t get_cellOffset(void) {
      double modulus = fmod(m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin, double(getEnv().GrainSize))/getEnv().CellSize ;
      return size_t(modulus) ;
    } ;

    bool calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne);

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) ;

  } ;

}//namespace
#endif // GUARD_SecondOrderAssocReaction_h
