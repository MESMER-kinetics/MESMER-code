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
      :AssociationReaction(pMoleculeManager, Env, Flags, id, isReactant) { Flags.bIsSystemSecondOrder = true; } ;

    // Destructor.
    virtual ~SecondOrderAssocReaction(){}

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) ;

    // return relative reactant, product and transition state zero-point energy
    virtual double get_relative_rctZPE() const { return m_rct1->getDOS().get_zpe() - getEnv().EMin; }

	// Reset zero point energy locations of the reactants such that
	// location of the pair is entirely on the pseudoisomer.
	virtual double resetZPEofReactants() ;

    // returns the reaction type
    virtual ReactionType getReactionType(){return SECONDORDERASSOCIATION;};

    // Get reactants grain ZPE
    virtual const int get_rctsGrnZPE(void);

    // Get cell offset for the reactants.
    virtual size_t get_cellOffset(void) {
      double modulus = fmod(m_rct1->getDOS().get_zpe() - getEnv().EMin, double(getEnv().GrainSize))/getEnv().CellSize ;
      return size_t(modulus) ;
    } ;

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) ;

	// The following method takes an effective unimolecular rate 
	// coefficient and normalizes it by concentration to obtain
	// a second order rate coefficient.
	virtual double normalizeRateCoefficient(const double rateCoefficient) const {return rateCoefficient/(4.0*m_ERConc) ; }  ;

  } ;

}//namespace
#endif // GUARD_SecondOrderAssocReaction_h
