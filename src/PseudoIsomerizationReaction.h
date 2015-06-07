#ifndef GUARD_PseudoIsomerizationReaction_h
#define GUARD_PseudoIsomerizationReaction_h

//-------------------------------------------------------------------------------------------
//
// PseudoIsomerizationReaction.h
//
// Author: Struan Robertson
// Date:   26/May/2013
//
// This header file contains the declaration of the PseudoIsomerizationReaction class.
//
// This class describes a linearized association reaction in which one reactant is in such
// excess that reaction does not significantly alter its concentration. The reactant with
// the smaller concentration is deemed to be the pseudo-isomer of the reaction. Following
// regular isomerization, a number of reaction properties are delegated to the pseudo-isomer,
// e.g. the zero point energy location of the associating pair. Other quantities, such as
// the combined density of states, are properties of the reaction and are held at that level.
//
//-------------------------------------------------------------------------------------------
#include "AssociationReaction.h"

using namespace Constants ;
using namespace mesmer;

namespace mesmer
{

  // Forward declatation of PseudoIsomerizationReaction class.

  class PseudoIsomerizationReaction ;

  //
  // Abstract base class for the calculation of fragment distribution on dissocistion.
  // SHR 16/Jun/2013: This class definition should probably the base class to a set of
  // plug-in classes.
  //
  class FragDist 
  {
  public: 

	// Constructors.
	FragDist(){} ;

    // Destructor.
    virtual ~FragDist(){} ;

	// Initialize the fragment distribution.
	virtual void initialize(PseudoIsomerizationReaction *pReaction) = 0 ;

	// Calculate distribution.
	virtual void calculate(double excessEnergy, std::vector<double>& dist, size_t size) = 0 ;

	// Return resources
	virtual void clear() = 0 ;
	
  } ;

  //
  // Implementation class for the calculation of fragment distribution on dissocistion
  // for the prior model.
  // SHR 16/Jun/2013: This class definition should probably be a plug-in class.
  //
  class priorDist : public FragDist 
  {
  public: 

	// Constructors.
	priorDist(){} ;

    // Destructor.
    virtual ~priorDist(){} ;

	// Initialize the fragment distribution.
	virtual void initialize(PseudoIsomerizationReaction *pReaction) ;

	// Calculate distribution
	virtual void calculate(double excessEnergy, std::vector<double>& dist, size_t size) ;

	// Return resources
	virtual void clear() {
	  m_rctDOS.clear() ;
	  m_upperConv.clear() ;
	  m_lowerConv.clear() ;
	};

  protected: 

	PseudoIsomerizationReaction *m_pReaction ;

	vector<double> m_rctDOS;

	vector<double> m_upperConv;
    
	vector<double> m_lowerConv;

  } ;

  class modPriorDist : public priorDist 
  {
  public: 

	// Constructors.
	modPriorDist(){} ;

    // Destructor.
    virtual ~modPriorDist(){} ;

	// Initialize the fragment distribution.
	virtual void initialize(PseudoIsomerizationReaction *pReaction) ;

  } ;


  class PseudoIsomerizationReaction : public AssociationReaction
  {
  public:

    // Constructors.
    PseudoIsomerizationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant)
      :AssociationReaction(pMoleculeManager, Env, Flags, id, isReactant) {} ;

    // Destructor.
    virtual ~PseudoIsomerizationReaction(){}

	virtual void updateSourceMap(molMapType& sourcemap) {/* This is NULL operation as source is treated as an isomer. */ } ;

	  // Get unimolecular species information:
    virtual int get_unimolecularspecies(std::vector<Molecule *> &unimolecularspecies) const
    {
      unimolecularspecies.push_back(m_pdt1) ;
      unimolecularspecies.push_back(m_rct1) ;
      return unimolecularspecies.size() ;
    } ;

	// Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) {
	 
	  // Determine fragmentation model.
	  PersistPtr ppDistribution = ppReac->XmlMoveTo("me:FragmentDist");
	  if (!ppDistribution)
		 m_fragDist = new priorDist() ;
	  else {
		const char* txt = ppReac->XmlReadValue("me:FragmentDist") ;
		string fragMod(txt) ;
		if (fragMod == "modPrior") 
		  m_fragDist = new modPriorDist() ;
		else 
		  m_fragDist = new priorDist() ;
	  }
	  return AssociationReaction::InitializeReaction(ppReac) ;
	};

	// returns the reaction type
	virtual ReactionType getReactionType(){return PSEUDOISOMERIZATION;};

	// Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) ;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) {
	  throw std::runtime_error("Contracted basis Set not yet implemeneted for pseudoisomerization.");
	};

  private:

	FragDist *m_fragDist ;

  } ;

}//namespace
#endif // GUARD_PseudoIsomerizationReaction_h
