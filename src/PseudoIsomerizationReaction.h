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
#include "PriorDistFragmentation.h"

using namespace Constants ;
using namespace mesmer;

namespace mesmer
{

  class PseudoIsomerizationReaction : public AssociationReaction
  {
  public:

    // Constructors.
    PseudoIsomerizationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char *id, bool isReactant)
      :AssociationReaction(pMoleculeManager, Env, Flags, id, isReactant), m_fragDist(NULL) {} ;

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
      PersistPtr ppDstbn = ppReac->XmlMoveTo("me:FragmentDist");
      if (ppDstbn) {
        const char* ptxt = ppDstbn->XmlReadValue("xsi:type",optional); ;
        string fragMod(ptxt) ;
        if (fragMod == "me:modPrior") 
          m_fragDist = new modPriorDist(ppDstbn, this->getName()) ;
        else 
          m_fragDist = new priorDist() ;
      } else {
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
