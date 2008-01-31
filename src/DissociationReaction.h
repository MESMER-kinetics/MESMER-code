#ifndef GUARD_DissociationReaction_h
#define GUARD_DissociationReaction_h

//-------------------------------------------------------------------------------------------
//
// DissociationReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the DissociationReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

namespace mesmer
{

  class DissociationReaction : public Reaction
  {
  public:

    // Constructors.
    DissociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
      Reaction(pMoleculeManager, Env, id){} ;

      // Destructor.
      virtual ~DissociationReaction(){} ;

      // Get unimolecular species information:
      virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const
      {
        unimolecularspecies.push_back(m_rct1) ;
        return 1;
      } ;

  private:

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

  } ;


}//namespace
#endif // GUARD_DissociationReaction_h
