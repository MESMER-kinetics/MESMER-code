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

#include "MoleculeManager.h"
#include "MicroRate.h"
#include "Reaction.h"

namespace mesmer
{

  class DissociationReaction : public Reaction
  {
  public:

    // Constructors.
	  DissociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
		Reaction(pMoleculeManager, Env, id)
		{m_reactiontype = DISSOCIATION ; } ;

    // Destructor.
		~DissociationReaction(){} ;

	virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

	private:

  } ;


}//namespace
#endif // GUARD_DissociationReaction_h
