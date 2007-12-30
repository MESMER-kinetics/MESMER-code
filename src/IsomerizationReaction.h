#ifndef GUARD_IsomerizationReaction_h
#define GUARD_IsomerizationReaction_h

//-------------------------------------------------------------------------------------------
//
// IsomerizationReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the IsomerizationReaction class.
//
//-------------------------------------------------------------------------------------------

#include "MoleculeManager.h"
#include "MicroRate.h"
#include "Reaction.h"

namespace mesmer
{

  class IsomerizationReaction : public Reaction
  {
  public:

    // Constructors.
	  IsomerizationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
		Reaction(pMoleculeManager, Env, id)
		{m_reactiontype = ISOMERIZATION; } ;

    // Destructor.
		~IsomerizationReaction() {} ;

	virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

	private:
 } ;


}//namespace
#endif // GUARD_IsomerizationReaction_h
