#ifndef GUARD_AssociationReaction_h
#define GUARD_AssociationReaction_h

//-------------------------------------------------------------------------------------------
//
// AssociationReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the AssociationReaction class.
//
//-------------------------------------------------------------------------------------------

#include "MoleculeManager.h"
#include "MicroRate.h"
#include "Reaction.h"

namespace mesmer
{

	class AssociationReaction : public Reaction
	{
	public:

		// Constructors.
		AssociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
		  Reaction(pMoleculeManager, Env, id),
			  m_sourceMap(NULL) {m_reactiontype = ASSOCIATION ; } ;

		  // Destructor.
		  ~AssociationReaction();

		  void putSourceMap(sourceMap *sourcemap){m_sourceMap = sourcemap ; } ;

		  // Get unimolecular species information:
		  virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const 
		  {        
			  unimolecularspecies.push_back(m_pdt1) ;
			  return 1;
		  } ;

	private:

		// Add reaction terms to collision matrix.
		  virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

		  sourceMap *m_sourceMap ; 

	} ;


}//namespace
#endif // GUARD_AssociationReaction_h
