#ifndef GUARD_ExchangeReaction_h
#define GUARD_ExchangeReaction_h

//-------------------------------------------------------------------------------------------
//
// ExchangeReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the ExchangeReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

namespace mesmer
{

  class ExchangeReaction : public Reaction
  {
  public:

    // Constructors.
    ExchangeReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id):
      Reaction(pMoleculeManager, Env, id)
      {m_reactiontype = EXCHANGE ; } ;

    // Destructor.
    ~ExchangeReaction(){} ;

    // Get unimolecualr species information:
    virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const 
    {return 0 ;} ;

  private:

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

  } ;


}//namespace
#endif // GUARD_ExchangeReaction_h
