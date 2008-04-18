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
        Reaction(pMoleculeManager, Env, id), m_pdt2(NULL){} ;

        // Destructor.
        virtual ~DissociationReaction(){} ;

        // Get unimolecular species information:
        virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const
        {
          unimolecularspecies.push_back(m_rct1) ;
          return 1;
        } ;

        // Initialize reaction.
        virtual bool InitializeReaction(PersistPtr ppReac) ;

        virtual void grainRateCoeffDetailedBalance(const int dir) {} ;

  private:

    // Grain average microcanonical rate coefficients.
    virtual bool grnAvrgMicroRateCoeffs();

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() ;

    // Read parameters requires to determine reaction heats and rates.
    virtual bool ReadRateCoeffParameters(PersistPtr ppReac) ;

     ModelledMolecule    *m_pdt2 ;                 // Subsidiary product molecule.

  } ;


}//namespace
#endif // GUARD_DissociationReaction_h
