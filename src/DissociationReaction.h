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
        DissociationReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id)
            : Reaction(pMoleculeManager, Env, id),       
            m_rct1(NULL),
            m_pdt1(NULL), 
            m_pdt2(NULL) {} ;

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

        // return relative reactant, product and transition state zero-point energy
        virtual double get_relative_rctZPE() const {return m_rct1->get_zpe() - getEnv().EMin;}
        virtual double get_relative_pdtZPE() const {return m_pdt1->get_zpe() + m_pdt2->get_zpe() - getEnv().EMin;}
        virtual double get_relative_TSZPE(void) const {return m_TransitionState->get_zpe() - getEnv().EMin;};

    private:

        // Add reaction terms to collision matrix.
        virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

        // Calculate reaction equilibrium constant.
        virtual double calcEquilibriumConstant() ;

        // Read parameters requires to determine reaction heats and rates.
        virtual bool ReadRateCoeffParameters(PersistPtr ppReac) ;

        // Grain averaged microcanonical rate coefficients.
        virtual void calcGrainRateCoeffs();

        CollidingMolecule   *m_rct1 ;                 // Reactant Molecule.
        ModelledMolecule    *m_pdt1 ;                 // Product Molecule.
        ModelledMolecule    *m_pdt2 ;                 // Subsidiary product molecule.

    } ;


}//namespace
#endif // GUARD_DissociationReaction_h
