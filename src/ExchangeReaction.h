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
        ExchangeReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id)
            : Reaction(pMoleculeManager, Env, id),
            m_rct1(NULL),
            m_rct2(NULL), 
            m_pdt1(NULL), 
            m_pdt2(NULL) {} ;

        // Destructor.
        virtual ~ExchangeReaction(){} ;

        // Get unimolecualr species information:
        virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const 
        {return 0 ;} ;

        // Initialize reaction.
        virtual bool InitializeReaction(PersistPtr ppReac) ;

        // return relative reactant, product and transition state zero-point energy
        virtual double get_relative_rctZPE() const {return m_rct1->get_zpe() + m_rct2->get_zpe() - getEnv().EMin;}
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

        ModelledMolecule    *m_rct1 ;                 // Reactant Molecule.
        ModelledMolecule    *m_rct2 ;                 // Subsidiary reactant molecule. 
        ModelledMolecule    *m_pdt1 ;                 // Product Molecule.
        ModelledMolecule    *m_pdt2 ;                 // Subsidiary product molecule.

    } ;


}//namespace
#endif // GUARD_ExchangeReaction_h
