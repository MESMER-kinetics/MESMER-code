#ifndef GUARD_DissociationReaction_h
#define GUARD_DissociationReaction_h

//-------------------------------------------------------------------------------------------
//
// IrreversibleReaction.h
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This header file contains the declaration of the IrreversibleReaction class.
//
//-------------------------------------------------------------------------------------------

#include "Reaction.h"

namespace mesmer
{

    class IrreversibleReaction : public Reaction
    {
    public:

        // Constructors.
        IrreversibleReaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id)
            : Reaction(pMoleculeManager, Env, id),
            m_rct1(NULL),
            m_pdt1(NULL),
            m_pdt2(NULL) {} ;

        // Destructor.
        virtual ~IrreversibleReaction(){} ;

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
        virtual double get_relative_pdtZPE() const {
            double zpe = m_pdt1->get_zpe() - getEnv().EMin;
            if (m_pdt2)
                zpe += m_pdt2->get_zpe();
            return zpe;
        }
        virtual double get_relative_TSZPE(void) const {return m_TransitionState->get_zpe() - getEnv().EMin;};

        // Calculate reaction equilibrium constant.
        virtual double calcEquilibriumConstant() ;

        // return the colloptrsize of the reactants
        int getRctColloptrsize(){return m_rct1->get_colloptrsize();}

        // Return products
        int get_products(std::vector<ModelledMolecule *> &product) const
        {
            product.push_back(m_pdt1) ;
            if(m_pdt2){
                product.push_back(m_pdt2) ;
                return 2;
            }
            return 1;
        } ;

        // Return reactants
        ModelledMolecule* get_reactants() const{return m_rct1;} ;

    private:

        // Add reaction terms to collision matrix.
        virtual void AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

        // Read parameters requires to determine reaction heats and rates.
        virtual bool ReadRateCoeffParameters(PersistPtr ppReac) ;

        // Grain averaged microcanonical rate coefficients.
        virtual void calcGrainRateCoeffs();

        // Test k(T)
        virtual void testRateConstant();

        CollidingMolecule   *m_rct1 ;                 // Reactant Molecule.
        ModelledMolecule    *m_pdt1 ;                 // Product Molecule.
        ModelledMolecule    *m_pdt2 ;                 // Subsidiary product molecule.

    } ;


}//namespace
#endif // GUARD_DissociationReaction_h
