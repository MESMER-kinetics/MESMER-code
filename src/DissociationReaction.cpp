//-------------------------------------------------------------------------------------------
//
// DissociationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the DissociationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "DissociationReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
    //
    // Read the Molecular data from input stream.
    //
    bool DissociationReaction::InitializeReaction(PersistPtr ppReac)
    {
        m_ppPersist = ppReac;

        // Read reactant details.

        PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
        Molecule* pMol1 = GetMolRef(ppReactant1);
        if(!pMol1){
            string errorMsg = "Dissociation reaction " + getName() + " has no reactant.";
            meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
            return false;
        }        
        PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        if(ppReactant2)
        {
            string errorMsg = "Dissociation reaction " + getName() + " has more than one reactant.";
            meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
            return false;
        }

        // Save reactant as CollidingMolecule.

        CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
        if(pColMol){
            m_rct1 = pColMol;
        } else {
            meErrorLog.ThrowError(__FUNCTION__, string("Dissociating reactant must be a colliding molecule"), obError);
            return false;
        }

        // Read product details. The detail of products may be absent or, may be needed
        // to calculate the microcanonical rates. If there are products, save them as 
        // type Molecule.

        PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
        if (ppProduct1) {
            pMol1 = GetMolRef(ppProduct1);
            if (pMol1) {
                m_pdt1 = dynamic_cast<CollidingMolecule*>(pMol1) ;
            } else {
                string errorMsg = "Dissociation reaction" + getName() + " has no products defined.";
                meErrorLog.ThrowError(__FUNCTION__, errorMsg, obWarning);
            }

            Molecule* pMol2 = NULL ;
            PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
            if (ppProduct2) {
                pMol2 = GetMolRef(ppProduct1);
                if (pMol2) {
                    m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2) ;
                } else {
                    string errorMsg = "Dissociation reaction " + getName() + " has only one product defined.";
                    meErrorLog.ThrowError(__FUNCTION__, errorMsg, obWarning);
                }
            }
        }

        // Read the transition state (if present).

        PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
        if (ppTransitionState)
        {
            TransitionState* pTrans = dynamic_cast<TransitionState*>(GetMolRef(ppTransitionState));
            if(pTrans) 
                m_TransitionState = pTrans;
        }

        // Read heat of reaction and rate parameters.

        return ReadRateCoeffParameters(ppReac) ;
    }

    //
    // Calculate reaction equilibrium constant.
    //
    double DissociationReaction::calcEquilibriumConstant() {

        double Keq(0.0) ;

        return Keq ;
    }

    //
    // Add dissociation reaction terms to collision matrix.
    //
    void DissociationReaction::AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

        // Locate reactant in system matrix.
        const int rctLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_rct1)] ;
        const int colloptrsize = dynamic_cast<CollidingMolecule*>(m_rct1)->get_colloptrsize();

        if (collisionOperatorCheck){
            ctest << "\nSystem collision operator check before adding " << getName() << " microrates:\n";
            (*CollOptr).showFinalBits(8);
        }

        for ( int i = 0 ; i < colloptrsize; ++i ) {
            int ii(rctLocation + i) ;
            (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[i] ;   // Forward loss reaction.
        }

        if (collisionOperatorCheck){
            ctest << "\nSystem collision operator check after adding " << getName() << " microrates:\n";
            (*CollOptr).showFinalBits(8);
        }

    }


}//namespace
