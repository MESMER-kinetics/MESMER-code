//-------------------------------------------------------------------------------------------
//
// ExchangeReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the ExchangeReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "ExchangeReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
    //
    // Read the Molecular data from input stream.
    //
    bool ExchangeReaction::InitializeReaction(PersistPtr ppReac)
    {
        m_ppPersist = ppReac;

        // Read reactant details.

        PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
        Molecule* pMol1 = GetMolRef(ppReactant1);
        if(!pMol1){
            string errorMsg = "Exchange reaction " + getName() + " has no reactant.";
            meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
            return false;
        }        
        PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        Molecule* pMol2 = GetMolRef(ppReactant2);
        if(!pMol2)
        {
            string errorMsg = "Exchange reaction " + getName() + " has only one reactant.";
            meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
            return false;
        }

        // Save first reactant as CollidingMolecule. SHR 6/Apr/2008 this incorrect,
        // but we are forced at present to do this because of the types declared in
        // the reaction class. 

        CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
        if(pColMol){
            m_rct1 = pColMol;
        } else {
            meErrorLog.ThrowError(__FUNCTION__, string("Exchange reactant missing"), obError);
            return false;
        }
        m_rct2 = dynamic_cast<ModelledMolecule*>(pMol2) ;
        if (m_rct2) {
            meErrorLog.ThrowError(__FUNCTION__, string("Exchange reactant missing"), obError);
            return false;
        }

        // Read product details. Save them as type Molecule.

        PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
        if (ppProduct1) {
            pMol1 = GetMolRef(ppProduct1);
            if (pMol1) {
                m_pdt1 = dynamic_cast<CollidingMolecule*>(pMol1) ;
            } else {
                string errorMsg = "Exchange reaction" + getName() + " has no products defined.";
                meErrorLog.ThrowError(__FUNCTION__, errorMsg, obWarning);
            }

            PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
            if (ppProduct2) {
                pMol2 = GetMolRef(ppProduct1);
                if (pMol2) {
                    m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2) ;
                } else {
                    string errorMsg = "Exchange reaction " + getName() + " has only one product defined.";
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
    double ExchangeReaction::calcEquilibriumConstant() {
        // equilibrium constant:

        double Keq(0.0) ;

        // Get Canonical partition functions.

        double Qrcts = 1.0;
        if (m_rct2)
            Qrcts = m_srct->grnCanPrtnFn();
        else
            Qrcts = m_rct1->grnCanPrtnFn() ;

        double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

        double mass_rct1 = m_rct1->getMass();
        double mass_rct2 = (m_rct2)? m_rct2->getMass() : 0.0;
        double mass_srct = mass_rct1 + mass_rct2;

        // Calculate the equilibrium constant.
        double beta = getEnv().beta ;

        Keq = Qpdt1/Qrcts ;
        if(debugFlag) ctest << "Keq = " << Keq << endl;

        /* Electronic degeneracies were already accounted for in DOS calculations */
        // Heat of reaction
        Keq *= exp(-beta * getHeatOfReaction() * kJPerMolInRC) ;
        if(debugFlag) ctest << "Keq = " << Keq << endl;

        // Translational contribution
        // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
        // double tau = 2.0593e19 * pow(2. * M_PI ,1.5);
        if (m_rct2)
            Keq /= (tp_C * pow(mass_rct1 * mass_rct2 / (mass_srct * beta), 1.5));

        if(debugFlag) ctest << "Keq = " << Keq << endl;
        return Keq ;
    }

    //
    // Add exchange reaction terms to collision matrix.
    //
    void ExchangeReaction::AddReactionTerms(dMatrix      *CollOptr,
        isomerMap    &isomermap,
        const double rMeanOmega)
    {

    }

}//namespace
