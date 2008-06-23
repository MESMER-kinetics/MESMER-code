//-------------------------------------------------------------------------------------------
//
// AssociationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the AssociationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "AssociationReaction.h"
#include <math.h>

using namespace Constants ;
using namespace std;

namespace mesmer
{
    //
    // Read the Molecular for association reaction data from input stream.
    // Note: the convention adopted here is that there are two reactants
    // and one product (adduct).
    //
    bool AssociationReaction::InitializeReaction(PersistPtr ppReac)
    {
        m_ppPersist = ppReac;

        // Read reactant details.

        PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
        Molecule* pMol1 = GetMolRef(ppReactant1);
        if(!pMol1){
            cerr << "Cannot find 1st reactant molecule definition for association reaction " << getName() << ".";
            return false;
        }
        PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        Molecule* pMol2 = GetMolRef(ppReactant2);
        if(!pMol2)
        {
            cerr << "Cannot find 2nd reactant molecule definition for association reaction " << getName() << ".";
            return false;
        }

        // Put the CollidingMolecule into m_rct1, even if it is second in datafile

        CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
        if(pColMol){
            m_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);
        } else {
            pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
            if(!pColMol){
                cerr << "At least one of the reactants has to be a colliding molecule";
                return false;
            }
            m_rct2 = dynamic_cast<ModelledMolecule*>(pMol1);
        }
        m_rct1 = pColMol;

        // SuperMolecule section --- to be removed. ------------------------------
        string id; //shoud not set any name for it.
        SuperMolecule* pSupMol = NULL;
        while(m_pMoleculeManager->GetNextMolecule(id, pSupMol)){ // get next SuperMolecule
            // if found a SuperMolecule
            ModelledMolecule* rm1 = pSupMol->getMember1();
            ModelledMolecule* rm2 = pSupMol->getMember2();
            if (!rm1 && !rm2){// there is no data inside, occupy it!
                pSupMol->setMembers(m_rct1, m_rct2);
                m_srct = pSupMol;
                break;
            }
        }
        if (!pSupMol){
            cerr << "No SuperMolecule was found.";
            // there will always at least one SuperMolecule in m_molmap, check the end of addmol()
            // in MoleculeManager.cpp.
        }
        // SuperMolecule section --- to be removed. ------------------------------

        //Read product details.

        PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
        pMol1 = GetMolRef(ppProduct1);
        if (!pMol1) {
            cerr << "Cannot find product molecule definition for association reaction " << getName() << ".";
            return false;
        }

        // Save product as CollidingMolecule.

        pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
        if(pColMol){
            m_pdt1 = pColMol;
        } else {
            cerr << "Isomer product must be a colliding molecule";
            return false;
        }

        // Read the transition state (if present)
        PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
        if (ppTransitionState)
        {
            TransitionState* pTrans = dynamic_cast<TransitionState*>(GetMolRef(ppTransitionState));
            if(pTrans) m_TransitionState = pTrans;
        }

        // Read heat of reaction and rate parameters.

        return ReadRateCoeffParameters(ppReac) ;
    }

    //
    // Add (REVERSIBLE) association reaction terms to collision matrix.
    //
    void AssociationReaction::AddReactionTerms(dMatrix      *CollOptr,
        isomerMap    &isomermap,
        const double rMeanOmega)
    {
        // Get densities of states of the adduct for detailed balance.
        vector<double> pdtDOS;
        m_pdt1->getGrainDensityOfStates(pdtDOS) ;

        // Locate isomers in system matrix.
        const int pdtLoc =      isomermap[m_pdt1] ;
        const int jj     = (*m_sourceMap)[m_srct] ;

        // Get equilibrium constant.
        const double Keq = calcEquilibriumConstant() ;

        // Get Boltzmann distribution for detailed balance.
        const int MaximumGrain = getEnv().MaxGrn ;
        vector<long double> adductPopFrac ; // Population fraction of the adduct

        const int pdtGrainZPE = m_pdt1->get_grnZpe();
        m_pdt1->normalizedBoltzmannDistribution(adductPopFrac, MaximumGrain - pdtGrainZPE) ;

        double DissRateCoeff(0.0) ;

        const int TSFluxGrainZPE  = getTSFluxGrnZPE();

        for ( int i = TSFluxGrainZPE, j = 0; i < MaximumGrain; ++i, ++j) {
            int ll(i - pdtGrainZPE);
            int ii(pdtLoc + ll) ;

            (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainTSFlux[j] / pdtDOS[ll];                                // Loss of the adduct to the source
            (*CollOptr)[jj][ii]  = rMeanOmega * m_GrainTSFlux[j] * sqrt(adductPopFrac[ll] * Keq) / pdtDOS[ll];// Reactive gain of the source
            (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii] ;                                                      // Reactive gain (symmetrization)
            DissRateCoeff       += m_GrainTSFlux[j] * adductPopFrac[ll] / pdtDOS[ll];
        }
        (*CollOptr)[jj][jj] -= (rMeanOmega * DissRateCoeff * Keq);       // Loss of the source from detailed balance.
    }

    //
    // Calculate reaction equilibrium constant for the general reaction
    //        A + B  <===> C
    //
    double AssociationReaction::calcEquilibriumConstant() {

        // equilibrium constant:
        long double Keq(0.0) ;
        const long double beta = getEnv().beta ;

        // partition function for each reactant
        long double Qrcts = m_srct->rovibronicGrnCanPrtnFn();

        // rovibronic partition function for reactants multiplied by translation contribution
        Qrcts *= translationalContribution(m_rct1->getMass(), m_rct2->getMass(), beta);

        // rovibronic partition function for product
        const long double Qpdt1 = m_pdt1->rovibronicGrnCanPrtnFn() ;

        Keq = Qpdt1 / Qrcts;

        /* Electronic degeneracies were already accounted for in DOS calculations */

        // Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells
        const double HeatOfReaction = getHeatOfReaction() ;
        const long double _expon = -beta * HeatOfReaction;
        Keq *= exp(_expon) ;

        const double excess = m_srct->getExcessReactantConc();
        Keq *= excess ;
        //
        // K_eq = ( [C]/[A][B] ) * [A] = [C]/[B]
        //
        // where [A] is the reactant what is in excess (seen as constant).
        // Therefore, the K_eq here is essentially the pseudo-first-order equilibrium constant.

        return (double) Keq ;

    }

    // Read parameters requires to determine reaction heats and rates.
    bool AssociationReaction::ReadRateCoeffParameters(PersistPtr ppReac) {

        // Read the heat of reaction (if present).
        const char* pHeatRxntxt = ppReac->XmlReadValue("me:HeatOfReaction",false);
        if (pHeatRxntxt){
            stringstream s1(pHeatRxntxt);
            double value = 0.0; s1 >> value ; setHeatOfReaction(value);
        } else { // Calculate heat of reaction.
            setHeatOfReaction(get_relative_pdtZPE(), get_relative_rctZPE());
        }

        const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy",false);
        if (pActEnetxt)
        {
            double value = 0.0;
            stringstream s2(pActEnetxt); s2 >> value; m_ActivationEnergy = value;
        }
        const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential",false);
        if (pPreExptxt)
        {
            double value = 0.0;
            stringstream s2(pPreExptxt); s2 >> value ; set_PreExp(value);
        }
        const char* pNInftxt   = ppReac->XmlReadValue("me:nInfinity",false);
        if (pNInftxt)
        {
            double value = 0.0; stringstream s3(pNInftxt); s3 >> value ; set_NInf(value);
        }
        const char* pERConctxt   = ppReac->XmlReadValue("me:excessReactantConc",false);
        if (pERConctxt)
        {
            stringstream s3(pERConctxt);
            double erconc = 0.0;
            s3 >> erconc;
            m_srct->setExcessReactantConc(erconc);
        }

        // Determine the method of MC rate coefficient calculation.
        const char* pMCRCMethodtxt = ppReac->XmlReadValue("me:MCRCMethod") ;
        if(pMCRCMethodtxt)
        {
            m_pMicroRateCalculator = MicroRateCalculator::Find(pMCRCMethodtxt);
            if(!m_pMicroRateCalculator)
            {
                cerr << "Unknown method " << pMCRCMethodtxt
                    << " for the determination of Microcanonical rate coefficients in reaction "
                    << getName();
                return false;
            }
        }

        // Determine the method of estimating tunneling effect.
        const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling") ;
        if(pTunnelingtxt)
        {
            m_pTunnelingCalculator = TunnelingCalculator::Find(pTunnelingtxt);
            if(!m_pTunnelingCalculator)
            {
                cerr << "Unknown method " << pTunnelingtxt
                    << " for the determination of tunneling coefficients in reaction "
                    << getName();
                return false;
            }
        }
        else{
            cinfo << "No tunneling method was found for " << getName() << endl;
        }

        return true ;
    }

    //
    // Calculate grained forward and reverse k(E)s from trainsition state flux
    //
    void AssociationReaction::calcGrainRateCoeffs(){

        vector<double> rctGrainDOS;
        vector<double> pdtGrainDOS;
        m_srct->getGrainDensityOfStates(rctGrainDOS) ;
        m_pdt1->getGrainDensityOfStates(pdtGrainDOS) ;

        const int TSFluxGrainZPE  = getTSFluxGrnZPE();
        const int rctGrainZPE = m_srct->get_grnZpe();
        const int pdtGrainZPE = m_pdt1->get_grnZpe();

        const int MaximumGrain = getEnv().MaxGrn;
        m_GrainKfmc.clear();
        m_GrainKfmc.resize(MaximumGrain , 0.0);
        m_GrainKbmc.clear();
        m_GrainKbmc.resize(MaximumGrain , 0.0);

        // For AssociationReaction, TSFlux is calculated from ILT
        for (int i = TSFluxGrainZPE - pdtGrainZPE, j = 0; i < MaximumGrain; ++i, ++j){
            m_GrainKbmc[i] = m_GrainTSFlux[j] / pdtGrainDOS[i];
        }
        for (int i = TSFluxGrainZPE - rctGrainZPE, j = 0; i < MaximumGrain; ++i, ++j){
            m_GrainKfmc[i] = m_GrainTSFlux[j] / rctGrainDOS[i];
        }

        // the code that follows is for printing of the f & r k(E)s
        if (getEnv().kfEGrainsEnabled){
            ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
            for (int i = 0; i < MaximumGrain; ++i){
                ctest << m_GrainKfmc[i] << endl;
            }
            ctest << "}\n";
        }
        if (getEnv().kbEGrainsEnabled){
            ctest << "\nk_b(e) grains for " << getName() << ":\n{\n";
            for (int i = 0; i < MaximumGrain; ++i){
                ctest << m_GrainKbmc[i] << endl;
            }
            ctest << "}\n";
        }
    }


    //this function retrieves the activation energy for an association reaction
    double AssociationReaction::get_ThresholdEnergy(void) const {
        if (!m_ActivationEnergy){
            if (!m_TransitionState) {
                cinfo << "No TransitionState for " << getName() << ", threshold energy = 0.0" << endl;
                return 0.0;
            }
            double ThresholdEnergy = get_relative_TSZPE() - get_relative_rctZPE();
            // Activation energy should be defined by user, otherwise return zero.
            return ThresholdEnergy;
        }
        return m_ActivationEnergy;
    } ;

}//namespace
