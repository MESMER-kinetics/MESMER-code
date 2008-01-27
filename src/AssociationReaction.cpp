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
using namespace mesmer;

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
			string errorMsg = "Association reaction " + getName() + " has no reactants.";
			meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
			return false;
		}        
		PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
		Molecule* pMol2 = GetMolRef(ppReactant2);
		if(!pMol2)
		{
			string errorMsg = "Association reaction " + getName() + " does not have two reactants.";
			meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
			return false;
		}

		// Put the CollidingMolecule into m_rct1, even if it is second in datafile

		CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
		if(pColMol){
			m_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);
		} else {
			pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
			if(!pColMol){
				meErrorLog.ThrowError(__FUNCTION__, string("At least one of the reactants has to be a colliding molecule"), obError);
				return false;
			}
			m_rct2 = dynamic_cast<ModelledMolecule*>(pMol1);
		}
		m_rct1 = pColMol;

		if (m_rct1 && m_rct2){ // the reactant side has two molecules
			{stringstream errorMsg;
			errorMsg << "Reaction " << getName() << " has two reactants. ";
			meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

			// check whether there is any SuperMolecule in m_molmap contains pMol1 & pMol2
			string id; //shoud not set any name for it.
			SuperMolecule* pSupMol = NULL;
			while(m_pMoleculeManager->GetNextMolecule(id, pSupMol)){ // get next SuperMolecule
				// if found a SuperMolecule
				ModelledMolecule*  rm1 = pSupMol->getMember1();
				ModelledMolecule*  rm2 = pSupMol->getMember2();
				if (!rm1 && !rm2){// there is no data inside, occupy it!
					pSupMol->setMembers(m_rct1, m_rct2);
					m_srct = pSupMol;
					{stringstream errorMsg;
					errorMsg << "Set members of the SuperMolecule: " << m_srct->getName();
					meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}
					break;
				}
			}
			if (!pSupMol){
				meErrorLog.ThrowError(__FUNCTION__, string("No SuperMolecule was found."), obInfo);
				// there will always at least one SuperMolecule in m_molmap, check the end of addmol()
				// in MoleculeManager.cpp.
				/* need to create one (mark _2007_12_10__17_10_18_)
				write a SuperMolecule creator that acquire a position in the XML
				*/
			}
		}
		else{
			string errorMsg = "Reaction " + getName() + " has only one reactant";
			meErrorLog.ThrowError(__FUNCTION__, errorMsg, obInfo);
		}

		//Read product details.

		PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
		pMol1 = GetMolRef(ppProduct1);
		if (!pMol1) {
			string errorMsg = "Association reaction " + getName() + " has no product.";
			meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
			return false;
		}           
		PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
		if(ppProduct2)
		{
			string errorMsg = "Association reaction " + getName() + " has more than one product.";
			meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
			return false;
		}

		// Save product as CollidingMolecule.

		pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
		if(pColMol){
			m_pdt1 = pColMol;
		} else {
			meErrorLog.ThrowError(__FUNCTION__, string("Isomer product must be a colliding molecule"), obError);
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
		// Locate isomers in system matrix.
		const int pdtLoc =      isomermap[dynamic_cast<CollidingMolecule*>(m_pdt1)] ;
		const int sL     = (*m_sourceMap)[dynamic_cast<SuperMolecule    *>(m_srct)] ;

		// Get equilibrium constant.
		double Keq = calcEquilibriumConstant() ;

		// Multiply equilibrum constant by concentration of excess reactant.
		// concentration of excess reactant should be in molec/cm3. This gives
		// a dimensionless pseudo-isomerization equilibrium constant.

        Keq *= m_excessReactantConc ;

		// Get Boltzmann distribution for detailed balance.

		const int MaximumGrain = getEnv().MaxGrn ;
		vector<double> adductBoltz(MaximumGrain, 0.0) ;
		m_pdt1->grnBoltzDist(adductBoltz) ;

		const int colloptrsize = dynamic_cast<CollidingMolecule*>(m_pdt1)->get_colloptrsize() ;

		double DissRateCoeff(0.0) ;

		const int idx = m_pdt1->get_grnZpe() - m_rct1->get_grnZpe() ;
		for ( int i = max(0,-idx) ; i < min(colloptrsize,(colloptrsize-idx)) ; ++i ) {
			int ll = i + idx ;
			int pL(pdtLoc + ll) ;

			(*CollOptr)[pL][pL] -= rMeanOmega * m_GrainKfmc[ll] ;                           // Forward loss reaction.
			(*CollOptr)[pL][sL]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(adductBoltz[ll]*Keq) ; // Reactive gain.
			(*CollOptr)[sL][pL]  = (*CollOptr)[pL][sL] ;                                    // Reactive gain.
			DissRateCoeff       += m_GrainKfmc[ll]*adductBoltz[ll] ;
		}
		(*CollOptr)[sL][sL] -= DissRateCoeff*Keq ;       // Backward loss reaction from detailed balance.
	}

	//
	// Calculate reaction equilibrium constant for the general reaction
	//        A + B  <===> C
	//
	double AssociationReaction::calcEquilibriumConstant() {

		double Keq(0.0) ;

		// Get Canonical rovibronic partition functions.

		const double Qrct1 = m_rct1->grnCanPrtnFn() ;
		const double Qrct2 = m_rct2->grnCanPrtnFn() ;
		const double Qpdt1 = m_pdt1->grnCanPrtnFn() ;

		{stringstream errorMsg;
		errorMsg << "Qrct1 = " << Qrct1 << ", Qpdt1 = " << Qpdt1 << ", Qrct2 = " << Qrct2 ;
		meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

		// Get mass.

		const double mrct1 = m_rct1->getMass() ;
		const double mrct2 = m_rct2->getMass() ;
		const double mpdt1 = mrct1 + mrct2 ;

		// Calculate the equilibrium constant.

		const double beta = getEnv().beta ;

		const double Tau = 2.0593e+19*pow( (2.0*acos(-1.0)), 1.5) ; // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> J)/(h2*Na)))^3/2

		Keq  = Qpdt1/(Qrct1*Qrct2) ;                             // Rovibronic contribution.
        Keq *= pow( ((mpdt1*beta)/(mrct1*mrct2)), 1.5)/Tau ;     // Translational contribution.
		Keq *= exp(-beta*m_HeatOfReaction) ;

		return Keq ;
	}

}//namespace
