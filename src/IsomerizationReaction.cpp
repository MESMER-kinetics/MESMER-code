//-------------------------------------------------------------------------------------------
//
// IsomerizationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IsomerizationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IsomerizationReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
    // Read the Molecular data from input stream.
    //
    bool IsomerizationReaction::InitializeReaction(PersistPtr ppReac)
    {
        m_ppPersist = ppReac;

        // Read reactant details.

		PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
        Molecule* pMol1 = GetMolRef(ppReactant1);
        if(!pMol1){
            string errorMsg = "Isomerizaton reaction " + m_Name + " has no reactant.";
            meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
            return false;
		}        
		PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        if(ppReactant2)
        {
            string errorMsg = "Isomerizaton reaction " + m_Name + " has more than one reactant.";
            meErrorLog.ThrowError(__FUNCTION__, errorMsg, obError);
            return false;
        }

        // Save reactant as CollidingMolecule.

        CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
        if(pColMol){
			m_rct1 = pColMol;
        } else {
            meErrorLog.ThrowError(__FUNCTION__, string("Isomer reactant must be a colliding molecule"), obError);
            return false;
        }

        //Read product details.

        PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
		pMol1 = GetMolRef(ppProduct1);
        if (!pMol1) {
            string errorMsg = "Isomerizaton reaction " + m_Name + " has no product.";
            meErrorLog.ThrowError(__FUNCTION__, string(errorMsg), obError);
            return false;
		}           
        PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
        if(ppReactant2)
        {
            string errorMsg = "Isomerizaton reaction " + m_Name + " has more than one product.";
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

        // Read the transition state (if present).

        PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
        if (ppTransitionState)
        {
            TransitionState* pTrans = dynamic_cast<TransitionState*>(GetMolRef(ppTransitionState));
            if(pTrans) 
				m_TransitionState = pTrans;

            /* It would be better to use the ZPEs rather than threshold
            if(m_TransitionState)
            {
            const char* pthreshtxt = ppReac->XmlReadValue("me:threshold",false);
            if(pthreshtxt)
            {
            stringstream ss(pthreshtxt);
            ss >> m_E0;
            }
            */
        }

        // Read heat of reaction and rate parameters.

        return ReadRateCoeffParameters(ppReac) ;
    }

	//
  // Add isomer reaction terms to collision matrix.
  //
  void IsomerizationReaction::AddReactionTerms(dMatrix         *CollOptr,
                                               isomerMap       &isomermap,
                                               const double    rMeanOmega)
  {
    // Locate isomers in system matrix.
    const int rctLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_rct1)] ;
    const int pdtLocation = isomermap[dynamic_cast<CollidingMolecule*>(m_pdt1)] ;

    // Get densities of states for detailed balance.
    const int MaximumGrain = getEnv().MaxGrn;
    vector<double> rctDOS(MaximumGrain, 0.0) ;
    vector<double> pdtDOS(MaximumGrain, 0.0) ;

    m_rct1->grnDensityOfStates(rctDOS) ;
    m_pdt1->grnDensityOfStates(pdtDOS) ;

    const int idx = m_pdt1->get_grnZpe() - m_rct1->get_grnZpe() ;
    for ( int i = max(0,-idx) ; i < min(MaximumGrain,(MaximumGrain-idx)) ; ++i ) {
      int ll = i + idx ;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation + i) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[ll] ;                            // Forward loss reaction.
      (*CollOptr)[jj][jj] -= rMeanOmega * m_GrainKfmc[ll]*rctDOS[ll]/pdtDOS[i] ;       // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(rctDOS[ll]/pdtDOS[i]) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                     // Reactive gain.
    }
  }


}//namespace
