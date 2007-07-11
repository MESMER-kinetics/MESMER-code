#ifndef GUARD_SimpleRRKM_h
#define GUARD_SimpleRRKM_h

//-------------------------------------------------------------------------------------------
//
// SimpleRRKM.h 
//
// Author: Struan Robertson 
// Date:   1/Jul/2007
//
// This header file contains the declaration of the SimpleRRKM class.
//
//-------------------------------------------------------------------------------------------

#include <vector>
#include "Molecule.h"
#include "MoleculeManager.h"
#include "Persistence.h"
#include "Reaction.h"

namespace mesmer
{
	class SimpleRRKM : public MicroRateCalculator
	{

	public:

		SimpleRRKM(MoleculeManager *pMoleculeManager, CollidingMolecule *pReactant): 
		  m_pMoleculeManager(pMoleculeManager),
			  m_Reactant(pReactant),
			  m_TransitionState(NULL),
			  m_E0(0.0) {}

		  virtual ~SimpleRRKM() {}

		  virtual void calculateMicroRateCoeffs(std::vector<double> &kfmc)  ; 
		  bool Initialize(PersistPtr &ppReac) ;

	private:

		void testMicroRateCoeffs(std::vector<double> &kfmc) const;

		MoleculeManager   *m_pMoleculeManager ; // Pointer to molecule manager.
		CollidingMolecule *m_Reactant ;         // Reactant Molecule.
		TransitionState   *m_TransitionState ;  // Transition State.
		double             m_E0 ;               // Reaction Threshold energy (measured from zero point of reactant).
	} ;


}//namespace
#endif // GUARD_SimpleRRKM_h
