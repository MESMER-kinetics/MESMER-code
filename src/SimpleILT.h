#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

//-------------------------------------------------------------------------------------------
//
// SimpleILT.h 
//
// Author: Struan Robertson 
// Date:   1/Jul/2007
//
// This header file contains the declaration of the SimpleILT class.
//
//-------------------------------------------------------------------------------------------

#include <vector>
#include "Molecule.h"
#include "MoleculeManager.h"
#include "Persistence.h"
#include "MicroRateCalculator.h"

namespace mesmer
{
	class SimpleILT : public MicroRateCalculator
	{

	public:

		SimpleILT(MoleculeManager *pMoleculeManager, CollidingMolecule *pReactant): 
		  m_pMoleculeManager(pMoleculeManager),
			  m_Reactant(pReactant),
			  m_ActEne(0.0),
			  m_PreExp(0.0) {}

		  virtual ~SimpleILT() {}

		  virtual void calculateMicroRateCoeffs(std::vector<double> &kfmc)  ; 
		  bool Initialize(PersistPtr &ppReac) ;

	private:

		MoleculeManager   *m_pMoleculeManager ; // Pointer to molecule manager.
		CollidingMolecule *m_Reactant ;         // Reactant Molecule.
		double             m_ActEne ;           // Arrhenius activation energy.
		double             m_PreExp ;           // Arrhenius Pre-exponential factor.
	} ;


}//namespace
#endif // GUARD_SimpleILT_h
