#ifndef GUARD_MicroRateCalculator_h
#define GUARD_MicroRateCalculator_h

//-------------------------------------------------------------------------------------------
//
// MicroRateCalculator.h 
//
// Author: Struan Robertson 
// Date:   21/Jul/2007
//
// This header file contains the declaration of the MicroRateCalculator class - Base class
// for Microcanonical rate calculators
//
//-------------------------------------------------------------------------------------------

namespace mesmer
{
	//  

	class MicroRateCalculator {

	public:

		MicroRateCalculator(){} ;

		virtual ~MicroRateCalculator(){} ;

		virtual void calculateMicroRateCoeffs(std::vector<double> &kfmc) = 0 ;

	protected:

		void testMicroRateCoeffs(std::vector<double> &kfmc, CollidingMolecule *m_Reactant) const;

	};


}//namespace
#endif // GUARD_MicroRateCalculator_h
