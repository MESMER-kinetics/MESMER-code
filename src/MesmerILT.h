#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

#include "System.h"
#include "ClassicalRotor.h"

namespace mesmer
{
	class MesmerILT : public MicroRateCalculator
	{
	public:

		///Constructor which registers with the list of MicroRateCalculators in the base class
		MesmerILT(const std::string& id) : MicroRateCalculator(id) { 
			m_pDensityOfStatesCalculator = dynamic_cast<ClassicalRotor*>(DensityOfStatesCalculator::Find("Classical rotors")) ;
		}

		virtual ~MesmerILT() {}

		virtual bool calculateMicroRateCoeffs(Reaction* pReact) ;

	private:

		ClassicalRotor *m_pDensityOfStatesCalculator ;

	};
}//namespace

#endif // GUARD_SimpleILT_h
