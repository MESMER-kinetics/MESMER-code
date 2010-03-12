//-------------------------------------------------------------------------------------------
//
// EckartCoefficients.h
//
// Author: Dave Glowacki, based on fortran code written by Jeremy Harvey
// Date:   28-8-2009
//
// Produces WKB spin forbidden crossing coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../System.h"

namespace mesmer
{
	class WKBCrossingCoeff : public CrossingCalculator
	{
	public:

		///Constructor which registers with the list of CrossingCalculators in the base class
		WKBCrossingCoeff(const std::string& id) : CrossingCalculator(id){}

		virtual ~WKBCrossingCoeff() {}

		virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability);

		virtual bool ThereIsTunnellingWithCrossing(void) {return true;};

		bool ReadDoubleAndUnits(double& element, PersistPtr pp, const std::string identifier, const std::string units);
	};
}//namespace


