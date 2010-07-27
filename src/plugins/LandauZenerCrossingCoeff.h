//-------------------------------------------------------------------------------------------
//
// EckartCoefficients.h
//
// Author: Dave Glowacki, based on fortran code written by Jeremy Harvey
// Date:   28-8-2009
//
// Produces Landau-Zener spin forbidden crossing coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../System.h"

namespace mesmer
{
  class LandauZenerCrossingCoeff : public CrossingCalculator
  {
  public:
  
    ///Constructor which registers with the list of CrossingCalculators in the base class
    LandauZenerCrossingCoeff(const std::string& id) : CrossingCalculator(id){}
  
    virtual ~LandauZenerCrossingCoeff() {}
  
    virtual bool calculateCellCrossingCoeffs(Reaction* pReact, std::vector<double>& CrossingProbability);

		virtual bool ThereIsTunnellingWithCrossing(void) {return false;};

		bool ReadDoubleAndUnits(double& element, PersistPtr pp, const std::string identifier, const std::string units);
  };
}//namespace


