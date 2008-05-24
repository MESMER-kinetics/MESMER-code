//-------------------------------------------------------------------------------------------
//
// BoltzmannDistribution.h
//
// Author: Chi-Hsiu Liang
// Date:   _2008_05_15_
//
// Produces Boltzmann distribution
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "Distribution.h"

namespace mesmer
{
  class BoltzmannDistribution : public DistributionCalculator
  {
  public:
  
    ///Constructor which registers with the list of TunnelingCalculators in the base class
    BoltzmannDistribution(const std::string& id) : DistributionCalculator(id){}
  
    virtual ~BoltzmannDistribution() {}
  
    virtual bool calculateDistribution(const std::vector<double>& DOS, const std::vector<double>& ene, const double& beta, std::vector<double>& distribution);
  };
}//namespace


