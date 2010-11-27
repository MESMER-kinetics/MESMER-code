//-------------------------------------------------------------------------------------------
//
// BoltzmannDistribution.cpp
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
#include "../Distribution.h"

namespace mesmer
{
  class BoltzmannDistribution : public DistributionCalculator
  {
  public:
  
    ///Constructor which registers with the list of DistributionCalculators in the base class
    BoltzmannDistribution(const std::string& id) : DistributionCalculator(id){}
  
    virtual ~BoltzmannDistribution() {}
  
    virtual bool calculateDistribution(std::vector<double> DOS, std::vector<double> ene, const double& beta, std::vector<double>& distribution);
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  BoltzmannDistribution theBoltzmannDistribution("Boltzmann");
  //************************************************************

  bool BoltzmannDistribution::calculateDistribution(std::vector<double> DOS,
                                                    std::vector<double> Ene,
                                                    const double& beta,
                                                    std::vector<double>& dist)
  {
    // The Boltzmann distribution applies only to particles at a high enough temperature
    // and low enough density that quantum effects can be ignored, and the particles are
    // obeying Maxwell–Boltzmann statistics.

    dist.clear();
    int vsize = static_cast<int>(DOS.size());
    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    for (int i = 0; i < vsize; ++i) {
      dist.push_back(exp(log(DOS[i]) - beta * Ene[i] + 10.0));
    }
    return true;
  }
}//namespace

