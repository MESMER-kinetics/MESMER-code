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
    BoltzmannDistribution(const char* id) :m_id(id){ Register(); }
  
    virtual ~BoltzmannDistribution() {}
    virtual const char* getID() override { return m_id; }

    virtual bool calculateDistribution(std::vector<double> DOS, std::vector<double> ene, const double& beta, std::vector<double>& distribution);
  private:
    const char* m_id;
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
    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    dist.clear();
    for (size_t i(0) ; i < DOS.size() ; ++i) {
      dist.push_back(exp(log(DOS[i]) - beta * Ene[i] + 10.0)) ;
    }

    // Calculate the normalization coefficient. Reverse for numerical accuracy.
    double sum(0.0) ;
    for (size_t i(0), j(dist.size()-1) ; i < dist.size() ; ++i,--j) {
      sum += dist[j] ;
    }
    
    // Normalize distribution. 
    for (size_t i(0) ; i < dist.size() ; ++i) {
      dist[i] /= sum ;
    }
    
    return true;
  }
}//namespace

