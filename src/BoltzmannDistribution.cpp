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
#include "BoltzmannDistribution.h"

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  BoltzmannDistribution theBoltzmannDistribution("Boltzmann");
  //************************************************************

  bool BoltzmannDistribution::calculateDistribution(const std::vector<double>& DOS,
                                                    const std::vector<double>& Ene,
                                                    const double& beta,
                                                    std::vector<double>& dist,
                                                    double& prtfn)
  {
    int v_size = static_cast<int>(DOS.size());
    dist.resize(v_size, 0.0);
    // Calculate the Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.

    prtfn = 0.;
    for (int i = 0; i < v_size; ++i) {
      double tmp = log(DOS[i]) - beta * Ene[i] + 10.0 ;
      tmp = exp(tmp) ;
      prtfn += tmp ;
      dist[i] = tmp ;
    }
    return true;
  }
}//namespace

