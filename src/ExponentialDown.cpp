
#include "ExponentialDown.h"
#include <cmath>

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance).
  ExponentialDown exponentialDown("ExponentialDown");
  //************************************************************

  double ExponentialDown::calculateTransitionProbability(double Ei, double Ej) {
    // return exp(-(Ei -Ej)/m_deltaEdown) ;
    return exp(-m_alpha*(Ei -Ej)) ;
  }

  void ExponentialDown::ReadParameters(double deltaEdown) { 
    m_deltaEdown =  deltaEdown;
    m_alpha      =  1.0/deltaEdown;
  }

}//namespace

