
#include "ExponentialDown.h"
#include <cmath>

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance).
  ExponentialDown exponentialDown("ExponentialDown");
  //************************************************************

  double ExponentialDown::calculateTransitionProbability(double Ei, double Ej) {
    return exp(-(Ei -Ej)/m_deltaEdown) ;
  }

  bool ExponentialDown::ReadParameters() { 
    return true ;
  }

}//namespace

