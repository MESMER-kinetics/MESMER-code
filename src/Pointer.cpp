//----------------------------------------------------------------------------
// pointer class
// Author: Chi-Hsiu Liang
// A class to cast double values into a dynamic vector for grouped optimization
//----------------------------------------------------------------------------

#include "Pointer.h"
#include "oberror.h"

namespace mesmer
{
  vector<double> allP;
  vector<DPoint> fitDP;

  void DPoint::set_range(const double valueL, const double valueU, const double stepsize_){
    allP.push_back(valueL);
    p_ = int(allP.size()) - 1;

    if (valueU == valueL){
      cerr << "Upper value cannot equal to lower value.";
      return;
    }
    else if (valueU > valueL){
      lower = valueL;
      upper = valueU;
    }
    else{
      lower = valueU;
      upper = valueL;
    }

    stepsize = abs(stepsize_);

    int numStep = int((upper - lower)/stepsize);
    if (numStep < 1){
      cerr << "There is only one point for this variable. Remove range setting to clear this error.";
    }

    fitDP.push_back(*this);
  }

  void DPoint::get_range(double& lower_, double& upper_, double& stepsize_){
    lower_ = lower;
    upper_ = upper;
    if (stepsize) stepsize_ = stepsize;
    else stepsize_ = (upper - lower) / 30.;
  }
}
