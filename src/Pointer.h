//----------------------------------------------------------------------------
// pointer class
// Author: Chi-Hsiu Liang
// A class to cast double values into a dynamic vector for grouped optimization
//----------------------------------------------------------------------------
#ifndef GUARD_Pointer_h
#define GUARD_Pointer_h

#include <vector>
#include <cmath>

using namespace std;

namespace mesmer
{
  class DPoint; // pre-declaration

  extern vector<double> allP;
  extern vector<DPoint> fitDP;

  class DPoint{
  private:
    int    p_;        // index to the double allP vector (if p_ = -1 then this object is not a point)
    double value;     // the value used to calculate the current DOS vector
    double lower;     // lower bound of this variable
    double upper;     // upper bound of this variable
    double stepsize;  // step size used for grid search or fitting

  public:
    //----------------------------------------------
    //Constructors produce an object from ground up
    DPoint(): p_(-1), lower(0.0), upper(0.0), stepsize(0.0) {}
    DPoint(const DPoint& p): p_(p.p_), lower(p.lower), upper(p.upper), stepsize(p.stepsize) {}
    DPoint& operator=(const DPoint& p)
    {
      if (p_ == p.p_) return *this;
      else{ p_ = p.p_; return *this; }
    }

    DPoint(double x) { allP.push_back(x); p_ = int(allP.size()) - 1; }

    DPoint& operator= (const double& r) { allP[p_] =  r;          return *this; }
    DPoint& operator+=(const DPoint& b) { allP[p_] += allP[b.p_]; return *this; }
    DPoint& operator-=(const DPoint& b) { allP[p_] -= allP[b.p_]; return *this; }

    DPoint& operator*=(const double& b) { allP[p_] *= b;          return *this; }
    DPoint& operator/=(const double& b) { allP[p_] /= b;          return *this; }

    const double& get_value(void)       { return allP[p_];                      }
    
    // check if the "value" used to calculate the previous data is identical to the current value in allP vector.
    bool isConstant(void){
      if (value == get_value()){ 
        return true;
      }
      return false;
    }
    
    // This function explicitly tell DPoint to update the value, means that the DOS is re-calculated using the current value.
    void updateValue(void){
      value = get_value();
    }

    void set_range(const double valueL, const double valueU, const double stepsize_);

    void get_range(double& lower_, double& upper_, double& stepsize_);

  };
}

#endif // GUARD_Pointer_h
