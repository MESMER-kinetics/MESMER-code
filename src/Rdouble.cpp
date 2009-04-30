#include "Rdouble.h"

using namespace std;
namespace mesmer
{

//global variable 
std::vector<Rdouble*> ActiveRdoubles;

Rdouble& Rdouble::operator = (const double& val)
{
  value = val;
  if(IsNan(val))
    ActiveRdoubles.push_back(this);
  return *this;
}

void Rdouble::set_range(const double valueL, const double valueU, const double valueS)
{
  double numStep = (valueU - valueL)/valueS;
  if(valueU<valueL || valueS<0)
  {
    cerr << "The upper value should be larger than lower and the stepsize should be positive"<<endl;
    return;
  }
  cinfo << "Range with " << numStep << " steps set" << endl;
  if (numStep < 1)
    cerr << "There is only one point for this variable. Remove range setting to clear this error."<<endl;

  //Push only when set_range is called directly on a Rdouble object
  //Not needed with the setting NaN route.
  if(ActiveRdoubles.back()!=this)
    ActiveRdoubles.push_back(this);
  value    = lower = valueL;
  upper    = valueU;
  stepsize = valueS;
  prev     = NaN;
}

bool Rdouble::get_range(double& lower_, double& upper_, double& stepsize_)const
{
  if(IsNan(lower))
    return false;
  lower_   = lower;
  upper_   = upper;
  stepsize_= stepsize;
  return true;
}

}//namespace