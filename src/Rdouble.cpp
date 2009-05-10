#include "Rdouble.h"

using namespace std;
namespace mesmer
{
  
//static variable
Rdouble*  Rdouble::pendingVar;

Rdouble& Rdouble::operator = (const double& val)
{
  value = val;
  if(IsNan(val))
    pendingVar =this;
  return *this;
}

void Rdouble::set_range(const double valueL, const double valueU, const double valueS, const char* txt)
{
  double numStep = (valueU - valueL)/valueS;
  if(valueU<valueL || valueS<0)
  {
    cerr << "The upper value should be larger than lower and the stepsize should be positive"<<endl;
    return;
  }
  if (numStep < 1)
    cerr << "There is only one point for this variable. Remove range setting to clear this error."<<endl;

  //Push on to the vector of Rdouble objects which have a range 
  withRange().push_back(this);
  value    = lower = valueL;
  upper    = valueU;
  stepsize = valueS;
  prev     = NaN;

  cinfo << txt <<" was given a range with " << (int)numStep+1 << " steps " << endl;
}

void Rdouble::set_range_indirect
                (const double valueL, const double valueU, const double valueS,const char* txt)
{
  if(pendingVar)
    pendingVar->set_range(valueL, valueU, valueS, txt);
  else
    cerr << "Indirect setting of a range " << txt 
         << " without previously set variable to NaN. Range setting ignored." <<endl;
  pendingVar = NULL;
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

