#ifndef GUARD_Rdouble_h
#define GUARD_Rdouble_h

#include <vector>
#include "Constants.h"
#include "MesmerConfig.h"
#include "oberror.h"

namespace mesmer
{
class Rdouble;
extern std::vector<Rdouble*> ActiveRdoubles;

class Rdouble
{
private:
  double value, lower, upper, stepsize, prev;
public:
  Rdouble(double val=0.0):value(val),lower(NaN),upper(NaN),stepsize(NaN), prev(NaN){}
  
  operator double() const  { return value; }
  Rdouble& operator = (const double& val);

  void set_range(const double valueL, const double valueU, const double valueS);
  bool get_range(double& lower_, double& upper_, double& stepsize_)const;

  //Returns true if the value is the same as when setUnchanged was last called.
  bool isUnchanged() { return (prev==value); }

  void setUnchanged() { prev = value; }

  // Increment the current value by stepsize if the result will be <= upper
  // and return the result. If not incremented return NaN.
  //This is pre-increment: it is the new value after incrementing that is returned. 
  double& operator++()
  {
    if(value + stepsize > upper)
      return NaN;
    return value += stepsize;
  }
  double& operator--()
  {
    if(value - stepsize < upper)
      return NaN;
    return value -= stepsize;
  }
};

/*
Rdouble is a variable of type double that can hold a range of values,
and the one actually in use currently can be set externally. 

All that is necessary is for the variable type to be changed from
double to Rdouble; it behaves as it did before, like a double.

The range can be set from any external code by setting the variable's
value to NaN. This causes the address of the Rdouble to be pushed
to the quasi-global vector ActiveRdoubles. It is then retrieved by
the code from the back of ActiveRdoubles and set_range called.

  pReact->set_EInf(NaN);
  ActiveRdoubles.back()->set_range(valueL,valueU,valueS);

Alternatively, there is a macro which does some error checking and
which makes a log entry.

Functions like gridSearch() can retrieve the variable's address from
ActiveRdoubles and change its current value to anything. Alternatively,
they can increment or decrement the value using ++ or --, when the end of
the range is indicated by NaN being returned.  
*/

}//namespace
#endif

