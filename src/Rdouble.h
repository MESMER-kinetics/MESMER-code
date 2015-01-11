#ifndef GUARD_Rdouble_h
#define GUARD_Rdouble_h

#include <vector>
#include <string>
#include "Constants.h"
#include "MesmerConfig.h"
#include "error.h"
#include "Persistence.h"

namespace mesmer
{

  class Rdouble
  {
  private:
    double value, lower, upper, stepsize, prev;
	double m_unitConversion ;
    std::string varname ;
    PersistPtr m_XMLPtr ;
    //linked variable
    const char* linkedname;
    Rdouble* link; //NULL if this is an independent variable
    double factor, addand;

    static Rdouble* pendingVar;
    static const double eps;
  public:
    Rdouble(double val=0.0):value(val), lower(NaN), upper(NaN), stepsize(NaN), prev(NaN), m_unitConversion(1.0),
      varname(), m_XMLPtr(NULL), linkedname(NULL), link(NULL), factor(1.0), addand(0.0){}

    operator double() const;
    Rdouble& operator = (const double& val);

    double set_to_lower() { return value = lower; }
    void test(){return;}

    // Vector of Rdoubles that have a range. 
    // A function, rather than static member variable, for proper initialization.   
    static std::vector<Rdouble*>& withRange()
    {
      static std::vector<Rdouble*> wr;
      return wr;
    }

    void set_range(const double valueL, const double valueU, const double valueS,const char* txt=NULL);
    static void set_range_indirect
      (const double valueL, const double valueU, const double valueS,const char* txt=NULL);

    bool get_range(double& lower_, double& upper_, double& stepsize_)const;

    const char* get_varname(){ return varname.c_str(); }
	double get_lower(){ return lower; }
	double get_upper(){ return upper; }
	double get_stepsize(){ return stepsize; }

    void set_varname(const std::string& name) { varname = name; }

    int get_numsteps(){ return (int)(1 + eps + (upper - lower) / stepsize); } 

    // Map of Rdouble names
    // c.f. withRange above. 
    typedef std::map<std::string, Rdouble* > labelmap ; 
    static labelmap& withLabel()
    {
      static labelmap lr;
      return lr;
    }

    void set_label(const std::string& label) ;
    void set_link_params(const char* name, double fac, double add)
    {
      linkedname = name;
      if(!IsNan(fac))
        factor = fac;
      if(!IsNan(add))
        addand = add;
    }

    static bool SetUpLinkedVars();

    void set_XMLPtr(PersistPtr pp) { m_XMLPtr = pp ;}

    void XmlWriteValue() {
      std::ostringstream s; 
      s << *(this) ;
      m_XMLPtr->XmlWrite(s.str());
    }

	void store_conversion(double conversion_factor) { m_unitConversion = conversion_factor ; }

	double originalUnits() const { return value/((m_unitConversion > 0.0) ? m_unitConversion : 1.0) ; } ;

	double originalUnits(double quantity) const { return quantity/((m_unitConversion > 0.0) ? m_unitConversion : 1.0) ; } ;

	void XmlWriteAttribute(const std::string& name, const std::string& value) 
    {
      m_XMLPtr->XmlWriteAttribute(name, value) ;
    }
     
    static void UpdateXMLLabelVariables() ;

    //Returns true if the value is the same as when setUnchanged was last called.
    bool isUnchanged() { return (prev==value); }

    void setUnchanged() { prev = value; }

    // Increment the current value by stepsize if the result will be <= upper
    // and return the result. If not incremented return NaN and the value is reset to lower.
    //This is pre-increment: it is the new value after incrementing that is returned. 
    const double& operator++()
    {
      value += stepsize;
      if(value-upper < eps*stepsize)
        return value; 
      value = lower;
      return NaN;
    }
    const double& operator--()
    {
      value -= stepsize;
      if(lower-value < eps*stepsize)
        return value; 
      value = upper;
      return NaN;
    }
  };

  /*
  Rdouble is a variable which looks like a double but can hold a range
  of values, and the one actually in use currently can be controlled externally. 

  All that is necessary is for the variable type to be changed from
  double to Rdouble, At this stage it behaves as it did before, like a double.

  To give it a range, and if there is access to the variable, use set_range().
  If there is not, for example, if a member variable of Reaction is being set
  from an ILT class, the trick used is to set the value to NaN, which stores
  its address in a static variable. The static function Rdouble::set_rangeindirect()
  is then called to set the range, e.g.,

  pReact->set_EInf(NaN);
  Rdouble::set_range_indirect(valueL, valueU, valueS, "EInf");

  Both the set_range functions have an optional text parameter which is
  intended to contain the name of the variable and is used in a log message.

  RDouble::withRange() returns a reference to the std::vector<Rdouble*>
  that contains pointers to all the Rdoubles with a range set.
  This allows the current value of any of these to be set to anything.
  Alternatively, an Rdouble can be incremented using ++. Normally the return
  is the value after incrementing but at the end of the range NaN is returned
  and the value is reset to the lower value.

  A Rdouble variable can be linked to the value of another Rdouble variable
  by adding a me:derivedFrom attribute to its definition, e.g.

  <molecule id="IM2">
 ...
    <property dictRef="me:deltaEDown">
     <scalar me:derivedFrom="IM1:DeltaEDown">161</scalar>
    </property>
 ...
  </molecule>
  
  IM2:DeltaEDown will always have the same value as IM1:DeltaEDown (the 161 value
  is ignored but must be present). Optionally a factor and addand attribute can
  be included: 
   <scalar derivedFrom="IM1:DeltaEDown" factor="1.0 addand="0.0>161</scalar>
  when (in this case)
   IM2:DeltaEDown = IM1:DeltaEDown * factor  +  (addand)

  */

  std::istringstream& operator>>(std::istringstream& iss, Rdouble& rdouble) ;

  // Utility function to read parameter range. 
  //
  bool ReadRdoubleRange(const std::string& name, PersistPtr pp, Rdouble& rdouble, 
    bool& rangeSet, double cnvrsnFctr = 1.0, double shift = 0.0) ;

}//namespace
#endif

