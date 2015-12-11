#include <map>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include "unitsConversion.h"
#include "Rdouble.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  //static variables
  Rdouble*  Rdouble::pendingVar;
  const double Rdouble::eps = 1e-7;

  Rdouble::operator double() const
  {
    return value;
  }

  Rdouble& Rdouble::operator = (const double& val)
  {
    prev = value ;
    value = val ;

    if(IsNan(val))
      pendingVar =this;
    return *this;
  }

  void Rdouble::set_range(const double valueL, const double valueU, const double valueS, const char* txt)
  {
    if(valueU<valueL || valueS<0)
    {
      cerr << "in " << txt << " The upper value should be larger than lower and the stepsize should be positive"<<endl;
      return;
    }

    //On a ReadRdoubleRange() call, set_label always adds to withLabel.
    //If there is a range specified (lower, upper, etc.), set_range adds to withRange.
    //value is NOT set
    lower    = valueL; 
    upper    = valueU;
    stepsize = valueS;
    varname  = txt;
    prev     = NaN;

    if(!IsNan(lower))
      //Push on to the vector of Rdouble objects with a valid range
      withRange().push_back(this);

    if (get_numsteps() <= 1)
      cwarn << "There is only one point for the range variable " << txt << endl;

    cinfo << txt <<" was given a range with " << get_numsteps() << " steps " << endl;
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


  void Rdouble::set_label(const std::string& label) {

    // Check to see if label already used.
    labelmap::iterator itrlabel = withLabel().find(label) ;

    if (itrlabel != withLabel().end()) {
      cerr << "Error: Parameter label " << label << " redefined." << endl;
    } else {
      withLabel()[label] = this ;
      varname = label;
    }
  }

  // Calculate the current value of a derivedFrom variable in cm-1
  void Rdouble::update_value()
  {
    value = *link * factor + addand;
  }

  // This static method updates the value of all the derivedFrom variables in an XML file
  // with their current Rdouble value.
  void Rdouble::UpdateXMLDerivedVariables() {
    labelmap::iterator itrlabel = withLabel().begin() ;
    for ( ;itrlabel != withLabel().end() ; itrlabel++) {
      if (itrlabel->second->link)
      {
        itrlabel->second->update_value();
        itrlabel->second->XmlWriteValue();
      }
    }
  }

   // Called after all the xml data has been read
  bool Rdouble::SetUpLinkedVars()
  {
    bool ok=true;
    labelmap::iterator it;
    for(it=withLabel().begin();it!=withLabel().end();++it)
    {
      Rdouble* depvar = it->second;
      if(depvar->linkedname)
      {
        labelmap::iterator pos = withLabel().find(depvar->linkedname);
        if(pos!=withLabel().end())
        {
          depvar->link = pos->second; //pointer to independent variable
          cinfo << it->first << " is derived from " << depvar->linkedname;
          if (depvar->link && IsNan(depvar->link->lower)) //i.e. NOT a range variable
          {
            depvar->m_unitConversion = depvar->link->m_unitConversion;
            depvar->update_value();
            cinfo << ", and has a fixed value = " << depvar->value/ depvar->m_unitConversion
                  << " " << depvar->m_units;
          }
          cinfo << endl;
         
          //The units and conversion factor must be the same for both variables; checked below.
          depvar->m_unitConversion = depvar->link->m_unitConversion;

          string depunits = depvar->m_units;
          if (depunits.empty())
            depvar->m_units = depvar->link->m_units;
          else
            //check units are the same
            if (depvar->m_units != depvar->link->m_units)
            {
              cerr << depvar->varname << " and " << depvar->linkedname
                << " must have the same units, or both must have no units." << endl;
              return false;
            }

          // Convert addand of dependent variable to the units of the independent variable
          // Currently this works only with energy units (ZPE, deltaEDown),
          // but could be generalised.
          if (!depvar->m_units.empty() && !depvar->link->m_units.empty())
            depvar->set_addand(ConvertEnergy(depvar->m_units, depvar->link->m_units, depvar->addand));
        }
        else
        {
          cerr << "In linking to " <<  it->first << ", " 
               <<depvar->linkedname << " could not be found" << endl;
          ok = false;
        }
      }
    }
    return ok;
  }

  bool Rdouble::set_link_params(const char* name, double fac, double add, const char* units)
  {
    linkedname = name;
    if (!IsNan(fac))
      factor = fac;
    if (!IsNan(add))
      addand = add * m_unitConversion;
    if (units)
      m_units = std::string(units);
    return true;
  }

  // Utility function to read parameter range. 
  //   The parameter cnvrsnFctr applies any conversion factor that is required.
  //   The parameter rangeSet indicates if a range has actually been set for the rdouble variable.
  //   If there is a derivedFrom attribute, the cnvrsnFctr and shift parameters are ignored.
  //
  bool ReadRdoubleRange(const std::string& name, PersistPtr pp, Rdouble& rdouble, 
    bool& rangeSet, double cnvrsnFctr, double shift)
  {
    // Store the names of all the potential range variables, their units and values,
    // in case they are referred to in a derivedFrom attribute of another linked variable.
    rdouble.set_label(name);
    rdouble.set_units(pp->XmlReadValue("units", optional));

    // Save a pointer to XML location for result update.
    rdouble.set_XMLPtr(pp) ;
    rdouble.store_conversion(cnvrsnFctr);

    const char* pLowertxt = pp->XmlReadValue("lower", optional);
    const char* pUppertxt = pp->XmlReadValue("upper", optional);
    const char* pStepStxt = pp->XmlReadValue("stepsize", optional);

    if (pLowertxt && pUppertxt){
      rangeSet = true ;
      double valueL(0.0), valueU(0.0), stepsize(0.0);
      stringstream strLower(pLowertxt), strUpper(pUppertxt), strStepSize(pStepStxt);
      strLower    >> valueL; 
      strUpper    >> valueU; 
      strStepSize >> stepsize;
      valueL   = cnvrsnFctr*valueL   + shift; 
      valueU   = cnvrsnFctr*valueU   + shift; 
      stepsize = cnvrsnFctr*stepsize + shift;


      rdouble.set_range(valueL, valueU, stepsize, name.c_str());	  
	 
    } else {
      rdouble.set_value(); //reads and stores converted value
      rangeSet = false ;
    }

    // If there is a derivedFrom attribute, a pointer to its value is put in
    const char* pDerivedtxt = pp->XmlReadValue("me:derivedFrom", optional);
    if(!pDerivedtxt)
      pDerivedtxt = pp->XmlReadValue("derivedFrom", optional);
    if(pDerivedtxt)
    {
      rdouble.set_link_params(
        pDerivedtxt,
        pp->XmlReadDouble("factor", optional),
        pp->XmlReadDouble("addand", optional),
        pp->XmlReadValue("units", optional));

      cinfo << name << " will be derived from " << pDerivedtxt
            << ". Check below that this has been found." << endl;
    }
    return true;
  }

  std::istringstream& operator>>(std::istringstream& iss, Rdouble& rdouble) {

	double value(0.0);
	iss >> value;
	rdouble = value;

	return iss ;
  }

}//namespace


