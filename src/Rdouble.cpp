#include <map>
#include <stdexcept>
#include "Rdouble.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  //static variable
  Rdouble*  Rdouble::pendingVar;
  const double Rdouble::eps = 1e-7;

  Rdouble::operator double() const
  {
    return link ? *link * factor + addand : value;
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

    //Push on to the vector of Rdouble objects which have a range
    withRange().push_back(this);
    //value is NOT set
    lower    = valueL; 
    upper    = valueU;
    stepsize = valueS;
    varname  = txt;
    prev     = NaN;

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
    }
  }


  // This method updates the value of a labelled variable in an XML file
  // with its current Rdouble value.
  void Rdouble::UpdateXMLLabelVariables() {
    labelmap::iterator itrlabel = withLabel().begin() ;
    for ( ;itrlabel != withLabel().end() ; itrlabel++) {
      itrlabel->second->XmlWriteValue() ;
    }
  }

  // Called after all the xml data has been read
  // Look up the linked variable names
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
          cinfo << it->first << " is derived from " << depvar->linkedname << endl;
          //Debug:
          cinfo << it->first << " = " << *(it->second) << endl;
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


  // Utility function to read parameter range. 
  //   The parameter cnvrsnFctr applies any conversion factor that is required.
  //   The parameter rangeSet indicates if a range has actually been set for the rdouble variable.
  //
  bool ReadRdoubleRange(const std::string& name, PersistPtr pp, Rdouble& rdouble, 
    bool& rangeSet, double cnvrsnFctr, double shift)
  {
    // Store the names of all the potential range variables in case they are referred to
    // in a derivedFrom attribute of another linked variable
    rdouble.set_label(name);

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
      rangeSet = false ;
    }

    // If there is a derivedFrom attribute, a pointer to its value is put in
    // 
    const char* pDerivedtxt = pp->XmlReadValue("me:derivedFrom", optional);
    if(pDerivedtxt)
    {
      rdouble.set_link_params(
        pDerivedtxt,
        pp->XmlReadDouble("factor", optional),
        pp->XmlReadDouble("addand", optional));

      cinfo << name << " will be derived from " << pDerivedtxt
            << ". Check below that this has been found." << endl;
    }

/*    // Check for constraint label.

    const char* pLabel = pp->XmlReadValue("label", optional);
    if (pLabel) {
      rdouble.set_label(string(pLabel), pp);
    if (!rangeSet) rdouble.set_varname(name) ;
    }
*/
    // Save a pointer to XML location for result update.
    rdouble.set_XMLPtr(pp) ;

    return true;
  }

}//namespace


