//-------------------------------------------------------------------------------------------
//
// AnalyticalRepresentation.cpp
//
// Author: Struan Robertson
// Date:   21/Jul/2013
//
// This class implements the methods to that determine and analytical representation of a 
// master equation representation.
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"
#include "../dMatrix.h"
#include "FittingUtils.h"

namespace mesmer
{
  class AnalyticalRepresentation : public CalcMethod, private FittingUtils
  {
  public:

    AnalyticalRepresentation(const char* id) : FittingUtils(), 
      m_id(id),
      m_format(),
      m_reactionRef(),
      m_TMax(0.0),
      m_TMin(0.0),
      m_CMax(0.0),
      m_CMin(0.0),
      m_NTpt(0),
      m_NCpt(0)
    { Register(); }

    virtual ~AnalyticalRepresentation() {}
    virtual const char* getID()  { return m_id; }

    //Read in data for this method from XML
    virtual bool ParseData(PersistPtr pp);

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  protected:

  private:

    const char* m_id;
    string m_format ;
    string m_reactionRef ;

    double m_TMax ;
    double m_TMin ;
    double m_CMax ;
    double m_CMin ;

    double m_TSft ;
    double m_TDiv ;
    double m_CSft ;
    double m_CDiv ;

    int    m_NTpt ;
    int    m_NCpt ;

  };

  ////////////////////////////////////////////////
  //Global instance
  AnalyticalRepresentation theAnalyticalRepresentation("analyticalRepresentation");
  ///////////////////////////////////////////////

  bool AnalyticalRepresentation::ParseData(PersistPtr pp)
  {
    /* Typical data
    <me:control>
    ...
    <me:calcMethod xsi:type="me:analyticalRepresentation">
    <me:fittingTolerance>0.1</me:fittingTolerance>
    <me:fittingIterations>5</me:fittingIterations>
    </me:calcMethod>
    ...
    </me:control>
    */
    //Read in fitting parameters, or use values from defaults.xml.

    m_reactionRef = pp->XmlReadValue("me:analyticalRepRef");

    m_NTpt = pp->XmlReadInteger("me:chebNumTemp");
    m_NCpt = pp->XmlReadInteger("me:chebNumConc");

    m_TMax = pp->XmlReadDouble("me:chebMaxTemp");
    m_TMin = pp->XmlReadDouble("me:chebMinTemp");
    m_CMax = pp->XmlReadDouble("me:chebMaxConc");
    m_CMin = pp->XmlReadDouble("me:chebMinConc");

    double RpTMin = 1.0 / m_TMin ;
    double RpTMax = 1.0 / m_TMax ;
    double lgCMin = log10(m_CMin);
    double lgCMax = log10(m_CMax);

    m_TSft = - RpTMin - RpTMax;
    m_TDiv = 1.0 / (RpTMax - RpTMin);
    m_CSft = - lgCMin - lgCMax;
    m_CDiv = 1.0 / (lgCMax - lgCMin);

    return true;
  }

  bool AnalyticalRepresentation::DoCalculation(System* pSys)
  {
    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are caluclated).
    pSys->m_Flags.useTheSameCellNumber = true;

    // Warnings and less not sent to console.
    ChangeErrorLevel e(obError); 

    // First gets some points.

    vector<double> Temperature;
    vector<double> Concentration;
    vector<double> RateCoefficients;

    pSys->calculate(Temperature, Concentration, RateCoefficients) ;

    return true;
  }


}//namespace

