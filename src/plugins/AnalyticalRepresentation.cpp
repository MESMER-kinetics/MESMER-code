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

    AnalyticalRepresentation(const char* id) : FittingUtils(), m_id(id)
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

	m_format      = pp->XmlReadValue("me:analyticalRepType");
    m_reactionRef = pp->XmlReadValue("me:analyticalRepRef");

    // Read in parameter constraints.
    //ReadParameterConstraints(ppControl) ;

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

	double chiSquare(0.0) ;

	// First gets some points.

	vector<double> Temperature;
	vector<double> Concentration;
	vector<qdMatrix *> RateCoefficients;

    pSys->calculate(Temperature, Concentration, RateCoefficients) ;

	// Release rate coefficient data.

	for (size_t i(0) ; i < RateCoefficients.size() ; i++) {
	  delete RateCoefficients[i] ;
	}

    return true;
  }


}//namespace

