//-------------------------------------------------------------------------------------------
//
// AnalyticalRepresentation.cpp
//
// Author: Struan Robertson
// Date:   21/Jul/2013
//
// This class implements the methods to that determinne and analytical representation of a 
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
  class AnalyticalReperesentation : public CalcMethod, private FittingUtils
  {
  public:

    AnalyticalReperesentation(const char* id) : FittingUtils(), m_id(id)
    { Register(); }

    virtual ~AnalyticalReperesentation() {}
    virtual const char* getID()  { return m_id; }

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  protected:

  private:

    const char* m_id;

  };

  ////////////////////////////////////////////////
  //Global instance
  AnalyticalReperesentation theAnalyticalReperesentation("Analytical Reperesentation");
  ///////////////////////////////////////////////

  bool AnalyticalReperesentation::DoCalculation(System* pSys)
  {
    //Read in fitting parameters, or use values from defaults.xml.
    PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");
    unsigned maxIterations= ppControl->XmlReadInteger("me:fittingIterations");
    double tol = ppControl->XmlReadDouble("me:fittingTolerance");

	// Read in parameter constraints.
//	ReadParameterConstraints(ppControl) ;

    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are caluclated).
    pSys->m_Flags.useTheSameCellNumber = true;

    // Warnings and less not sent to console.
    ChangeErrorLevel e(obError); 

	double chiSquare(0.0) ;

    pSys->calculate(chiSquare) ;

    return true;
  }


}//namespace

