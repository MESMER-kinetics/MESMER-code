//-------------------------------------------------------------------------------------------
//
// SensitivityAnalysis.cpp
//
// Author: Struan Robertson
// Date:   25/Nov/2012
//
// This class implements sensitivity analysis algorithm used to identify this most important
// parameters from a fit. 
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"

namespace mesmer
{
  class SensitivityAnalysis : public CalcMethod
  {
  public:

    SensitivityAnalysis(const char* id) : m_id(id) {}

    virtual ~SensitivityAnalysis() {}

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

	virtual const char* getID()  { return m_id; }

  private:

    const char* m_id;
  };

  ////////////////////////////////////////////////
  //Global instance
  SensitivityAnalysis theSensitivityAnalysis("SensitivityAnalysis");
  ///////////////////////////////////////////////

  bool SensitivityAnalysis::DoCalculation(System* pSys)
  {

    return true;
  }


} //namespace

