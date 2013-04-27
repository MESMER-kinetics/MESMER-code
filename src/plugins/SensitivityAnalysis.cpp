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

    SensitivityAnalysis(const char* id) : m_id(id), m_CovMtx(NULL) { Register(); }

    virtual ~SensitivityAnalysis() { if (m_CovMtx) delete m_CovMtx ; }

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

	virtual const char* getID()  { return m_id; }

  private:

    const char* m_id;

	dMatrix *m_CovMtx ;      // Covariance matrix.

  };

  ////////////////////////////////////////////////
  //Global instance
  SensitivityAnalysis theSensitivityAnalysis("SensitivityAnalysis");
  ///////////////////////////////////////////////

  bool SensitivityAnalysis::DoCalculation(System* pSys)
  {

	// Locate covarinance matrix and, if present, read it in. 

	PersistPtr ppMtx = pSys->getPersistPtr()->XmlMoveTo("me:analysis");
	if (ppMtx) {
	  if (ppMtx = ppMtx->XmlMoveTo("me:hessian")) { 
    	  ppMtx = ppMtx->XmlMoveTo("matrix") ;
	  }
	}

	if (!ppMtx)
	  return false ;

	m_CovMtx = ReadMatrix<double>(ppMtx) ;

	return true;
  }


} //namespace

