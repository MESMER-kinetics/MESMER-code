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
#include <random>

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

  const int nrolls=10000;  // number of experiments
  const int nstars=100;    // maximum number of stars to distribute

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(5.0,2.0);

  int p[10]={};

  for (int i=0; i<nrolls; ++i) {
    double number = distribution(generator);
    if ((number>=0.0)&&(number<=10.0)) ++p[int(number)];
  }

  std::cout << "normal_distribution (5.0,2.0):" << std::endl;

  for (int i=0; i<10; ++i) {
    std::cout << i << "-" << (i+1) << ": ";
    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
  }
    return true;
  }


} //namespace

