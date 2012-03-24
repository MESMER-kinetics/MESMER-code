//-------------------------------------------------------------------------------------------
//
// UnitTests.cpp
//
// Author: Struan Robertson
// Date:   24/Mar/2012
//
// This class implements a number of unit tests. 
//
//-------------------------------------------------------------------------------------------

#include "../System.h"
#include "../calcmethod.h"
#include "../MesmerMath.h"

namespace mesmer
{
  class UnitTests : public CalcMethod
  {
  public:
	UnitTests(const std::string& id) : CalcMethod(id) {}
	virtual ~UnitTests() {}
	//Function to do the work
	virtual bool DoCalculation(System* pSys);

  private:
	bool Test_Chi2SignificanceTest() const;

  };

  ////////////////////////////////////////////////
  //Global instance
  UnitTests theUnitTests("UnitTests");
  ///////////////////////////////////////////////

  bool UnitTests::DoCalculation(System* pSys)
  {

	bool status(true) ; 

	// Test Chi-squared test function.
	status = Test_Chi2SignificanceTest() ;

	return status ;
  }

  bool UnitTests::Test_Chi2SignificanceTest() const {

    size_t nChi2(6) ;  // Number of Ch^2 values.
	size_t nNDOF(10) ; // Number of " Number of Degrees of Freedom" [sic].

    // Test Chi2 values and write header.
	vector<double> Chi2 ;
	for (size_t j(1) ; j <= nChi2 ; j++ ) {
	  double tmp = pow(10.0, (double(j)/2.0) ) ;
	  Chi2.push_back(tmp);
	  ctest << formatFloat(tmp, 6, 15) ;
	}
    ctest << endl ;

 	for (size_t i(0) ; i < nNDOF ; i++ ) {
	  size_t NoDegFreedom((i+1)*10) ;
	  ctest << setw(10) << NoDegFreedom ;
	  for (size_t j(0) ; j < nChi2 ; j++ ) {
		double probChi2 = ChiSquaredPrbFn(Chi2[j]/2.0, double(NoDegFreedom)/2.0) ;
		ctest << formatFloat(probChi2, 6, 15) ;
	  }
	  ctest << endl ;
	}

	return true ;

  }



}//namespace

