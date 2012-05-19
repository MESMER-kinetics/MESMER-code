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
#include "../Spline.h"

namespace mesmer
{
  namespace {

	// An anonymous namespace to containing structure definitions to be used by tests.

	struct chi2Data {
	  chi2Data(double NoDegFreedom, double Chi2, double Gammaq){
		m_NoDegFreedom = NoDegFreedom ;
		m_Chi2         = Chi2 ;
		m_Gammaq       = Gammaq ;
	  }

	  double m_NoDegFreedom ;
	  double m_Chi2 ;
	  double m_Gammaq ;

	} ;

	static const double testCriterion(0.00001) ;
  } 

  class UnitTests : public CalcMethod
  {
  public:
	UnitTests(const std::string& id) : CalcMethod(id) {}
	virtual ~UnitTests() {}
	//Function to do the work
	virtual bool DoCalculation(System* pSys);

  private:

	// Tests:

	bool Test_Chi2SignificanceTest() const;

	bool Test_Spline() const ;

	// Support methods:

	void underlineText(const string& text) const ;

  };

  ////////////////////////////////////////////////
  //Global instance
  UnitTests theUnitTests("UnitTests");
  ///////////////////////////////////////////////

  bool UnitTests::DoCalculation(System* pSys)
  {

	bool status(true) ; 

	ctest << endl ;
	underlineText("MESMER Unit Tests.") ;

	// Test Chi-squared test function.
	status = ( status && Test_Chi2SignificanceTest()) ;

	// Test Spline class.
	status = ( status && Test_Spline()) ;

	ctest << endl ;
	if (status) {
	  ctest << "  All tests pass." ;
	} else {
	  ctest << "  One or more tests failed" ;
	}
	ctest << endl ;

	return status ;
  }

  bool UnitTests::Test_Chi2SignificanceTest() const {

	ctest << endl ;
	underlineText("Test: Chi2 Significance test.") ;

	ctest << endl ;
	ctest << "  Values for the incomplete Gamma function taken from" << endl ;
	ctest << "  M. Abramowitz and I.A. Stegun, Handbook of Mathematical" << endl ;
	ctest << "  Functions, Dover, 1972, pp. 978-983." << endl ;
	ctest << endl ;

	ctest << endl ;
	underlineText("Deg. of F.    Chi2      Published     Calculated") ;
	ctest << endl ;

	vector<chi2Data> TableValues ;

	TableValues.push_back(chi2Data( 5.0,  2.0, 0.84915));
	TableValues.push_back(chi2Data(10.0,  2.0, 0.99634));
	TableValues.push_back(chi2Data(15.0,  2.0, 0.99997));

	TableValues.push_back(chi2Data( 5.0,  4.0, 0.54942));
	TableValues.push_back(chi2Data(10.0,  4.0, 0.94735));
	TableValues.push_back(chi2Data(15.0,  4.0, 0.99774));
	TableValues.push_back(chi2Data(20.0,  4.0, 0.99995));

	TableValues.push_back(chi2Data( 5.0,  6.0, 0.30622));
	TableValues.push_back(chi2Data(10.0,  6.0, 0.81526));
	TableValues.push_back(chi2Data(15.0,  6.0, 0.97975));
	TableValues.push_back(chi2Data(20.0,  6.0, 0.99890));
	TableValues.push_back(chi2Data(25.0,  6.0, 0.99997));

	TableValues.push_back(chi2Data( 5.0,  8.0, 0.15624));
	TableValues.push_back(chi2Data(10.0,  8.0, 0.62884));
	TableValues.push_back(chi2Data(15.0,  8.0, 0.92378));
	TableValues.push_back(chi2Data(20.0,  8.0, 0.99187));
	TableValues.push_back(chi2Data(25.0,  8.0, 0.99949));
	TableValues.push_back(chi2Data(30.0,  8.0, 0.99998));

	TableValues.push_back(chi2Data( 5.0, 10.0, 0.07524));
	TableValues.push_back(chi2Data(10.0, 10.0, 0.44049));
	TableValues.push_back(chi2Data(15.0, 10.0, 0.81974));
	TableValues.push_back(chi2Data(20.0, 10.0, 0.96817));
	TableValues.push_back(chi2Data(25.0, 10.0, 0.99665));
	TableValues.push_back(chi2Data(30.0, 10.0, 0.99977));

	TableValues.push_back(chi2Data( 5.0, 20.0, 0.00125));
	TableValues.push_back(chi2Data(10.0, 20.0, 0.02925));
	TableValues.push_back(chi2Data(15.0, 20.0, 0.17193));
	TableValues.push_back(chi2Data(20.0, 20.0, 0.45793));
	TableValues.push_back(chi2Data(25.0, 20.0, 0.74683));
	TableValues.push_back(chi2Data(30.0, 20.0, 0.91654));

	TableValues.push_back(chi2Data( 5.0, 30.0, 0.00002));
	TableValues.push_back(chi2Data(10.0, 30.0, 0.00086));
	TableValues.push_back(chi2Data(15.0, 30.0, 0.01192));
	TableValues.push_back(chi2Data(20.0, 30.0, 0.06985));
	TableValues.push_back(chi2Data(25.0, 30.0, 0.22429));
	TableValues.push_back(chi2Data(30.0, 30.0, 0.46565));

	TableValues.push_back(chi2Data(10.0, 40.0, 0.00002));
	TableValues.push_back(chi2Data(15.0, 40.0, 0.00045));
	TableValues.push_back(chi2Data(20.0, 40.0, 0.00500));
	TableValues.push_back(chi2Data(25.0, 40.0, 0.02916));
	TableValues.push_back(chi2Data(30.0, 40.0, 0.10486));

	TableValues.push_back(chi2Data(15.0, 50.0, 0.00001));
	TableValues.push_back(chi2Data(20.0, 50.0, 0.00022));
	TableValues.push_back(chi2Data(25.0, 50.0, 0.00213));
	TableValues.push_back(chi2Data(30.0, 50.0, 0.01240));

	TableValues.push_back(chi2Data(20.0, 60.0, 0.00001));
	TableValues.push_back(chi2Data(25.0, 60.0, 0.00011));
	TableValues.push_back(chi2Data(30.0, 60.0, 0.00092));

	bool status(true) ;
	double writeBlankline(0.0) ;
	for (size_t i(0) ; i < TableValues.size() ; i++ ) {
	  double NoDegFreedom = TableValues[i].m_NoDegFreedom ;
	  double Chi2         = TableValues[i].m_Chi2 ;
	  double Gammaq       = TableValues[i].m_Gammaq ;
	  double probChi2     = ChiSquaredPrbFn(Chi2/2.0, double(NoDegFreedom)/2.0) ;
	  if (NoDegFreedom < writeBlankline) {
		ctest << endl ;
	  }
	  writeBlankline = NoDegFreedom ;
	  ctest << formatFloat(NoDegFreedom, 2, 10) ;
	  ctest << formatFloat(Chi2, 2, 10) ;
	  ctest << formatFloat(Gammaq, 5, 15) ;
	  ctest << formatFloat(probChi2, 5, 15) ;
	  if (abs(Gammaq - probChi2) > testCriterion) {
		status = false ;
		ctest << "*";
	  }
	  ctest << endl ;
	}

	return status ;

  }

  bool UnitTests::Test_Spline() const {

	bool status(true) ; 

	ctest << endl ;
	underlineText("Test: Spline class.") ;

	ctest << endl ;
	ctest << "  Values for a cubic are calculated on a grid" << endl ;
	ctest << "  and interpolated using the Spline class." << endl ;
	ctest << endl ;

	const size_t N(10) ;
	const size_t nKnots(2*N+1) ;
	const double width(2.0) ;

	vector<double> x(nKnots, 0.0) ;
	vector<double> y(nKnots, 0.0) ;

	for (size_t i(0) ; i < nKnots ; i++) {
	  double xx = 2.0*width*(double(i)/double(N) - 1.0) ;
	  x[i] = xx ;
	  y[i] = xx * (xx - width) * (xx + width) ;
	}

	Spline spline ;

	double lower = 3.0 *x[0]*x[0]               - width*width ;
	double upper = 3.0 *x[nKnots-1]*x[nKnots-1] - width*width ;
	spline.Initialize(x, y, lower, upper) ;

	ctest << endl ;
	underlineText("      x(grid)        y(grid)      y(intplt)") ;

	for (size_t i(0) ; i < nKnots ; i++) {
	  double xx = width*(double(i)/double(N) - 1.0) + 0.1 ;
	  double yy = xx * (xx - width) * (xx + width) ;
	  double ys = spline.Calculate(xx) ;
	  ctest << formatFloat(xx, 5, 15) ;
	  ctest << formatFloat(yy, 5, 15) ;
	  ctest << formatFloat(ys, 5, 15) ;
	  if (abs(ys - yy) > testCriterion) {
		status = false ;
		ctest << "*";
	  }
	  ctest << endl ;
	}

 	return status ;

 }


  void UnitTests::underlineText(const string& text) const {

	ctest << "  " << text << endl ;
	ctest << "  " ;
	for (size_t i(0) ; i < text.size() ; i++ ) 
	  ctest << "-" ;
	ctest << endl ;

  }

}//namespace

