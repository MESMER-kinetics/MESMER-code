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
#include "../dMatrix.h"
#include "../Sobol.h"
#include "../qd_test.h"

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
    UnitTests(const char* id) : m_id(id) { Register(); }
    virtual ~UnitTests() {}
    virtual const char* getID()  { return m_id; }
    
    virtual bool DoesOwnParsing() { return true; }

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:

    // Tests:

    bool Test_Chi2SignificanceTest() const;

    bool Test_Spline() const ;

    bool Test_Sobol() const ;

    template<class T> 
    bool Test_LinearAlgebra(string precision) const ;

    // MEIC tests:

    bool Test_MEIC_1(Molecule* pMol) const ;

    bool Test_MEIC_Anharmonic(Molecule *pMol) const ;

    bool Test_MEIC_Rotors(Molecule *pMol) const ;

    bool Test_MEIC_HinderedRotor(Molecule *pMol) const ;

    // Support methods:

    bool Test_MEIC_formGrainSOS(vector<double> &cellDOS) const ;

    bool parseInput(System* pSys) const ;

    void underlineText(const string& text) const ;

  private:

    const char* m_id;

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

    // Parse input.
    parseInput(pSys) ;

    // Test Chi-squared test function.
    status = ( status && Test_Chi2SignificanceTest()) ;

    // Test Spline class.
    status = ( status && Test_Spline()) ;

    // Test Sobol random number class.
    status = ( status && Test_Sobol()) ;

    // Test Linear Algebra ;
    status = ( status && Test_LinearAlgebra<double>(string("double"))) ;
    status = ( status && Test_LinearAlgebra<long double>(string("long-double"))) ;
    status = ( status && Test_LinearAlgebra<dd_real>(string("double-double"))) ;
    status = ( status && Test_LinearAlgebra<qd_real>(string("quad-double"))) ;

    MesmerEnv& Env = pSys->getEnv() ;
    Env.GrainSize  = 100 ; 
    Env.MaxGrn     = 500 ;
    Env.MaxCell    = Env.GrainSize * Env.MaxGrn ;

    // MEIC test: Harmonic oscillator.
    Molecule *pMol = pSys->getMoleculeManager()->find("Test_Molecule");
    if(!pMol)
      pMol = pSys->getMoleculeManager()->find("Test Molecule");
    status = ( status && Test_MEIC_1(pMol)) ;

    pMol = pSys->getMoleculeManager()->find("AcO2_Harmonic");
    if(!pMol)
      pMol = pSys->getMoleculeManager()->find("AcO2 Harmonic");
    status = ( status && Test_MEIC_1(pMol)) ;

    pMol = pSys->getMoleculeManager()->find("AcO2_Anharmonic");
    if(!pMol)
      pMol = pSys->getMoleculeManager()->find("AcO2 Anharmonic");
    status = ( status && Test_MEIC_Anharmonic(pMol)) ;

    status = ( status && Test_MEIC_Rotors(pMol)) ;

    pMol = pSys->getMoleculeManager()->find("AcO2_Asymmetric_Top");
    if(!pMol)
      pMol = pSys->getMoleculeManager()->find("AcO2 Asymmetric Top");
    status = ( status && Test_MEIC_Rotors(pMol)) ;

    pMol = pSys->getMoleculeManager()->find("Ethane");
    status = ( status && Test_MEIC_HinderedRotor(pMol)) ;

    // Execute QD tests for extended precision.

    ctest << endl ;
    underlineText("QD Unit Tests.") ;

    ctest << endl ;
    underlineText("double-double tests:") ;

    TestSuite<dd_real> *testdd = new TestSuite<dd_real> ;

    status = ( status && testdd->testall() ) ;

    delete testdd ;

    ctest << endl ;
    underlineText("quad-double tests:") ;

    TestSuite<qd_real> *testqd = new TestSuite<qd_real>  ;

    status = ( status && testqd->testall() ) ;

    delete testqd ;

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

  bool UnitTests::Test_Sobol() const {

    bool status(true) ; 

    ctest << endl ;
    underlineText("Test: Sobol random number class.") ;

	Sobol sobol ;
	const size_t rnsize(3) ;

	ctest << endl ;
    underlineText("Sobol: single precision") ;

	vector<float> rndmf(rnsize,0.0) ;
	int seed(0) ;
    for (size_t i(0) ; i < 20 ; i++) {
      sobol.sobol(rndmf.size(), &seed, rndmf) ;
      for (size_t j(0) ; j < rndmf.size() ; j++) {
	    ctest << formatFloat(rndmf[j], 5, 15) ;
	  }
	  ctest << endl ;
	}

	ctest << endl ;
    underlineText("Uniform: single precision") ;

	seed = 1;
    for (size_t i(0) ; i < 20 ; i++) {
      ctest << "   " << sobol.i4_uniform(0, 1000000, &seed) << endl ;
	}

	ctest << endl ;
    underlineText("Sobol: double precision") ;

	vector<double> rndmd(rnsize,0.0) ;
	long long seed2(0) ;
    for (size_t i(0) ; i < 20 ; i++) {
      sobol.sobol(rndmd.size(), &seed2, rndmd) ;
      for (size_t j(0) ; j < rndmd.size() ; j++) {
	    ctest << formatFloat(rndmd[j], 5, 15) ;
	  }
	  ctest << endl ;
	}

	ctest << endl ;
    underlineText("Uniform: double precision") ;

	seed = 1;
    for (size_t i(0) ; i < 20 ; i++) {
      ctest << "   " << sobol.i8_uniform(0, 1000000, &seed) << endl ;
	}

	ctest << endl ;

    return status ;

  }

  template<class T> 
  bool UnitTests::Test_LinearAlgebra(string precision) const {

    bool status(true) ; 

    ctest << endl ;
    const string Title = string("Test: Linear Algebra: ") + precision ; 
    underlineText(Title) ;

    // Equilibrium matrix.

    TMatrix<T> Mtx1(3,0.0) ;
    Mtx1[0][0] = -3.0/2.0 ;
    Mtx1[0][1] = 1.0 ;
    Mtx1[1][1] = -10.0/6.0 ;
    Mtx1[1][2] = 1.0 ;
    Mtx1[2][0] = 1.0 ;
    Mtx1[2][1] = 1.0 ;
    Mtx1[2][2] = 1.0 ;

    string Heading("Test matrix 1") ;
    Mtx1.print(Heading, ctest) ;

    Mtx1.invertGaussianJordan() ;

    Heading.clear() ;
    Heading = "Inverse of Test matrix 1" ;
    Mtx1.print(Heading, ctest) ;

    // Calculate Equilibrum population.

    size_t msize = Mtx1.size() ;
    vector<T> rhs(msize,0.0) ;
    rhs[2] = 1.0 ;
    rhs *= Mtx1 ;

    // Rate coefficient matrix.

    TMatrix<T> Mtx2(msize,0.0) ;

    Mtx2[0][0] = -13.0 ;
    Mtx2[0][1] =   2.0 ;
    Mtx2[0][2] =   4.0 ;
    Mtx2[1][0] =   3.0 ;
    Mtx2[1][1] = -12.0 ;
    Mtx2[1][2] =   6.0 ;
    Mtx2[2][0] =  10.0 ;
    Mtx2[2][1] =  10.0 ;
    Mtx2[2][2] = -10.0 ;

    Heading.clear() ;
    Heading = "Test matrix 2" ;
    Mtx2.print(Heading, ctest) ;

    // Symmetrize rate coefficient matrix.

    TMatrix<T> Mtx3(msize, 0.0) ;
    TMatrix<T> Mtx4(Mtx3) ;
    for (size_t i(0) ; i < msize ; i++) {
      Mtx3[i][i] = sqrt(rhs[i]) ;
      Mtx4[i][i] = 1.0/Mtx3[i][i] ;
      rhs[i] = 0.0 ;
    }

    TMatrix<T> Mtx5 = Mtx4*Mtx2*Mtx3 ;

    Heading.clear() ;
    Heading = "Symmetrized test matrix 2" ;
    Mtx5.print(Heading, ctest) ;

    // Diagonalize symmetrize rate coefficient matrix.

    for (size_t i(0) ; i < msize ; i++) {
      for (size_t j(i+1) ; j < msize ; j++) {
        Mtx3[i][j] = Mtx3[j][i] ;
      }
    }

    Mtx5.diagonalize(&rhs[0]) ;

    Heading.clear() ;
    Heading = "Eigenvectors of symmetrized test matrix 2" ;
    Mtx5.print(Heading, ctest) ;

    ctest << endl ;
    underlineText("Eigenvalues:") ;
    for (size_t i(0) ; i < msize ; i++) {
      ctest << formatFloat(rhs[i], 5, 15) << endl ;
    }

    TMatrix<T> Mtx6 = Mtx3*Mtx5 ;
    Heading.clear() ;
    Heading = "Right eigenvectors of test matrix 2" ;
    Mtx6.print(Heading, ctest) ;

    Mtx5.Transpose() ;
    TMatrix<T> Mtx7 = Mtx5*Mtx4 ;
    Heading.clear() ;
    Heading = "Left eigenvectors of test matrix 2" ;
    Mtx7.print(Heading, ctest) ;

    TMatrix<T> Mtx8 = Mtx6*Mtx7 ;
    Heading.clear() ;
    Heading = "Check: Left * Right = Identity" ;
    Mtx8.print(Heading, ctest) ;

    return status ;

  }

  bool UnitTests::Test_MEIC_1(Molecule *pMol) const {

    bool status(true) ; 

    // First, extract the total, i.e. rovibrational, densities of states.

    ctest << endl ;

    vector<double> totalCellDOS ;
    pMol->getDOS().getCellDensityOfStates(totalCellDOS, 0, false) ;

    underlineText(string("MEIC Test: \"") + pMol->getName() + string("\" Harmonic oscillators + rotors.") ) ;

    status = status && Test_MEIC_formGrainSOS(totalCellDOS) ;

    ctest << endl ;

    underlineText(string("MEIC Test: \"") + pMol->getName() + string("\" Harmonic oscillators only.") ) ;

    DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find("BeyerSwinehart");

    // Calculate vibrational densities of states.

    size_t MaximumCell(50000) ;
    vector<double> cellDOS(MaximumCell, 0.0) ;
    cellDOS[0] = 1.0 ;
    pMol->getDOS().setCellDensityOfStates(cellDOS) ; 

    status = status && pDOSCalculator->countCellDOS(&(pMol->getDOS()), MaximumCell) ;

    // Retrieve the DOS vector without recalculating.

    pMol->getDOS().getCellDensityOfStates(cellDOS, 0, false) ;

    return (status && Test_MEIC_formGrainSOS(cellDOS)) ;

  }

  bool UnitTests::Test_MEIC_Anharmonic(Molecule *pMol) const {

    bool status(true) ; 

    // First, extract the total, i.e. rovibrational, densities of states.

    ctest << endl ;

    vector<double> totalCellDOS ;
    pMol->getDOS().getCellDensityOfStates(totalCellDOS, 0, false) ;

    underlineText(string("MEIC Test: \"") + pMol->getName() + string("\" Anharmonic oscillators + rotors.") ) ;

    status = status && Test_MEIC_formGrainSOS(totalCellDOS) ;

    ctest << endl ;
    underlineText(string("MEIC Test: \"") + pMol->getName() + string("\" Anharmonic oscillators only.") ) ;

    size_t MaximumCell(50000) ;
    vector<double> cellDOS(MaximumCell, 0.0) ;
    cellDOS[0] = 1.0 ;
    pMol->getDOS().setCellDensityOfStates(cellDOS) ;

    //Use the instance of Morse plugin initialized earlier.
    DensityOfStatesCalculator* pDOSCalculator = pMol->getDOS().GetDOSCalculator("Morse");
    if(!pDOSCalculator)
      return false;

    status = status && pDOSCalculator->countCellDOS(&(pMol->getDOS()), MaximumCell) ;

    // Retrieve the DOS vector without recalculating.

    pMol->getDOS().getCellDensityOfStates(cellDOS, 0, false) ;

    return Test_MEIC_formGrainSOS(cellDOS) ;

  }

  bool UnitTests::Test_MEIC_Rotors(Molecule *pMol) const {

    bool status(true) ; 

    ctest << endl ;
    underlineText(string("MEIC Test: \"") + pMol->getName() + string("\" Rotors only.") ) ;

    DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find("QMRotors");

    // Calculate rotational densities of states.

    size_t MaximumCell(50000) ;
    vector<double> cellDOS(MaximumCell, 0.0) ;
    pMol->getDOS().setCellDensityOfStates(cellDOS) ;

    status = status && pDOSCalculator->countCellDOS(&(pMol->getDOS()), MaximumCell) ;

    // Retrieve the DOS vector without recalculating.

    pMol->getDOS().getCellDensityOfStates(cellDOS, 0, false) ;

    return Test_MEIC_formGrainSOS(cellDOS) ;

  }

  bool UnitTests::Test_MEIC_HinderedRotor(Molecule *pMol) const {

    bool status(true) ; 

    ctest << endl ;
    underlineText(string("MEIC Test: \"") + pMol->getName() + string("\" Hindered Rotor only.") ) ;

    ctest << endl ;
    ctest << "  State values should be compared with:" << endl ;
    ctest << "  Waage & Rabinovitch, IJCK, 3, 105 (1971)" << endl ;
    ctest << endl ;
	
    DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find("HinderedRotorQM1D");

    // Calculate internal rotational densities of states.

    size_t MaximumCell(50000) ;
    vector<double> cellDOS(MaximumCell, 0.0) ;
	cellDOS[0] = 1.0 ;
    pMol->getDOS().setCellDensityOfStates(cellDOS) ;

    PersistPtr ppMol = pMol->get_PersistentPointer() ;
	ppMol = ppMol->XmlMoveTo("me:ExtraDOSCMethod") ;
    status = status && pDOSCalculator->ReadParameters(&(pMol->getDOS()), ppMol) ;
    status = status && pDOSCalculator->countCellDOS(&(pMol->getDOS()), MaximumCell) ;

    return status ;

  }

  // Support methods.
  //----------------------------------------------------------

  bool UnitTests::Test_MEIC_formGrainSOS(vector<double> &cellDOS) const {

    bool status(true) ; 

    // Calculate grain numbers and averages.

    const int GrainSize    = 10 ;
    const int MaximumGrain = cellDOS.size()/GrainSize ;
    vector<double> cellEne(cellDOS.size(), 0.0) ; 
    vector<double> grainDOS(MaximumGrain, 0.0) ;
    vector<double> grainEne(MaximumGrain, 0.0) ;
    size_t i(0) ;
    for (i = 0 ; i < cellDOS.size() ; i++) {
      cellEne[i] += double(i) ;
    }

    calcGrainAverages(MaximumGrain, GrainSize, cellDOS, cellEne, grainDOS, grainEne) ;

    // Calculate cell and grain sums of states.

    for (i = 1 ; i < cellDOS.size() ; i++) {
      cellDOS[i] += cellDOS[i-1] ;
    }

    for (i = 1 ; i < grainDOS.size() ; i++) {
      grainDOS[i] += grainDOS[i-1] ;
    }

    // Now write out results. The selected energies are those defined by the MEIC test.

    ctest << endl ;
    underlineText("   Energy/cm-1                  SoS Cell                 SoS Grain") ;

    ctest << formatFloat(cellEne[0],   5, 15) << "," ;
    ctest << formatFloat(cellDOS[0],  13, 25) << "," ;
    ctest << formatFloat(grainDOS[0], 13, 25) ;
    ctest << endl ;

    const double tolerance = 0.05 ;
    size_t idx(0), jdx(0) ;
    for (i = 1 ; i < 40 ; i++) {
      idx += (idx < 1000) ? 100 : 1000  ;
      const double tcDOS = cellDOS[idx] ;
      double tgDOS(0.0) ;
      const double energy = cellEne[idx] ;

      // Because grains with no content are elimated there is no simply
      // mapping between cells and grains. Consequently the grain whose 
      // mean energy does not exceed that of the specified energy is used.
      // The following loop searches for that grain.

      while (grainEne[jdx] <= energy ) {
        tgDOS = grainDOS[jdx] ;
        jdx++ ;
      }
      ctest << formatFloat(energy, 5, 15) << "," ;
      ctest << formatFloat(tcDOS, 13, 25) << "," ;
      ctest << formatFloat(tgDOS, 13, 25) ;
      if (abs(tcDOS - tgDOS)/tcDOS > tolerance) {
        status = false ;
        ctest << "*";
      }
      ctest << endl ;
    }

    return status ;

  }

  bool UnitTests::parseInput(System* pSys) const {

    // Parse molecule data. 

    MoleculeManager* pMoleculeManager = pSys->getMoleculeManager() ;

    PersistPtr ppMolList = pMoleculeManager->get_PersistPtr();
    if(!ppMolList)	{
      cerr << "No molecules have been specified." << endl;
      return false;
    }

    PersistPtr ppmol = ppMolList ;
    while (ppmol = ppmol->XmlMoveTo("molecule")) { 

      Molecule *pMol = NULL ;
      // Get the name of the molcule.
      const char* reftxt = ppmol->XmlReadValue("id");
      if (reftxt) {
        pMol = pMoleculeManager->addmol(string(reftxt), string(""), pSys->getEnv(), pSys->m_Flags);
      }
    }

    return true ;
  }


  void UnitTests::underlineText(const string& text) const {

    ctest << "  " << text << endl ;
    ctest << "  " ;
    for (size_t i(0) ; i < text.size() ; i++ ) 
      ctest << "-" ;
    ctest << endl ;

  }

}//namespace

