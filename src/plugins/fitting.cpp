//-------------------------------------------------------------------------------------------
//
// fitting.cpp
//
// Author: Struan Robertson
// Date:   29/Nov/2009
//
// This class implements the methods to determine the minimum in the chi-squared surface.
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
  class Fitting : public CalcMethod, private FittingUtils
  {
  public:

    Fitting(const char* id) : FittingUtils(), m_id(id), 
      m_nVar(0), m_A(), m_B(), m_C(), m_chi2a(0.0), m_chi2b(0.0), m_chi2c(0.0)
    { Register(); }

    virtual ~Fitting() {}
    virtual const char* getID()  { return m_id; }

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  protected:

  private:

    // Perform a golden search for a specified variable.
    void LineSearch(System* pSys, const vector<double> &direction, double &currentChi2, double tol) ;

    // Bracket minimum along a given line direction.    
    bool BracketMinimum(System* pSys, const vector<double> &direction) ;

    // Golden section search on a bracketed minimum.
    void GoldenSectionSearch(System* pSys, double &currentChi2, double tol) ;

    // Calculate the weighted sum of two vectors.
    vector<double> VectorAdd(const double a, const vector<double> &A, const double b, const vector<double> &B) const ;

    // Calculate the cartesian length of a vector.
    double VectorLength(const vector<double> &A) const ;

    // Normalize a vector.
    void VectorNormalize(vector<double> &A) const ;

    // Write out current variable values.
    void WriteVarVals(const double chiSquare) const ;

    // Check for line search convergence.
    bool CheckLineSearchConvergence(const vector<double> &X) const ;

    // Initialize the direcion vectors.
    void initializeDirections(dMatrix &A) const ;

    // Update direction matrix in accord with the Powell algorithm.
    void cycleDirections(dMatrix &A, const vector<double> &X) const ;

    const char* m_id;

    // Constants used in bracketing and Golden section search in a line minimizaton.
    // Note that the number tol should be an estimate of the square root of machine precision.
    static const double m_Gold ;
    static const double m_GRatio ;
    static const double m_tol ;

    // Dimension of fit.
    size_t m_nVar ;

    // Vectors to hold the position of various points during bracketing and Golden section search.
    vector<double> m_A, m_B, m_C ;

    // Values of the chi2 surface corresponding to the locations above.
    double m_chi2a, m_chi2b, m_chi2c ;
  };

  ////////////////////////////////////////////////
  //Global instance
  Fitting theFitting("fitting");
  ///////////////////////////////////////////////

  //
  // Standard constants used in line searches.
  //
  const double Fitting::m_Gold   = (3.0 - sqrt(5.0))/2.0 ;
  const double Fitting::m_GRatio = (1.0 - m_Gold)/m_Gold ;
  const double Fitting::m_tol    = 1.0e-8 ;

  bool Fitting::DoCalculation(System* pSys)
  {
    m_nVar = Rdouble::withRange().size() ;

    if (m_nVar < 1) { 
      cerr << "Fitting requires at least one range variable to be set." << endl;
      return false ;
    }

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

    //Default is to disable ctest during fitting. Restored when leaving this function.
    //StopCTestOutput stop(!ppControl->XmlReadBoolean("me:ctestOutputWhenFitting")) ;

    //
    // Initialize position vectors.
    //
    m_A.resize(m_nVar,0.0) ;
    m_B.resize(m_nVar,0.0) ;
    m_C.resize(m_nVar,0.0) ;

    //
    // Begin by finding the starting point chi-squared value.
    //
    double chiSquare(0.0) ;

    vector<double> initialLocation(m_nVar,0.0) ; 

    GetLocation(initialLocation) ;

	// Invoke SetLocation to catch any constrained parameters.
    SetLocation(initialLocation) ;

    pSys->calculate(chiSquare) ;

    double oldChiSquare = chiSquare ;

    WriteVarVals(chiSquare) ;

	//
    // The following implementation is loosely based on the Powell method. An initial
    // sweep is performed over all directions and from these a new direction for 
    // search is calculated and is substituted for the last vector in the direction 
    // matrix. The direaction matrix is periodically re-initialized to prevent linear
    // dependence problems.
    //

    // Setup initial search directions.

    dMatrix directions(m_nVar,0.0); 

    initializeDirections(directions) ;                   

    vector<double> direction(m_nVar,0.0) ;

    for (size_t itr(1) ; itr <= maxIterations ; itr++) {

      cinfo << "Iteration: " << itr << " of fitting. chiSquare = " << chiSquare << endl;

      // Perform an initial sweep across all vectors.

      for (size_t isweep(0); isweep < m_nVar ; isweep++) {

        cinfo << "Direction sweep:" << isweep << " of fitting. chiSquare = " << chiSquare << endl;

        // Determine direction of search.

        for (size_t i(0); i < m_nVar ; i++ ) {
          direction[i] = directions[i][isweep] ;
        }

        oldChiSquare = chiSquare ;
        LineSearch(pSys, direction, chiSquare, tol);

        WriteVarVals(chiSquare) ;
      }

      // Calculate new search direction.

      vector<double> currentLocation(m_nVar,0.0) ;

      GetLocation(currentLocation) ;

      direction = VectorAdd(1.0, currentLocation, -1.0, initialLocation) ;
      VectorNormalize(direction) ;

      oldChiSquare = chiSquare ;
      LineSearch(pSys, direction, chiSquare, tol);

      WriteVarVals(chiSquare) ;

      // Update direction vectors in accord with the modified Powell algorithm.

      if ((itr % m_nVar) == 0 ) { 
        cinfo << endl << "Direction Matrix Reset." << endl;

        initializeDirections(directions) ;

        // tol = max(m_tol, tol/10.) ;                   
      } else {
        cycleDirections(directions,direction);
      }

    }

    // Write meta data to the XML file.
    for (size_t i(0); i < m_nVar ; i++ ) {
      TimeCount events;
      std::string timeString;
      Rdouble::withRange()[i]->XmlWriteAttribute("fitted", events.setTimeStamp(timeString));
      stringstream cs;
      cs << chiSquare;
      Rdouble::withRange()[i]->XmlWriteAttribute("chiSquared", cs.str());
    }

    // Write out the results and the statisitics of the fit.
    vector<double> residuals ;
    pSys->calculate(chiSquare, residuals) ;

    vector<double> gradient(m_nVar,0.0) ;
    dMatrix hessian(m_nVar,0.0); 
    double delta(0.001) ;
    NumericalDerivatives(pSys, residuals, delta, gradient, hessian) ;

    ResultsAndStatistics(pSys, hessian) ;

    return true;
  }

  void Fitting::LineSearch(System* pSys, const vector<double> &direction, double &currentChi2, double tol) {

    cinfo << endl << "Begin line search" << endl ;

    m_chi2a = currentChi2 ;

    // First bracket minimum.

    if (!BracketMinimum(pSys, direction)) {
      // failed to bracket minimum within user defined limits. 
      // Simply return for now.
      currentChi2 = m_chi2a ;
      return ;
    }

    // At this point the minimum should be bracketed, so 
    // use golden section search to refine the minimum.

    GoldenSectionSearch(pSys, currentChi2, tol) ;

  }

  //
  // First catch your hare ... need to bracket the minimum. To do this
  // use the parameter limits supplied. 
  //
  bool Fitting::BracketMinimum(System* pSys, const vector<double> &direction) {

    // Get the current best estimate of the location of the chi2 minimum.

    GetLocation(m_A) ;

    //
    // SHR 13/Dec/2009: Need some criteria, probably based on range, that can be
    // used to determine how far in the direction of search to proceed.
    //

    m_B = VectorAdd(1.0, m_A, 1.0, direction) ;

    // Check bounds.    
    if (!CheckBounds(m_B)) {
      // Should throw, but will simply return for now.
      return false ;
    }

    // Calculate chi2 for new point.

    SetLocation(m_B) ;
    pSys->calculate(m_chi2b);

    // Alter the direction of search so that we are always going down hill.

    if (m_chi2a < m_chi2b) {
      vector<double> vtmp = m_A ;
      m_A = m_B ;
      m_B = vtmp ;
      double tmp = m_chi2a ;
      m_chi2a = m_chi2b ;
      m_chi2b = tmp ;
    }

    // Follow gradient down hill to estimate location of the next point.

    // c = b + m_GRatio*(b - a) ;
    m_C = VectorAdd( (1.0 + m_GRatio), m_B, (-m_GRatio), m_A) ;

    // Check bounds.    
    if (!CheckBounds(m_C)) {
      // Should throw, but will simply return for now.
      return false ;
    }

    // Calculate a new value of chi2 for the new parameter value. 

    SetLocation(m_C) ;
    pSys->calculate(m_chi2c);

    cinfo << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) << formatFloat(m_chi2c, 6, 15) << endl ;

    // Repeat the search until a minimum has been bracketed or
    // the search limit has been reached. 

    while (m_chi2c < m_chi2b) {

      // Shift values so as to maintain bracketing.
      m_A     = m_B ;
      m_chi2a = m_chi2b ;
      m_B     = m_C ;
      m_chi2b = m_chi2c ;

      // Determine next estimate of lower bracket point.

      // c = b + m_GRatio*(b - a) ;
      m_C = VectorAdd( (1.0 + m_GRatio), m_B, (-m_GRatio), m_A) ;

      // Check bounds.    
      if (!CheckBounds(m_C)) {
        cinfo << endl << "Bound check failed in bracket search." << endl ;
        return false ;
      }

      SetLocation(m_C) ;
      pSys->calculate(m_chi2c);

      cinfo << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) << formatFloat(m_chi2c, 6, 15) << endl ;

    }

    return true ;

  }

  //
  // Golden section search on a bracketed minimum.
  //
  void Fitting::GoldenSectionSearch(System* pSys, double &currentChi2, double tol) {

    static const int limit = 10 ;

    // x = c - m_Gold*(c - a) ;
    vector<double> X = VectorAdd( (1.0 - m_Gold), m_C, m_Gold, m_A) ;
    SetLocation(X) ;
    double chi2x ;
    pSys->calculate(chi2x);

    cinfo << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) 
      << formatFloat(  chi2x, 6, 15) << formatFloat(m_chi2c, 6, 15) << endl ;

    int count = 0 ;
    bool converged(false) ;
    while(count < limit && !converged) {
      count++ ;

      if (chi2x < m_chi2b) {
        m_A     = m_B ;
        m_chi2a = m_chi2b ;
        m_B     = X ;
        m_chi2b = chi2x ;

        // x = c - m_Gold*(c - a) ;
        X = VectorAdd( (1.0 - m_Gold), m_C, m_Gold, m_A) ;
        SetLocation(X) ;

        pSys->calculate(chi2x);

        converged = (fabs((chi2x/m_chi2b) - 1.0) < tol) || CheckLineSearchConvergence(X) ;

      } else {
        m_C     = X ;
        m_chi2c = chi2x ;
        X       = m_B ;
        chi2x   = m_chi2b ;

        // b = a + m_Gold*(c - a) ;
        m_B = VectorAdd( (1.0 - m_Gold), m_A, m_Gold, m_C) ;
        SetLocation(m_B) ;

        pSys->calculate(m_chi2b);

        converged = (fabs((m_chi2b/chi2x) - 1.0) < tol) || CheckLineSearchConvergence(X) ;

      }

      cinfo << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) 
        << formatFloat(  chi2x, 6, 15) << formatFloat(m_chi2c, 6, 15) << endl ;

    }

    // Save the value with the best Chi^2 value.

    if (chi2x < m_chi2b) {
      currentChi2 = chi2x ;
      SetLocation(X) ;
    } else {
      currentChi2 = m_chi2b ;
      SetLocation(m_B) ;
    }

  }

  //
  // Write out current variable values.
  //
  void Fitting::WriteVarVals(double chiSquare) const {

    cerr << endl << "Chi^2 = " << chiSquare << endl ;

    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar ; iVar++) {

      Rdouble var = *Rdouble::withRange()[iVar] ;
      cerr << var.get_varname() << "=" << setprecision(6) << double(var) << "  "; 

    }
    cerr << endl ;

  }


  //
  // Calculate the weighted sum of two vectors.
  //
  vector<double> Fitting::VectorAdd(const double a, const vector<double> &A, const double b, const vector<double> &B) const {

    if (A.size() != B.size()) {
      // Throw an error.
    }

    vector<double> sum(A.size(),0.0) ;
    vector<double>::size_type iVar ;
    for(iVar = 0 ; iVar < A.size() ; iVar++) {
      sum[iVar] += a*A[iVar] + b*B[iVar] ;
    }

    return sum ;

  }

  //
  // Calculate the cartesian length of a vector.
  //
  double Fitting::VectorLength(const vector<double> &A) const {

    if (A.size() == 0) {
      // Throw an error.
    }

    vector<double>::size_type iVar ;
    double sum(0.0) ;
    for(iVar = 0 ; iVar < A.size() ; iVar++) {
      sum += A[iVar]*A[iVar] ;
    }

    return sqrt(sum) ;

  }

  //
  // Normalize a vector.
  //
  void Fitting::VectorNormalize(vector<double> &A) const {

    if (A.size() == 0) {
      // Throw an error.
    }

    vector<double>::size_type iVar ;
    double norm = VectorLength(A) ;
    if ( norm > 0.0) {
      for(iVar = 0 ; iVar < A.size() ; iVar++) {
        A[iVar] /= norm ;
      }
    } else {
      // Throw an error.
    }

  }

  //
  // Check for line search convergence.
  // fabs(|c-a|) > m_tol*(fabs(b)+fabs(x)) )
  //
  bool Fitting::CheckLineSearchConvergence(const vector<double> &X) const {

    // bool converged(false) ;

    vector<double> vtmp = VectorAdd(1.0, m_C, -1.0, m_A) ;
    double interval = VectorLength(vtmp) ;

    vtmp = VectorAdd( 1.0, m_B, 1.0, X) ;
    double radius = VectorLength(vtmp) ;

    // converged = (interval < m_tol*radius) ;

    // return !converged ;

    return !(interval > m_tol*radius) ;
  }

  //
  // Initialize the direcion vectors.
  //
  void Fitting::initializeDirections(dMatrix &A) const {

    if (A.size() == 0) {
      // Throw an error.
    }

    size_t i, j ;
    for(i = 0 ; i < A.size() ; i++) {
      for(j = 0 ; j < A.size() ; j++) {
        A[j][i] = 0.0 ;
      }
      double lower(0.0), upper(0.0), stepsize(0.0) ;
      Rdouble::withRange()[i]->get_range(lower, upper, stepsize) ;
      A[i][i] = (stepsize > 0.0) ? stepsize : 1.0 ;
    }

  }

  //
  // Update direction matrix in accord with the Powell algorithm.
  //
  void Fitting::cycleDirections(dMatrix &A, const vector<double> &X) const {

    if (A.size() == 0) {
      // Throw an error.
    }

    size_t i, j, isize(A.size()-1) ;
    for(i = 0 ; i < isize ; i++) {
      for(j = 0 ; j < A.size() ; j++) {
        A[j][i] = A[j][i+1] ;
      }
    }

    for(j = 0 ; j < A.size() ; j++) {
      A[j][isize] = X[j] ;
    }

  }


}//namespace

/*
Writing the fitted result to the XML file
The fitted variable came originally from a structure like:
<property dictRef="me:ZPE">
<scalar units="cm-1" lower="10000" upper="14000" stepsize="10">
13000.0
</scalar>
</property>

The result should be like:
<property dictRef="me:ZPE">
<scalar units="cm-1" lower="10000" upper="14000" stepsize="10" fitted="20080705_104810">
11234.5
</scalar>
</property>

The result file could be used to:
repeat the fitting run,        leaving <me:calcMethod> at Fitting;
do a normal range calculation, after changing the <me:calcMethod> to GridSearch;
use only the fitted value,     after changing <me:calcMethod> to SimpleCalc.
This currently happens like this.

A PersistPtr to <scalar> (or equivalent element) needs to be stored. Then 
stringsteam ss;
ss << fittedval;
pp->XmlWrite(ss.str()); //needs to be written
pp->XmlWriteAttribute("fitted",events.setTimeStamp(timeString)); 
*/
