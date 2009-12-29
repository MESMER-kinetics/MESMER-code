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

#include "../calcmethod.h"
#include "../dMatrix.h"

namespace mesmer
{
  class Fitting : public CalcMethod
  {
  public:

    Fitting(const std::string& id) : CalcMethod(id), m_nVar(0), m_A(), m_B(), m_C(), m_chi2a(0.0), m_chi2b(0.0), m_chi2c(0.0) {}

    virtual ~Fitting() {}

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

    // Get the current location.
    void GetLocation(vector<double> &loc) ;

    // Set the current location.
    void SetLocation(vector<double> &loc) const ;

    // Calculate the weighted sum of two vectors.
    vector<double> VectorAdd(const double a, const vector<double> &A, const double b, const vector<double> &B) const ;

    // Calculate the cartesian length of a vector.
    double VectorLength(const vector<double> &A) const ;

    // Normalize a vector.
    void VectorNormalize(vector<double> &A) const ;

    // Check that the a point falls within the limits defined by the user.
    bool CheckBounds(const vector<double> &A) const ;

    // Write out current variable values.
    void WriteVarVals() const ;

    // Check for line search convergence.
    bool CheckLineSearchConvergence(const vector<double> &X) const ;

    // Initialize the direcion vectors.
    void Fitting::initializeDirections(dMatrix &A) const ;

    // Update direction matrix in accord with the Powell algorithm.
    void Fitting::cycleDirections(dMatrix &A, const vector<double> &X) const ;

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

      // Return error.

      return false ;

    } else {

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

      pSys->calculate(chiSquare) ;

      double oldChiSquare = chiSquare ;

      WriteVarVals() ;
      
      // Set initial tolerance to be 1%
      
      double tol(0.01) ;

      //
      // Next, loop sequentially over each variable, optimizing it in each direction.
      // SHR, 29/Nov/2009: tests with the isopropyl system show that this is not optimal
      // and the something like to Powell conjugate direction method would improve the 
      // rate of convergence.
      //

      //  // converged = ( (oldChiSquare - chiSquare)/oldChiSquare ) < tol ;


      // Setup initial search directions.

      dMatrix directions(m_nVar,0.0); 

      initializeDirections(directions) ;                   

      vector<double> direction(m_nVar,0.0) ;

      for (size_t itr(1), count(0) ; itr <= 10 ; itr++) {

        // Perform an initial sweep across all vectors.

        for (size_t isweep(0); isweep < m_nVar ; isweep++) {

          cout << "Step " << ++count << " of fitting. chiSquare = " << chiSquare << endl;

          // Determine direction of search.

          for (size_t i(0); i < m_nVar ; i++ ) {
            direction[i] = directions[isweep][i] ;
          }

          oldChiSquare = chiSquare ;
          LineSearch(pSys, direction, chiSquare, tol);

          WriteVarVals() ;
        }

        // Calculate new search direction.

        vector<double> currentLocation(m_nVar,0.0) ;

        GetLocation(currentLocation) ;

        direction = VectorAdd(1.0, currentLocation, -1.0, initialLocation) ;
        VectorNormalize(direction) ;

        oldChiSquare = chiSquare ;
        LineSearch(pSys, direction, chiSquare, tol);

        WriteVarVals() ;

        // Update direction vectors in accord with the modified Powell algorithm.

        if ((itr % m_nVar) == 0 ) { 
          cout << endl << "Direction Matrix Reset." << endl;

          initializeDirections(directions) ;
          
          // tol = max(m_tol, tol/10.) ;                   
        } else {
          cycleDirections(directions,direction);
        }

      }  

    }

    return true;
  }

  void Fitting::LineSearch(System* pSys, const vector<double> &direction, double &currentChi2, double tol) {

    cout << endl << "Begin line search" << endl ;

    m_chi2a = currentChi2 ;

    // First bracket minimum.

    if (!BracketMinimum(pSys, direction)) {
      // failed to bracket minimum within user defined limits. 
      // Simply return for now.
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

    cout << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) << formatFloat(m_chi2c, 6, 15) << endl ;

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
        // Should throw, but will simply return for now.
        return false ;
      }

      SetLocation(m_C) ;
      pSys->calculate(m_chi2c);

      cout << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) << formatFloat(m_chi2c, 6, 15) << endl ;

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

    cout << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) 
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

      cout << formatFloat(m_chi2a, 6, 15) << formatFloat(m_chi2b, 6, 15) 
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
  void Fitting::WriteVarVals() const {

    cout << endl ;
    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar ; iVar++) {

      double var = *Rdouble::withRange()[iVar] ;
      formatFloat(cout, var, 6, 15) ;

    }
    cout << endl ;

  }

  //
  // Get the current location.
  //
  void Fitting::GetLocation(vector<double> &loc) {

    if (loc.size() != m_nVar) {
      // Throw an error.
    }

    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar ; iVar++) {
      loc[iVar] = *Rdouble::withRange()[iVar] ;
    }

  }

  //
  // Set the current location.
  //
  void Fitting::SetLocation(vector<double> &loc) const {

    if (loc.size() != m_nVar) {
      // Throw an error.
    }

    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar ; iVar++) {
      *Rdouble::withRange()[iVar] = loc[iVar] ;
    }

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
  // Check that the a point falls within the limits defined by the user.
  //
  bool Fitting::CheckBounds(const vector<double> &A) const {

    bool check(true) ;

    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar && check ; iVar++) {

      double var = A[iVar] ;
      double lower(0.0), upper(0.0), stepsize(0.0) ;

      Rdouble::withRange()[iVar]->get_range(lower, upper, stepsize) ;

      check = ((var > lower) && (var < upper)) ;

    }

    return check ;

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
      A[i][i] = 1.0 ;
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

