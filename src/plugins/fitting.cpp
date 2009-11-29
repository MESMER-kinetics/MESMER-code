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

#include "../calcmethod.h"
#include <fstream>

namespace mesmer
{
  class Fitting : public CalcMethod
  {
  public:
    Fitting(const std::string& id) : CalcMethod(id) {}
    virtual ~Fitting() {}

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  protected:

  private:

    // Locate the starting point of the search as the middle of the allowed range. 
    double StartingPoint(System* pSys, size_t nVar) const ;

    // Perform a golden search for a specified variable.
    void LineSearch(System* pSys, const int varID, double &currentChi2) const ;
  };

  ////////////////////////////////////////////////
  //Global instance
  Fitting theFitting("fitting");
  ///////////////////////////////////////////////

  bool Fitting::DoCalculation(System* pSys)
  {
    size_t nVar = Rdouble::withRange().size() ;

    if (nVar < 1) { 

      // Return error.

    } else {

      //
      // Begin by finding the starting point of the search.
      //

      double chiSquare = StartingPoint(pSys, nVar) ;

      //
      // Next, loop sequentially over each variable, optimizing it in each direction.
      // SHR, 29/Nov/2009: tests with the isopropyl system show that this is not optimal
      // and the something like to Powell conjugate direction method would improve the 
      // rate of convergence.
      //
      int MaxNumSteps(5) ;
      bool converged(false) ;
      double tol(0.1) ;

      cout << endl ;
      size_t iVar ;
      for(iVar = 0 ; iVar < nVar ; iVar++) {

        double var = *Rdouble::withRange()[iVar] ;
        formatFloat(cout, var, 6, 15) ;

      }
      cout << endl ;

      for (int step(1); step <= MaxNumSteps && !converged ; step++ ){

        cout << "Step " << step << " of fitting. chiSquare = " << chiSquare << endl;

        int varID = step % nVar ;

        double oldChiSquare = chiSquare ;
        LineSearch(pSys, varID, chiSquare);

        cout << endl ;
        for(iVar = 0 ; iVar < nVar ; iVar++) {

          double var = *Rdouble::withRange()[iVar] ;
          formatFloat(cout, var, 6, 15) ;

        }
        cout << endl ;

        // converged = ( (oldChiSquare - chiSquare)/oldChiSquare ) < tol ;

      }

    }
    return true;
  }

  //
  // Locate the starting point of the search as the middle of the allowed range. 
  //
  double Fitting::StartingPoint(System* pSys, size_t nVar) const {

    // Locate starting point of search as the middle of the range.

    for(size_t varID(0) ; varID < nVar ; varID++) {

      double a, b, stepsize ;

      Rdouble::withRange()[varID]->get_range(a, b, stepsize);
      *Rdouble::withRange()[varID] = (a + b)/2.0 ;

    }

    double chi2 ;
    pSys->calculate(chi2) ;

    return chi2 ;

  }

  void Fitting::LineSearch(System* pSys, const int varID, double &currentChi2) const {

    static const double Gold = (3.0 - sqrt(5.0))/2.0 ;
    static const double GRatio = (1.0 - Gold)/Gold ;
    static const double tol = 1.0e-8 ;                 // This number should be an estimate of the square root of machine precision.
    static const int limit = 10 ;

    cout << endl << "Begin line search" << endl ;

    double chi2a = currentChi2 ;

    // First catch your hare ... need to bracket the minimum. To do this
    // use the parameter limits supplied. 

    double a, b, stepsize ;
    Rdouble::withRange()[varID]->get_range(a, b, stepsize);

    double diff = (a - b)/100.0 ;
    a = *Rdouble::withRange()[varID] ;
    b = a - diff ;

    // Calculate for new point.

    *Rdouble::withRange()[varID] = b ;
    double chi2b ;
    pSys->calculate(chi2b);

    // Alter the direction of search so that we are always going down hill.

    if (chi2a < chi2b) {
      double tmp = a ;
      a = b ;
      b = tmp ;
      tmp = chi2a ;
      chi2a = chi2b ;
      chi2b = tmp ;
    }

    // Follow gradient down hill to estimate location of the next point.

    double c = b + GRatio*(b - a) ;

    // Calculate a new value of chi2 for the new parameter value. 

    *Rdouble::withRange()[varID] = c ;
    double chi2c ;
    pSys->calculate(chi2c);

    formatFloat(cout, chi2a, 6, 15) ;
    formatFloat(cout, chi2b, 6, 15) ;
    formatFloat(cout, chi2c, 6, 15) ;
    cout << endl ;

    // Repeat the search until a minimum has been bracketed or
    // the search limit has been reached. 

    while (chi2c < chi2b) {

      // Shift values so as to maintain bracketing.
      a     = b ;
      chi2a = chi2b ;
      b     = c ;
      chi2b = chi2c ;

      // Determine next estimate of lower bracket point.

      c = b + GRatio*(b - a) ;
      *Rdouble::withRange()[varID]=c ;
      pSys->calculate(chi2c);

      formatFloat(cout, chi2a, 6, 15) ;
      formatFloat(cout, chi2b, 6, 15) ;
      formatFloat(cout, chi2c, 6, 15) ;
      cout << endl ;

    }

    // At this point the minimum should be bracketed, so 
    // use golden section search to refine the minimum.

    double x = c - Gold*(c - a) ;
    *Rdouble::withRange()[varID] = x ;
    double chi2x ;
    pSys->calculate(chi2x);

    formatFloat(cout, chi2a, 6, 15) ;
    formatFloat(cout, chi2b, 6, 15) ;
    formatFloat(cout, chi2x, 6, 15) ;
    formatFloat(cout, chi2c, 6, 15) ;
    cout << endl ;

    int count = 0 ;
    while(count < limit && fabs(c-a) > tol*(fabs(b)+fabs(x)) ){
      count++ ;

      if (chi2x < chi2b) {
        a = b ;
        b = x ;
        x = c - Gold*(c - a) ;
        chi2a = chi2b ;
        chi2b = chi2x ;
        *Rdouble::withRange()[varID] = x ;
        pSys->calculate(chi2x);

      } else {
        c = x ;
        x = b ;
        b = a + Gold*(c - a) ;
        chi2c = chi2x ;
        chi2x = chi2b ;
        *Rdouble::withRange()[varID] = b ;
        pSys->calculate(chi2b);

      }

      formatFloat(cout, chi2a, 6, 15) ;
      formatFloat(cout, chi2b, 6, 15) ;
      formatFloat(cout, chi2x, 6, 15) ;
      formatFloat(cout, chi2c, 6, 15) ;
      cout << endl ;

    }

    // Save the value with the best Chi^2 value.
    
    if (chi2x < chi2b) {
      currentChi2 = chi2x ;
      *Rdouble::withRange()[varID] = x ;
    } else {
      currentChi2 = chi2b ;
      *Rdouble::withRange()[varID] = b ;
    }
    
  }


}//namespace

