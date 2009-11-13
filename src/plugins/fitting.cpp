// fitting.cpp
// This class implements the Powell Convergent direction
// method to determine the minimum in the chi-squared surface.

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
  void LineSearch(System* pSys, const int varID);
};

////////////////////////////////////////////////
//Global instance
Fitting theFitting("fitting");
///////////////////////////////////////////////

bool Fitting::DoCalculation(System* pSys)
{
 //
  //
    size_t nVar = Rdouble::withRange().size() ;

    if (nVar < 1) { 

      // Return error.

    } else if (nVar == 1) {

      // Do a simple line search.
      LineSearch(pSys, 0);

    } else {

      int MaxNumSteps = 100 ;

      bool converged(false) ;

      double chiSquare(0.0);

      pSys->calculate(chiSquare);

      vector<double> currentPosition(nVar,0.0) ;

      dMatrix searchVectors(nVar);  // Create matrix of search directions.

      for (size_t i(0) ; i < nVar ; i++) { 
        searchVectors[i][i] = 1.0 ;
        currentPosition[i] = *Rdouble::withRange()[i];
      }

      for (int step(0); step < MaxNumSteps && !converged ; step++ ){

        for (size_t obj(0); obj < nVar ; ++obj){

        }

        ctest << "Step " << step << " of fitting. chiSquare = " << chiSquare << endl;
      }

    }
    return true;
  }

  void Fitting::LineSearch(System* pSys, const int varID) {

    static const double Gold = (3.0 - sqrt(5.0))/2.0 ;
    static const double GRatio = (1.0 - Gold)/Gold ;
    static const double tol = 1.0e-8 ;                 // This number should be an estimate of the square root of machine precision.
    static const int limit = 10 ;
    
    ofstream cfit("Fitting.res") ;

    cfit << endl << "Begin line search" << endl ;

    // First catch your hare ... need to bracket the minimum. To do this
    // use the parameter limits supplied. 

    double a, b, stepsize ;
    Rdouble::withRange()[varID]->get_range(a, b, stepsize);
    
    double diff = (a - b)/100.0  ;
    double mean = (a + b)/2.0 ;
    a = mean + diff ;
    b = mean - diff ;

    // Calculate chi2 at the upper and lower points.
    *Rdouble::withRange()[varID]=a ;
    double chi2a ;
    pSys->calculate(chi2a);

    *Rdouble::withRange()[varID]=b ;
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

    formatFloat(cfit, a,     6, 15) ;
    formatFloat(cfit, b,     6, 15) ;
    formatFloat(cfit, c,     6, 15) ;
    formatFloat(cfit, chi2a, 6, 15) ;
    formatFloat(cfit, chi2b, 6, 15) ;
    formatFloat(cfit, chi2c, 6, 15) ;
    cfit << endl ;

    // Repeat the search until a minimum has been bracketed or
    // the search limit has been reached. 
    
    int count(0) ;
    while (count < limit && chi2c < chi2b) {
      count++ ;

      // Shift values so as to maintain bracketing.
      a     = b ;
      chi2a = chi2b ;
      b     = c ;
      chi2b = chi2c ;

      // Determine next estimate of lower bracket point.

      c = b + GRatio*(b - a) ;
      *Rdouble::withRange()[varID]=c ;
      pSys->calculate(chi2c);

      formatFloat(cfit, a,     6, 15) ;
      formatFloat(cfit, b,     6, 15) ;
      formatFloat(cfit, c,     6, 15) ;
      formatFloat(cfit, chi2a, 6, 15) ;
      formatFloat(cfit, chi2b, 6, 15) ;
      formatFloat(cfit, chi2c, 6, 15) ;
      cfit << endl ;

    }

    // At this point the minimum should be bracketed, so 
    // use golden section search to refine the minimum.

    double x = c - Gold*(c - a) ;
    *Rdouble::withRange()[varID] = x ;
    double chi2x ;
    pSys->calculate(chi2x);

    count = 0 ;
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
      
      // Calculate chi2 for the new estimate of the minimum.
      
      formatFloat(cfit, a,     6, 15) ;
      formatFloat(cfit, b,     6, 15) ;
      formatFloat(cfit, c,     6, 15) ;
      formatFloat(cfit, chi2a, 6, 15) ;
      formatFloat(cfit, chi2b, 6, 15) ;
      formatFloat(cfit, chi2c, 6, 15) ;
      cfit << endl ;

    }
    
  }


}//namespace

