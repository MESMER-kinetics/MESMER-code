//-------------------------------------------------------------------------------------------
//
// FittingUtils.cpp
//
// Author: Struan Robertson
// Date:   16/Jul/2011
//
// Implementation of a utility class that is inherited by both Powell and Marquardt methods. 
//
//-------------------------------------------------------------------------------------------

#include "../System.h"
#include "FittingUtils.h"

namespace mesmer
{
  //
  // Get the current location.
  //
  void FittingUtils::GetLocation(vector<double> &loc) const {

    for(size_t iVar(0) ; iVar < loc.size() ; iVar++) {
      loc[iVar] = *Rdouble::withRange()[iVar] ;
    }

  }

  //
  // Set the current location.
  //
  void FittingUtils::SetLocation(vector<double> &loc) const {

    for(size_t iVar(0) ; iVar < loc.size() ; iVar++) {
      *Rdouble::withRange()[iVar] = loc[iVar] ;
      Rdouble::withRange()[iVar]->XmlWriteValue() ;
    }

  }

  //
  // Check that the a point falls within the limits defined by the user.
  //
  bool FittingUtils::CheckBounds(const vector<double> &A) const {

    bool check(true) ;
    for(size_t iVar(0) ; iVar < A.size() && check ; iVar++) {

      double var = A[iVar] ;
      double lower(0.0), upper(0.0), stepsize(0.0) ;

      Rdouble::withRange()[iVar]->get_range(lower, upper, stepsize) ;

      check = ((var > lower) && (var < upper)) ;

    }

    return check ;

  }

  //
  // The following methods calculates the numnerical derivatives and hessian of the 
  // Chi^2 surface. More specifically it calculates the -1/2 of the derivative and 
  // the 1/2 of the hessian as it these quantities can be used directly in the 
  // Marquardt algorithm to determine the next best guess of the the minimum.
  //
  void FittingUtils::NumericalDerivatives(System* pSys, vector<double> &residuals, double delta, vector<double> &gradient, qdMatrix &hessian) const {

    size_t nVar = gradient.size() ;
    vector<double> location(nVar,0.0), update(nVar,0.0) ;
    vector<double> derivatives ;
    GetLocation(location) ;
    for (size_t iVar(0) ; iVar < nVar ; iVar++) {

      update = location ;
      update[iVar] *= (1.0 + delta) ;
      SetLocation(update) ;

      double chiSquare(0.0) ;
      vector<double> newResiduals ;
      pSys->calculate(chiSquare, newResiduals) ;

      size_t sizeRes = residuals.size() ;
      if (newResiduals.size() != sizeRes) {
        cerr << "Error: residual vectors are of different size" ;
      }

      double grad(0.0), hess(0.0) ; 
      for (size_t i(0) ; i < sizeRes ; i++) {
        double deriv = (residuals[i] - newResiduals[i])/(delta*location[iVar]) ;
        grad += residuals[i]*deriv ;
        hess += deriv*deriv ;
        derivatives.push_back(deriv) ;
      }
      gradient[iVar] = grad ;
      hessian[iVar][iVar] = hess ;

      for (size_t jVar(0) ; jVar < iVar ; jVar++) {
        hess = 0.0 ;
        for (size_t i(0),ii(iVar*sizeRes),jj(jVar*sizeRes) ; i < sizeRes ; i++, ii++, jj++) {
          hess += derivatives[ii]*derivatives[jj] ;
        }
        hessian[iVar][jVar] = hessian[jVar][iVar] = qd_real(hess) ;
      }

    }

    SetLocation(location) ;

  }

  //
  // Write out the results and statistics of the fit. 
  //
  void FittingUtils::ResultsAndStatistics(System* pSys, qdMatrix &hessian) const {

    // Calculate model values with optimum parameters.

    double chiSquare(0.0) ;
    vector<double> residuals ;
    pSys->calculate(chiSquare, residuals, true) ;

    // Calculate covaraince matrix.

    hessian.invertLUdecomposition() ;

    cinfo << endl << "Chi^2 = " << chiSquare << endl << endl << "Best fit parameters:" << endl << endl ;

    // Best fit parameters.

    for(size_t iVar(0) ; iVar < hessian.size() ; iVar++) {

      Rdouble var = *Rdouble::withRange()[iVar] ;
      double sigma = to_double(sqrt(hessian[iVar][iVar])) ;
      cinfo << var.get_varname() << " = " << setprecision(6) << var.originalUnits() << " +/- " << var.originalUnits(sigma) << endl; 

    }

    // Correlation coefficients.

    cinfo << endl << "Correlation coefficients:" << endl << endl ;

    for(size_t iVar(0) ; iVar < hessian.size() ; iVar++) {

      Rdouble vara = *Rdouble::withRange()[iVar] ;
      double sigma = to_double(sqrt(hessian[iVar][iVar])) ;
      for(size_t jVar(0) ; jVar < iVar ; jVar++) {
        double corrlCoeff = to_double(hessian[iVar][jVar]/(sigma*sqrt(hessian[jVar][jVar]))) ;
        Rdouble varb = *Rdouble::withRange()[jVar] ;
        cinfo << vara.get_varname() << " , " << varb.get_varname() << " = " << setprecision(6) << corrlCoeff << endl; 
      }

    }

    // Goodness of fit.

    size_t NoDegFreedom = residuals.size() - hessian.size() ;
    cinfo << endl << "Goodness of Fit:" << endl << endl ;
    cinfo << "Number of degrees of Freedom = " <<  NoDegFreedom << endl ;
    cinfo << "Chi^2 probability = " << ChiSquaredPrbFn(chiSquare/2.0, double(NoDegFreedom)/2.0) << endl << endl ;

    return ;

  }

}//namespace
