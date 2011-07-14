//-------------------------------------------------------------------------------------------
//
// Marquardt.cpp
//
// Author: Struan Robertson
// Date:   19/Jun/2011
//
// This class implements the Levenberg-Marquardt non-linear least squares algorithm. 
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"
#include "../dMatrix.h"

namespace mesmer
{
  class Marquardt : public CalcMethod
  {
  public:

    Marquardt(const std::string& id) : CalcMethod(id), m_nVar(0), m_delta(0.001), m_lambdaScale(10.0) {}

    virtual ~Marquardt() {}

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  protected:

  private:

    // Numerical derivatives.
    void NumericalDerivatives(System* pSys, vector<double> &residuals, vector<double> &gradient, dMatrix &hessian) const ;

    // Get the current location.
    void GetLocation(vector<double> &loc) const ;

    // Set the current location.
    void SetLocation(vector<double> &loc) const ;

    // Check that the a point falls within the limits defined by the user.
    bool CheckBounds(const vector<double> &A) const ;

    // Write out current variable values.
    void WriteVarVals(double chiSquare, double lambda) const ;

    // Dimension of fit.
    size_t m_nVar ;

    // Numerical derivative delta.
    double m_delta ;

    // Factor for scalling lambda by during fitting.
    double m_lambdaScale ;

  };

  ////////////////////////////////////////////////
  //Global instance
  Marquardt theMarquardt("marquardt");
  ///////////////////////////////////////////////

  bool Marquardt::DoCalculation(System* pSys)
  {
    m_nVar = Rdouble::withRange().size() ;
    assert(m_nVar == RangeXmlPtrs.size());

    if (m_nVar < 1) { 
      cerr << "Marquardt requires at least one range variable to be set." << endl;
      return false ;
    }

    //Read in Marquardt parameters, or use values from defaults.xml.
    PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");
    m_delta = ppControl->XmlReadDouble("me:MarquardtDerivDelta");
    unsigned maxIterations= ppControl->XmlReadInteger("me:MarquardtIterations");
    double tol = ppControl->XmlReadDouble("me:MarquardtTolerance");

    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are caluclated).
    pSys->m_Flags.useTheSameCellNumber = true;

    //Default is to disable ctest during Marquardt. Restored when leaving this function.
    //StopCTestOutput stop(!ppControl->XmlReadBoolean("me:ctestOutputWhenMarquardt")) ;

    //
    // Begin by finding the starting point chi-squared value.
    //

    vector<double> currentLocation(m_nVar,0.0) ; 
    vector<double> newLocation(m_nVar,0.0) ; 

    GetLocation(currentLocation) ;

    double chiSquare(0.0), lambda(1.0) ;
    vector<double> residuals ;
    pSys->calculate(chiSquare, residuals) ;

    double bestChiSquare = chiSquare ;

    WriteVarVals(chiSquare, lambda) ;

    ChangeErrorLevel e(obError); // Warnings and less not sent to console.

    //
    // The following is slightly modified implementation of the Marquardt
    // algorithm. The modification is tha apllication of a bounds check on
    // each proposed new location of the minimum. If the the bounds check
    // fails lambda is increased with the consequence that the algorithm 
    // moves toward a short steepest decent algorithm. 
    //

    vector<double> gradient(m_nVar,0.0) ;

    dMatrix hessian(m_nVar,0.0); 

    NumericalDerivatives(pSys, residuals, gradient, hessian) ;

    bool converged(false) ;
    for (size_t itr(1) ; itr <= maxIterations && !converged ; itr++) {

      newLocation = currentLocation ;
      vector<double> deltaLocation = gradient;
      dMatrix invHessian = hessian ;

      for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
        invHessian[iVar][iVar] *= (1.0 + lambda) ;
      }

      invHessian.invertGaussianJordan() ;

      deltaLocation *= invHessian ; 

      for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
        newLocation[iVar] += deltaLocation[iVar] ;
      }

      // Check bounds.    
      if (CheckBounds(newLocation)) {
        SetLocation(newLocation) ;

        pSys->calculate(chiSquare, residuals) ;

        if (chiSquare > bestChiSquare) {
          lambda *= m_lambdaScale ;
          SetLocation(currentLocation) ;
        } else {
          double relativeChange = 1.0 - chiSquare/bestChiSquare ;
          cerr << relativeChange << endl ;
		  converged = (relativeChange < 0.001)? true : false ;
		  lambda /= m_lambdaScale ;
		  GetLocation(currentLocation) ;
		  bestChiSquare = chiSquare ;
		  NumericalDerivatives(pSys, residuals, gradient, hessian) ;
        }
      } else {
        lambda *= m_lambdaScale ;
      }

      WriteVarVals(bestChiSquare, lambda) ;

      cinfo << "Iteration: " << itr << " of Marquardt. chiSquare = " << chiSquare << lambda << endl;

    }

    // Write the optimized result to the XML file.
    for (size_t i(0); i < m_nVar ; i++ ) {
      RangeXmlPtrs[i]->XmlWrite(toString(*Rdouble::withRange()[i]));

      TimeCount events;
      std::string timeString;
      RangeXmlPtrs[i]->XmlWriteAttribute("fitted", events.setTimeStamp(timeString));
      stringstream cs;
      cs << chiSquare;
      RangeXmlPtrs[i]->XmlWriteAttribute("chiSquared", cs.str());
    }

    // Calculate model values with optimum parameters.

    pSys->calculate(chiSquare, true) ;

    // Calculate covaraince matrix.

    hessian.invertGaussianJordan() ;

    cinfo << endl << "Chi^2 = " << chiSquare << endl << endl << "Best fit parameters:" << endl << endl ;

    // Best fit parameters.

    for(size_t iVar(0) ; iVar < m_nVar ; iVar++) {

      Rdouble var = *Rdouble::withRange()[iVar] ;
      double sigma = sqrt(hessian[iVar][iVar]) ;
      cinfo << var.get_varname() << " = " << setprecision(6) << double(var) << " +/- " << sigma << endl; 

    }

    // Correlation coefficients.

    cinfo << endl << "Correlation coefficients:" << endl << endl ;

    for(size_t iVar(0) ; iVar < m_nVar ; iVar++) {

      Rdouble vara = *Rdouble::withRange()[iVar] ;
      double sigma = sqrt(hessian[iVar][iVar]) ;
      for(size_t jVar(0) ; jVar < iVar ; jVar++) {
        double corrlCoeff = hessian[iVar][jVar]/(sigma*sqrt(hessian[jVar][jVar])) ;
        Rdouble varb = *Rdouble::withRange()[jVar] ;
        cinfo << vara.get_varname() << " , " << varb.get_varname() << " = " << setprecision(6) << corrlCoeff << endl; 
      }

    }

    // Goodness of fit.

    cinfo << endl << "Goodness of Fit:" << endl << endl ;
	cinfo << "Number of degrees of Freedom = " << residuals.size() - m_nVar << endl ;
	// cinfo << "Chi^2 probability = " << ChiSquaredPrbFn(chiSquare, double(residuals.size() - m_nVar)) << endl ;

    return true;
  }

  //
  // The following methods calculates the numnerical derivatives and hessian of the 
  // Chi^2 surface. More specifically it calculates the -1/2 of the derivative and 
  // the 1/2 of the hessian as it these quantities can be used directly in the 
  // Marquardt algorithm to determine the next best guess of the the minimum.
  //
  void Marquardt::NumericalDerivatives(System* pSys, vector<double> &residuals, vector<double> &gradient, dMatrix &hessian) const {

    vector<double> location(m_nVar,0.0), update(m_nVar,0.0) ;
    vector<double> derivatives ;
    GetLocation(location) ;
    for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {

      update = location ;
      update[iVar] *= (1.0 + m_delta) ;
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
        double deriv = (residuals[i] - newResiduals[i])/(m_delta*location[iVar]) ;
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
        hessian[iVar][jVar] = hessian[jVar][iVar] = hess ;
      }

    }

    SetLocation(location) ;

  }

  //
  // Write out current variable values.
  //
  void Marquardt::WriteVarVals(double chiSquare, double lambda) const {

    cerr << endl << "Chi^2 = " << chiSquare << " Lambda = " << lambda << endl ;

    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar ; iVar++) {

      Rdouble var = *Rdouble::withRange()[iVar] ;
      cerr << var.get_varname() << "=" << setprecision(6) << double(var) << "  "; 

    }
    cerr << endl ;

  }

  //
  // Get the current location.
  //
  void Marquardt::GetLocation(vector<double> &loc) const {

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
  void Marquardt::SetLocation(vector<double> &loc) const {

    if (loc.size() != m_nVar) {
      // Throw an error.
    }

    size_t iVar ;
    for(iVar = 0 ; iVar < m_nVar ; iVar++) {
      *Rdouble::withRange()[iVar] = loc[iVar] ;
    }

  }

  //
  // Check that the a point falls within the limits defined by the user.
  //
  bool Marquardt::CheckBounds(const vector<double> &A) const {

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

} //namespace

