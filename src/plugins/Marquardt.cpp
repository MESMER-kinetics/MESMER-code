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
#include "FittingUtils.h"

namespace mesmer
{
  class Marquardt : public CalcMethod, private FittingUtils
  {
  public:

    Marquardt(const std::string& id) : CalcMethod(id), FittingUtils(), m_nVar(0), m_delta(0.001), m_lambdaScale(10.0) {}

    virtual ~Marquardt() {}

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:

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

    if (m_nVar < 1) { 
      cerr << "Marquardt requires at least one range variable to be set." << endl;
      return false ;
    }

    // Read in Marquardt parameters, or use values from defaults.xml.
    PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");
    m_delta = ppControl->XmlReadDouble("me:MarquardtDerivDelta");
    unsigned maxIterations= ppControl->XmlReadInteger("me:MarquardtIterations");
    double tol = ppControl->XmlReadDouble("me:MarquardtTolerance");

    // Read in parameter constraints.
//    ReadParameterConstraints(ppControl) ;

    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are caluclated).
    pSys->m_Flags.useTheSameCellNumber = true;

    // Warnings and less not sent to console.
    ChangeErrorLevel e(obError); 

    //Default is to disable ctest during Marquardt. Restored when leaving this function.
    //StopCTestOutput stop(!ppControl->XmlReadBoolean("me:ctestOutputWhenMarquardt")) ;

    //
    // Begin by finding the starting point chi-squared value.
    //

    vector<double> currentLocation(m_nVar,0.0) ; 
    vector<double> newLocation(m_nVar,0.0) ; 

    GetLocation(currentLocation) ;

    // Invoke SetLocation to catch any constrained parameters.
    SetLocation(currentLocation) ;

    double chiSquare(0.0), lambda(1.0) ;
    vector<double> residuals ;
    pSys->calculate(chiSquare, residuals) ;

    double bestChiSquare = chiSquare ;

    WriteVarVals(chiSquare, lambda) ;

    //
    // The following is slightly modified implementation of the Marquardt
    // algorithm. The modification is tha apllication of a bounds check on
    // each proposed new location of the minimum. If the the bounds check
    // fails lambda is increased with the consequence that the algorithm 
    // moves toward a short steepest decent algorithm. 
    //

    vector<double> gradient(m_nVar,0.0) ;
    dMatrix hessian(m_nVar,0.0); 
    NumericalDerivatives(pSys, residuals, m_delta, gradient, hessian) ;

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
          converged = (relativeChange < tol) ;
          lambda /= m_lambdaScale ;
          GetLocation(currentLocation) ;
          bestChiSquare = chiSquare ;
          NumericalDerivatives(pSys, residuals, m_delta, gradient, hessian) ;
        }
      } else {
        lambda *= m_lambdaScale ;
      }

      WriteVarVals(bestChiSquare, lambda) ;

      cinfo << "Iteration: " << itr << " of Marquardt. ChiSquare = " << bestChiSquare << ", Lambda = " << lambda << endl;

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
    ResultsAndStatistics(pSys, hessian) ;

    return true;
  }

  //
  // Write out current variable values.
  //
  void Marquardt::WriteVarVals(double chiSquare, double lambda) const {

    cerr << endl << "Chi^2 = " << chiSquare << " Lambda = " << lambda << endl ;
    for(size_t iVar(0) ; iVar < m_nVar ; iVar++) {

      Rdouble var = *Rdouble::withRange()[iVar] ;
      cerr << var.get_varname() << "=" << setprecision(6) << double(var) << "  "; 

    }
    cerr << endl ;

//    WriteConstrainedVals() ;

  }

} //namespace

