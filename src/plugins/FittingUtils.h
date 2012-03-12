#ifndef GUARD_FittingUtils_h
#define GUARD_FittingUtils_h

//-------------------------------------------------------------------------------------------
//
// FittingUtils.h
//
// Author: Struan Robertson
// Date:   16/Jul/2011
//
// Definition of a utility class that is inherited by both Powell and Marquardt methods. 
//
//-------------------------------------------------------------------------------------------

namespace mesmer
{

  class FittingUtils
  {
  public:

	FittingUtils() : m_parameterConstraints() {} ;
	~FittingUtils() {} ; 

  protected:

	// Read parameter constraints.
	void ReadParameterConstraints(PersistPtr ppControl) ;

	// Get the current location.
	void GetLocation(vector<double> &loc) const ;

	// Set the current location.
	void SetLocation(vector<double> &loc) const ;

	// Check that the a point falls within the limits defined by the user.
	bool CheckBounds(const vector<double> &A) const ;

	// Numerical derivatives.
	void NumericalDerivatives(System* pSys, vector<double> &residuals, double delta, vector<double> &gradient, dMatrix &hessian) const ;

	// Write out the results and statistics of the fit. 
	void ResultsAndStatistics(System* pSys, dMatrix &hessian) const ;

  private :

	//
	// Private class that maintains a parameter constraint.
	//
	class ParameterConstraint {
	public:
	  ParameterConstraint(std::string indpndPmtr, std::string dpndPmtr, double factor, double addand ): m_indpndPmtr(indpndPmtr), 
		m_dpndPmtr(dpndPmtr), 
		m_factor(factor), 
		m_addand(addand) {} 
	  ~ParameterConstraint() {} 

	  // Calculate dependent parameter values.

	  void calcDpndPmtrValues() ;

	private:
	  ParameterConstraint() {}  
	  std::string m_indpndPmtr ; 
	  std::string m_dpndPmtr ;
	  double m_factor, m_addand ;
	} ;

	std::vector<ParameterConstraint> m_parameterConstraints ;

  } ;

}  //namespace

#endif // GUARD_FittingUtils_h
