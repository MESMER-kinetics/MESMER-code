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

		~FittingUtils() {};

	protected:

		// Get the current location.
		void GetLocation(vector<double> &loc) const;

		// Set the current location.
		void SetLocation(vector<double> &loc) const;

		// Check that the a point falls within the limits defined by the user.
		bool CheckBounds(const vector<double> &A) const;

		// Numerical derivatives.
		void NumericalDerivatives(System* pSys, vector<double> &residuals, double delta, vector<double> &gradient, qdMatrix &hessian) const;

		// Write out the results and statistics of the fit. 
		void ResultsAndStatistics(System* pSys, qdMatrix &hessian) const;

		// Get correlated distribution from fitting hessian matrix
		void CorellatedDistribution(System* pSys, qdMatrix hessian, size_t &seqSize);

	};

}  //namespace

#endif // GUARD_FittingUtils_h
