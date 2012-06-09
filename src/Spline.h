#ifndef GUARD_Spline_h
#define GUARD_Spline_h

//-------------------------------------------------------------------------------------------
//
// Spline.h
//
// Author: Struan Robertson
// Date:   14/May/2012
//
// This header file contains the declaration of the Spline class.  
//
//-------------------------------------------------------------------------------------------

#include <vector> 

namespace mesmer
{

    // Criterion for determining if the spline is natural at the end points.

    static const double naturalLimit = 0.99e30 ;

	class Spline {

	public:

		Spline() : m_x(), m_y(), m_d2ydx2(), m_splineDefined(false) {} ;
		~Spline() {} ;

		bool Initialize(const std::vector<double> &x, const std::vector<double> &y, double lower = naturalLimit, double upper = naturalLimit) ; 
		bool Initialize(const std::vector<std::pair<double,double> >* data,         double lower = naturalLimit, double upper = naturalLimit) ; 

		double Calculate(double x) const ;

	private:

		Spline operator=(Spline &spline) ;
		Spline(Spline &spline) ;

		void Clear() { 
			m_x.clear(); 
			m_y.clear() ;
			m_d2ydx2.clear() ;
			m_splineDefined = false ;
		} ;

		std::vector<double> m_x ;
		std::vector<double> m_y ;
		std::vector<double> m_d2ydx2 ;

		bool m_splineDefined ;

	} ;


}//namespacer mesmer

#endif // GUARD_Spline_h
