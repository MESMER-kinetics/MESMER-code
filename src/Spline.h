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

	class Spline {

	public:

		Spline() : m_x(), m_y(), m_d2ydx2() {} ;
		~Spline() {} ;

		bool Initialize(std::vector<double> &x, std::vector<double> &y) ; 

		double Calculate(double x) const ;

		void Clear() { 
			m_x.clear(); 
			m_y.clear() ;
			m_d2ydx2.clear() ;
		} ;


	private:

		Spline operator=(Spline &spline) ;
		Spline(Spline &spline) ;

		std::vector<double> m_x ;
		std::vector<double> m_y ;
		std::vector<double> m_d2ydx2 ;

	} ;


}//namespacer mesmer

#endif // GUARD_Spline_h
