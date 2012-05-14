
//-------------------------------------------------------------------------------------------
//
// Spline.cpp
//
// Author: Struan Robertson
// Date:   14/May/2012
//
// This file contains the implementation of the Spline class.  
//
//-------------------------------------------------------------------------------------------

#include "Spline.h"

using namespace std ;

namespace mesmer
{

	bool Spline::Initialize(vector<double> &x, vector<double> &y) {

   // Ensure that that we are in a clean state before constructing a spline. 

    Clear() ;
		m_x = x ;
		m_y = y ;

		return true ;
	}

	double Spline::Calculate(double x) const {
		double y ;

		return y ;
	}

}