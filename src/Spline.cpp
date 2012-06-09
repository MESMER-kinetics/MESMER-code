
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
#include <iostream>

using namespace std ;

namespace mesmer
{

  bool Spline::Initialize(const vector<double> &x, const vector<double> &y, double lower, double upper) {

	// Ensure that that we are in a clean state before constructing a spline. 

	Clear() ;

	// Make sure arrays are the same size, have sufficient points and the points are unique.

	if (x.size() != y.size()) 
	  return false ;

	if (x.size() < 3) {
	  cerr << "There are insufficient points to construct a spline." ;
	  return false ;
	}

	for (size_t i(0) ; i < x.size() - 2 ; i++) {
	  if (x[i] == x[i+1]) {
		cerr << "There are knots at the same points." ;
		return false ;
	  }
	}

	m_x = x ;
	m_y = y ;

	size_t nsize = x.size() ;
	m_d2ydx2.resize(nsize,0.0) ;
	vector<double> wrk(nsize,0.0) ;

	// Lower boundary condition:

	m_d2ydx2[0] = wrk[0] = 0.0 ;
	if (lower < naturalLimit) {
	  m_d2ydx2[0] = -0.5 ;
	  wrk[0]      = (3.0/(m_x[1]-m_x[0]))*((m_y[1]-m_y[0])/(m_x[1]-m_x[0]) - lower) ;
	}

	// Decomposition of the tridiagonal matrix.

	for ( size_t i(1) ; i < nsize-1 ; i++) {
	  double sig  = (m_x[i] - m_x[i-1])/(m_x[i+1] - m_x[i-1]) ;
	  double p    = sig*m_d2ydx2[i-1] + 2.0 ; 
	  m_d2ydx2[i] = (sig - 1.0)/p ;
	  wrk[i]      = (m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]) - (m_y[i] - m_y[i-1])/(m_x[i] - m_x[i-1]) ;
	  wrk[i]      = (6.0*wrk[i]/(m_x[i+1]-m_x[i-1]) - sig*wrk[i-1])/p ;
	}

	// Upper boundary condition:

	double qn(0.0), un(0.0) ;
	if (upper < naturalLimit) {
	  qn = 0.5 ;
	  double h = m_x[nsize-1]-m_x[nsize-2] ;
	  un = (3.0/h)*(upper - (m_y[nsize-1]-m_y[nsize-2])/h) ;
	}

	m_d2ydx2[nsize-1] = (un-qn*wrk[nsize-2])/(qn*m_d2ydx2[nsize-2] + 1.0) ;

	// Determine second derivatives by back substitution.

	for (size_t k(nsize-1) ; k > 0 ; k--) {
	  m_d2ydx2[k-1] = m_d2ydx2[k-1]*m_d2ydx2[k] + wrk[k-1] ;
	}

	// Spline is defined.

	m_splineDefined = true ;

	return m_splineDefined ;
  }

  bool Spline::Initialize(const std::vector<std::pair<double,double> >* pdata,  double lower, double upper) {

	size_t nsize = pdata->size() ;
	vector<double> x(nsize, 0.0),  y(nsize, 0.0) ;
	for (size_t i(0) ; i < nsize ; i++ ) {
	  x[i] = (*pdata)[i].first ;
	  y[i] = (*pdata)[i].second ;
	}
	return Initialize(x, y, lower, upper) ;
  }

  double Spline::Calculate(double x) const {
	double y(0.0) ;

	size_t nsize = m_x.size() ;
	size_t klo   = 0 ;
	size_t kup   = nsize-1 ;

	// Check that we have a spline defined. 

	if (!m_splineDefined) {
	  cerr << "Spline is not defined" ;
	  return y ;
	}

	while (kup - klo > 1) {
	  size_t k = (kup + klo)/2 ;
	  if (m_x[k] > x) {
		kup = k ;
	  } else {
		klo = k ;
	  }
	}
	double h = m_x[kup] - m_x[klo] ;
	double a = (m_x[kup] - x)/h ;
	double b = (x - m_x[klo])/h ;

	y = a*m_y[klo] + b*m_y[kup]  + ((a*a*a -a)*m_d2ydx2[klo] + (b*b*b - b)*m_d2ydx2[kup])*h*h/6.0 ;

	return y ;
  }

}