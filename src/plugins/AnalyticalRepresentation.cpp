//-------------------------------------------------------------------------------------------
//
// AnalyticalRepresentation.cpp
//
// Author: Struan Robertson
// Date:   21/Jul/2013
//
// This class implements the methods to that determine and analytical representation of a 
// master equation representation.
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
  class AnalyticalRepresentation : public CalcMethod, private FittingUtils
  {
  public:

	AnalyticalRepresentation(const char* id) : FittingUtils(), 
	  m_id(id),
	  m_format(),
	  m_reactionRef1(),
	  m_reactionRef2(),
	  m_TMax(0.0),
	  m_TMin(0.0),
	  m_CMax(0.0),
	  m_CMin(0.0),
	  m_NTpt(0),
	  m_NCpt(0),
	  m_RpTMin(0),
	  m_RpTMax(0),
	  m_lgCMin(0),
	  m_lgCMax(0),
	  m_ExpanSizeT(0),
	  m_ExpanSizeC(0)
	{ Register(); }

	virtual ~AnalyticalRepresentation() {}
	virtual const char* getID()  { return m_id; }

	//Read in data for this method from XML
	virtual bool ParseData(PersistPtr pp);

	//Function to do the work
	virtual bool DoCalculation(System* pSys);

  private:

	// Function to calculation a chebyshev polynomial of a given order i.
	inline double Cheb_poly(int order, const double arg) const { return cos((double(++order) - 1.0) * acos(arg)); }

	// Function to give constant coefficient of chebyshev expansion.
	double Coefficient(int i, int j);

	// Converts the Chebyshev gridpoints to back to Temperature/Concentrations.
	vector<double> Transform(const vector<double> &gridpoints, const double &max, const double &min);

	// Write Chebyshev coefficients in Cantera format.
	void writeCanteraFormat(const vector<vector<vector<double> > > &ChebyshevCoeff) const ;

	// Test the  Chebyshev representation.
	void testRepresentation(
	  const vector<vector<vector<double> > > &ChebyshevCoeff,
	  const vector<vector<vector<double> > > &PTgrid,
	  const vector<double> &Concentration,
	  const vector<double> &Temperature,
	  const vector<double> &CGrid,
	  const vector<double> &TGrid ) const ;

	// Utility function to pad output with blank space.
	string padText(size_t len) const ;

	// Utility function to underline text.
	string underlineText(const string& text) const ;

	const char* m_id;
	string m_format ;
	vector<string> m_reactionRef1 ;
	vector<string> m_reactionRef2 ;

	double m_TMax ;
	double m_TMin ;
	double m_CMax ;
	double m_CMin ;

	double m_RpTMin ;
	double m_RpTMax ;
	double m_lgCMin ;
	double m_lgCMax ;

	int m_NTpt ;
	int m_NCpt ;

	int m_ExpanSizeT;
	int m_ExpanSizeC;

  };

  ////////////////////////////////////////////////
  //Global instance
  AnalyticalRepresentation theAnalyticalRepresentation("analyticalRepresentation");
  ///////////////////////////////////////////////

  bool AnalyticalRepresentation::ParseData(PersistPtr pp)
  {
	/* Typical data
	<me:control>
	...
	<me:calcMethod xsi:type="me:analyticalRepresentation">
	<me:fittingTolerance>0.1</me:fittingTolerance>
	<me:fittingIterations>5</me:fittingIterations>
	</me:calcMethod>
	...
	</me:control>
	*/
	//Read in fitting parameters, or use values from defaults.xml.

	m_NTpt = pp->XmlReadInteger("me:chebNumTemp");
	m_NCpt = pp->XmlReadInteger("me:chebNumConc");

	m_TMax = pp->XmlReadDouble("me:chebMaxTemp");
	m_TMin = pp->XmlReadDouble("me:chebMinTemp");
	m_CMax = pp->XmlReadDouble("me:chebMaxConc");
	m_CMin = pp->XmlReadDouble("me:chebMinConc");

	m_RpTMin = 1.0 / m_TMin ;
	m_RpTMax = 1.0 / m_TMax ;
	m_lgCMin = log10(m_CMin);
	m_lgCMax = log10(m_CMax);

	m_ExpanSizeT = pp->XmlReadInteger("me:chebTExSize");
	m_ExpanSizeC = pp->XmlReadInteger("me:chebPExSize");

	// Check expansion is consistent with grid:
	if (m_ExpanSizeT > m_NTpt || m_ExpanSizeC > m_NCpt )
	  throw(std::runtime_error("Analytical Represention: Requested expansion coefficients exceed grid specificaton.")) ;

	PersistPtr ppARepRef = pp->XmlMoveTo("me:analyticalRepRef");
	ppARepRef = ppARepRef->XmlMoveTo("me:reaction");
	while(ppARepRef) {
	  const char* txt1 = ppARepRef->XmlReadValue("me:ref1") ;
	  m_reactionRef1.push_back(string(txt1));
	  const char* txt2 = ppARepRef->XmlReadValue("me:ref2") ;
	  m_reactionRef2.push_back(string(txt2));
	  ppARepRef = ppARepRef->XmlMoveTo("me:reaction");
	}

	return true;
  }

  bool AnalyticalRepresentation::DoCalculation(System* pSys)
  {
	//Do not output all the intermediate results to XML
	pSys->m_Flags.overwriteXmlAnalysis = true;

	// Use the same grain numbers for for all calcuations regardless of 
	// temperature (i.e. reduce the number of times micro-rates are caluclated).
	pSys->m_Flags.useTheSameCellNumber = true;

	// Warnings and less not sent to console.
	ChangeErrorLevel e(obError); 

	// First gets some points. Need to get chebyshev grid points and transform back to (T,P) condition set.

	vector<double> TGrid(m_NTpt);
	for (int i = 0; i < m_NTpt; ++i){
	  TGrid[i] = cos(((2.0*(i+1.0))-1.0)*M_PI / (2.0 * double(m_NTpt))) ;
	}

	vector<double> CGrid(m_NCpt);
	for (int i = 0; i < m_NCpt; ++i){
	  CGrid[i] = cos(((2.0*(i+1.0))-1.0)*M_PI / (2.0 * double(m_NCpt)));
	}

	vector<double> Temperature=Transform(TGrid, m_RpTMin, m_RpTMax);
	for (int i = 0; i< m_NTpt; ++i){
	  Temperature[i] = 1/Temperature[i];
	}

	vector<double> Concentration=Transform(CGrid, m_lgCMin, m_lgCMax);
	for (int i = 0; i< m_NCpt; ++i){
	  Concentration[i] = pow(10,(Concentration[i]));
	}

	vector<double> RateCoefficients(m_reactionRef1.size());

	// Get rate coefficients
	vector<double> v(m_reactionRef1.size());
	vector<vector<double> >v2(m_NCpt,v);
	vector<vector<vector<double> > > PTgrid(m_NTpt, v2) ;
	for (int m=0; m < m_NTpt; ++m ){
	  for (int n=0; n < m_NCpt; ++n ){
		for(size_t k=0; k < m_reactionRef1.size(); ++k){
		  pSys->calculate(Temperature[m], Concentration[n], m_reactionRef1, m_reactionRef2, RateCoefficients) ;
		  PTgrid[m][n][k] = RateCoefficients[k];
		} 
	  }
	}

	// Calculate chebyshev coefficients. Three indicies are required in order 
	// to calculate Chebyshev coefficients for each specified BW rate.
	vector<double> v3(m_reactionRef1.size());
	vector<vector<double> >v4(m_ExpanSizeC,v3);	
	vector<vector<vector<double> > > ChebyshevCoeff(m_ExpanSizeT,v4);
	for (int i=0; i < m_ExpanSizeT; ++i ){
	  for(int j=0; j < m_ExpanSizeC ; ++j ){
		for(size_t k=0; k < m_reactionRef1.size(); ++k){
		  for(int m = 0; m <m_NTpt ; ++m){
			for(int n = 0; n <m_NCpt ; ++n){
			  ChebyshevCoeff[i][j][k] += log10(PTgrid[m][n][k])*Cheb_poly(i, TGrid[m])*Cheb_poly(j, CGrid[n]);
			}
		  }				
		  ChebyshevCoeff[i][j][k] *= Coefficient(i, j) / (double(m_NTpt) * double(m_NCpt)) ;
		}			
	  }
	}

	// Print out table of Chebyshev coefficients for each BW rate specified.

	writeCanteraFormat(ChebyshevCoeff) ;

	// Test expansion.

    testRepresentation(ChebyshevCoeff, PTgrid, Concentration, Temperature, CGrid, TGrid) ;

	return true;
  }

	// Converts the Chebyshev gridpoints to back to Temperature/Concentrations.
  vector<double> AnalyticalRepresentation::Transform(const vector<double> &gridpoints, const double &max, const double &min ){
	vector<double> conditions(gridpoints.size());
	for(size_t i=0; i < gridpoints.size(); ++i ){
	  conditions[i] = (gridpoints[i] * (max - min) + max + min) / 2.0;
	}
	return conditions;
  }

	// Function to give constant coefficient of chebyshev expansion.
  double AnalyticalRepresentation::Coefficient( int i , int j ){
	double coeff=0.;
	if(i != 0 && j != 0){
	  coeff = 4.0;
	} else if((i==0 && j!=0) || ( j == 0 && i != 0)){
	  coeff = 2.0;
	} else if( i == 0 && j == 0 ){
	  coeff = 1.0;
	} else {
	  // this branch should never be executed
	}

	return coeff;
  }

  // Write Chebyshev coefficients in Cantera format.
  void AnalyticalRepresentation::writeCanteraFormat(const vector<vector<vector<double> > > &ChebyshevCoeff) const {

	string header("chebyshev_reaction(") ;
	string coeffs("coeffs=[") ;
	for (size_t k=0; k < m_reactionRef1.size(); ++k ){
	  string indent = padText(header.size()) ;
	  ctest << header << "' <=> '," << endl ;
	  ctest << indent << "Tmin=" << setw(6) << m_TMin << ", Tmax=" << setw(6) << m_TMax << "," << endl  ;
	  ctest << indent << "Pmin=(" << setw(6) << m_CMin << ", 'units'), Pmax=(" << setw(6) << m_CMax << ", 'units'),"  << endl;
	  ctest << indent << coeffs << "[" ;
	  indent += padText(coeffs.size()) ;
	  for (int i=0; i < m_ExpanSizeT; ++i ) {
		for(int j=0; j < m_ExpanSizeC ; ++j ) {
		  ctest << formatFloat(ChebyshevCoeff[i][j][k], 6, 14) ;
		  if (j < m_ExpanSizeC-1 ) 
			ctest << "," ;
		}
		if (i < m_ExpanSizeT-1) {
		  ctest << "]," << endl << indent << "[";
		} else {  
		  ctest << "]])" << endl << endl ;
		}
	  }
	}
  }

	// Test the  Chebyshev representation.
	void AnalyticalRepresentation::testRepresentation(
	  const vector<vector<vector<double> > > &ChebyshevCoeff,
	  const vector<vector<vector<double> > > &PTgrid,
	  const vector<double> &Concentration,
	  const vector<double> &Temperature,
	  const vector<double> &CGrid,
	  const vector<double> &TGrid ) const {

	  for (size_t k=0; k < m_reactionRef1.size(); ++k ){
	  ctest << "Comparison of fitted rate coefficients for reaction " << k << endl;

     ostringstream concText ;	  
	 string indent = padText(10) ;
	  concText << indent << " | " ;
	  for(size_t n = 0; n < m_NCpt ; ++n){
		concText << formatFloat(Concentration[n], 6, 22);
	  }
	  ctest << underlineText(concText.str()) ;

	  for (int m = 0; m <m_NTpt ; ++m) {
		ctest << setw(10) << Temperature[m] << " | " ;
		for (int n = 0; n <m_NCpt ; ++n) {
		  double ChebRate(0.);
		  for (int i=0; i < m_ExpanSizeT; ++i ){    
			for(int j=0; j < m_ExpanSizeC ; ++j ){
			  ChebRate += ChebyshevCoeff[i][j][k]*Cheb_poly(i, TGrid[m])*Cheb_poly(j, CGrid[n]);
			}
		  }	
		  ChebRate = pow(10.0,ChebRate);
		  ctest << setw(14) << ChebRate <<"/" << PTgrid[m][n][k];
		}
		ctest << endl;
	  }
	  ctest << endl;
	}
	}

  // Utility function to pad output with blank space.
  string AnalyticalRepresentation::padText(size_t len) const {
	ostringstream sstrdatum ;
	for (size_t i(0) ; i < len ; i++ ) 
	  sstrdatum << " " ;
	return sstrdatum.str() ;
  }

  // Utility function to underline text.
  string AnalyticalRepresentation::underlineText(const string& text) const {

	ostringstream sstrdatum ;
	sstrdatum << text << endl ;
	for (size_t i(0) ; i < text.size() ; i++ ) 
	  sstrdatum << "-" ;
	sstrdatum << endl ;

	return sstrdatum.str() ;

  }


}//namespace

