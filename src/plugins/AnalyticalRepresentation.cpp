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
#include <sstream>

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
      m_TMax(0.0),
      m_TMin(0.0),
      m_CMax(0.0),
      m_CMin(0.0),
      m_RpTMin(0),
      m_RpTMax(0),
      m_lgCMin(0),
      m_lgCMax(0),
      m_NTpt(0),
      m_NCpt(0),
      m_PUnits(),
      m_ExpanSizeT(0),
      m_ExpanSizeC(0),
      m_reactions()
    { Register(); }

    virtual ~AnalyticalRepresentation() {}
    virtual const char* getID()  { return m_id; }

    //Read in data for this method from XML
    virtual bool ParseData(PersistPtr pp);

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:

    typedef pair<double, double> CTpoint ;

    // Function to calculation a chebyshev polynomial of a given order i.
    inline double Cheb_poly(int order, const double arg) const { return cos((double(++order) - 1.0) * acos(arg)); }

    // Function to give constant coefficient of chebyshev expansion.
    double Coefficient(int i, int j);

    // Converts the Chebyshev gridpoints to back to Temperature/Concentrations.
    vector<double> Transform(const vector<double> &gridpoints, const double &max, const double &min);

    // Write Chebyshev coefficients in Cantera format.
    void writeCanteraFormat(const vector<vector<vector<double> > > &ChebyshevCoeff, System* pSys) const ;

    // Write Chebyshev coefficients in Chemkin format.
    void writeChemkinFormat(const vector<vector<vector<double> > > &ChebyshevCoeff, System* pSys) const ;

    // Test the  Chebyshev representation.
    void testRepresentation(
      const vector<vector<vector<double> > > &ChebyshevCoeff,
      const vector<vector<double> > &RCGrid,
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

    double m_TMax ;
    double m_TMin ;
    double m_CMax ;
    double m_CMin ;

    double m_RpTMin ;
    double m_RpTMax ;
    double m_lgCMin ;
    double m_lgCMax ;

    size_t m_NTpt ;
    size_t m_NCpt ;
    string m_PUnits ;

    size_t m_ExpanSizeT;
    size_t m_ExpanSizeC;

    vector<string> m_reactions ;

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

    // Determine the required format, default to Cantera if not supplied.
    const char* txt = pp->XmlReadValue("me:format", optional) ;
    m_format = (txt) ? string(txt) : string("cantera") ;

    // Units can be an attribute on either <me:chebMaxConc> or <me:chebMinConc>;
    // if on neither the units in defaults.xml are used.
    PersistPtr pPUnits = pp->XmlMoveTo("me:chebMaxConc") ;
    txt = pPUnits->XmlReadValue("units", optional) ;
    if(!txt)
    {
      pPUnits = pp->XmlMoveTo("me:chebMinConc") ;
      txt = pPUnits->XmlReadValue("units") ; //or from defaults.xml
    }
    m_PUnits = string(txt) ;

	// Chemkin only supports pressure units of atm.
	if ( m_format == "chemkin" && m_PUnits != "atm")
	   throw(std::runtime_error("Chemkin only supports pressure units of atm.")) ;

    //Read in fitting parameters, or use values from defaults.xml.
    m_NTpt = pp->XmlReadInteger("me:chebNumTemp");
    m_NCpt = pp->XmlReadInteger("me:chebNumConc");

    m_TMax = pp->XmlReadDouble("me:chebMaxTemp");
    m_TMin = pp->XmlReadDouble("me:chebMinTemp");
    if (m_TMin > m_TMax)
      throw(std::runtime_error("Analytical Represention: Max. Temp. less than Min. Temp.")) ;
    m_CMax = pp->XmlReadDouble("me:chebMaxConc");
    m_CMin = pp->XmlReadDouble("me:chebMinConc");
    if (m_CMin > m_CMax)
      throw(std::runtime_error("Analytical Represention: Max. Pres. less than Min. Pres.")) ;

    m_RpTMin = 1.0 / m_TMin ;
    m_RpTMax = 1.0 / m_TMax ;
    m_lgCMin = log10(m_CMin);
    m_lgCMax = log10(m_CMax);

    m_ExpanSizeT = pp->XmlReadInteger("me:chebTExSize");
    m_ExpanSizeC = pp->XmlReadInteger("me:chebPExSize");

    // Check expansion is consistent with grid:
    if (m_ExpanSizeT > m_NTpt || m_ExpanSizeC > m_NCpt )
      throw(std::runtime_error("Analytical Represention: Requested expansion coefficients exceed grid specificaton.")) ;

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
    for (size_t i(0); i < m_NTpt; ++i){
      TGrid[i] = cos(((2.0*(i+1.0))-1.0)*M_PI / (2.0 * double(m_NTpt))) ;
    }

    vector<double> CGrid(m_NCpt);
    for (size_t i(0); i < m_NCpt; ++i){
      CGrid[i] = cos(((2.0*(i+1.0))-1.0)*M_PI / (2.0 * double(m_NCpt)));
    }

    // Create a grid of temperature and concentration (in ppcc) values.
    vector<double> Temperature=Transform(TGrid, m_RpTMin, m_RpTMax);
    for (size_t i(0); i < m_NTpt; ++i) {
      Temperature[i] = 1/Temperature[i] ;
    }

    vector<double> Concentration=Transform(CGrid, m_lgCMin, m_lgCMax);
    for (size_t i(0); i < m_NCpt; ++i) {
      Concentration[i] = pow(10,(Concentration[i])) ;
    }

    // Get rate coefficients.
    m_reactions.clear() ;
    vector<CTpoint> CTGrid ;
    bool flag(true) ;
    vector<vector<double> > RCGrid;
    for (size_t i(0); i < m_NTpt; ++i) {
      double Temp = Temperature[i] ;
      for (size_t j(0); j < m_NCpt; ++j) {
        double Conc = getConvertedP(m_PUnits, Concentration[j], Temp) ;
        CTGrid.push_back(CTpoint(TGrid[i],CGrid[j])) ;
        map<string, double> phenRates ;
        pSys->calculate(Temp, Conc, phenRates, m_TMax);
        vector<double> rate ;
        for (map<string, double>::const_iterator itr = phenRates.begin() ; itr != phenRates.end(); ++itr) {
          double concExcessReactant(0);
          if (flag) {
            //Expand the string in phenRates to include all the reactants and products
            pair<string,Reaction*> presult= pSys->getReactionManager()->getCompleteReactantsAndProducts
              (itr->first);
            Reaction* r = presult.second;
            m_reactions.push_back(presult.first) ;

            // concExcessReactant is zero for first order reactions and when reaction not found
            concExcessReactant = r ? 
              r->get_concExcessReactant()/getConvertedP(m_PUnits, 1.0, Temp) : 0.0;
          }
          //Adjust rate if a second order reaction
          rate.push_back(concExcessReactant>0 ? itr->second/concExcessReactant : itr->second) ;
        }
        flag = false ;
        RCGrid.push_back(rate) ;
      }
    }

    // Calculate chebyshev coefficients. Three indicies are required in order 
    // to calculate Chebyshev coefficients for each specified BW rate.
    vector<vector<double> > v(m_ExpanSizeC, vector<double>(m_reactions.size(), 0.0));	
    vector<vector<vector<double> > > ChebyshevCoeff(m_ExpanSizeT,v);
    for (size_t i(0); i < m_ExpanSizeT ; ++i ) {
      for (size_t j(0); j < m_ExpanSizeC ; ++j ) {
        for (size_t k(0); k < m_reactions.size() ; ++k) {
          for (size_t m(0); m < RCGrid.size() ; ++m ) {
            // The absolute value is taken below as small rate coefficients occasionally computed to be negative.
            ChebyshevCoeff[i][j][k] += log10(fabs(RCGrid[m][k]))*Cheb_poly(i, CTGrid[m].first)*Cheb_poly(j, CTGrid[m].second);
          }				
          ChebyshevCoeff[i][j][k] *= Coefficient(i, j) / (double(m_NTpt) * double(m_NCpt)) ;
        }			
      }
    }

    // Print out table of Chebyshev coefficients for each BW rate specified.

    if (m_format == string("cantera")) {
      writeCanteraFormat(ChebyshevCoeff, pSys) ;
    } else if (m_format == string("chemkin")) {
      writeChemkinFormat(ChebyshevCoeff, pSys) ;
    } else {
      writeCanteraFormat(ChebyshevCoeff, pSys) ;
    }

    // Test expansion.

    testRepresentation(ChebyshevCoeff, RCGrid, Concentration, Temperature, CGrid, TGrid) ;

    return true;
  }

  // Converts the Chebyshev gridpoints back to Temperature/Concentrations.
  vector<double> AnalyticalRepresentation::Transform(const vector<double> &gridpoints, const double &max, const double &min ){
    vector<double> conditions(gridpoints.size());
    for(size_t i(0); i < gridpoints.size(); ++i ){
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
  // See http://cantera.github.io/docs/sphinx/html/cti/reactions.html.
  void AnalyticalRepresentation::writeCanteraFormat(const vector<vector<vector<double> > > &ChebyshevCoeff, System* pSys) const {

    string header("chebyshev_reaction(") ;
    string coeffs("coeffs=[") ;
    cinfo << endl ;
    for (size_t k=0; k < m_reactions.size(); ++k ){
      string indent = padText(header.size()) ;
      cinfo << header << "'" << m_reactions[k] << "'," << endl ;
      cinfo << indent << "Tmin="  << setw(6) << m_TMin << ", Tmax=" << setw(6) << m_TMax << "," << endl  ;
      cinfo << indent << "Pmin=(" << setw(6) << m_CMin << ", '" << m_PUnits << "'), " 
        << "Pmax=(" << setw(6) << m_CMax << ", '" << m_PUnits << "'), "  << endl;
      cinfo << indent << coeffs << "[" ;
      indent += padText(coeffs.size()) ;
      for (size_t i(0); i < m_ExpanSizeT; ++i ) {
        for(size_t j(0); j < m_ExpanSizeC ; ++j ) {
          cinfo << formatFloat(ChebyshevCoeff[i][j][k], 6, 14) ;
          if (j < m_ExpanSizeC-1 ) 
            cinfo << "," ;
        }
        if (i < m_ExpanSizeT-1) {
          cinfo << "]," << endl << indent << "[";
        } else {  
          cinfo << "]])" << endl << endl ;
        }
      }

      /*XML output under <reaction>
      <rateConstant format="Cantera" units="PPCC">
        <![CDATA[chebyshev_reaction('CH3OCH2 => IM1',etc]]>
      </rateConstant>
      */
      std::stringstream ss;
      ss << cinfo.rdbuf();
      string s = ss.str();
    }
  }

  // Write Chebyshev coefficients in Chemkin format.
  void AnalyticalRepresentation::writeChemkinFormat(const vector<vector<vector<double> > > &ChebyshevCoeff, System* pSys) const {

    cinfo << endl ;
    for (size_t k=0; k < m_reactions.size(); ++k ) {
      cinfo << m_reactions[k] << " (+M)  1.00  0.0  0.0" << endl ;
	  cinfo << "! Data generated by MESMER " << MESMER_VERSION << " on " << __DATE__ << "." << endl ;
      cinfo << "    TCHEB/ " << formatFloat(m_TMin, 6, 14) << formatFloat(m_TMax, 6, 14) << "/" << endl ;
      cinfo << "    PCHEB/ " << formatFloat(m_CMin, 6, 14) << formatFloat(m_CMax, 6, 14) << "/" << endl ;
      cinfo << "    CHEB/"   << setw(6) << m_ExpanSizeT << setw(6) << m_ExpanSizeC << "/" << endl ;
      for (size_t i(0); i < m_ExpanSizeT; ++i ) {
		cinfo << "    CHEB/" ;
        for (size_t j(0); j < m_ExpanSizeC ; ++j ) {
          cinfo << formatFloat(ChebyshevCoeff[i][j][k], 6, 14) ;
        }
        cinfo << "/" << endl ;
      }
      cinfo << endl ;

      /*XML output under <reaction>
      <rateConstant format="Cantera" units="PPCC">
        <![CDATA[chebyshev_reaction('CH3OCH2 => IM1',etc]]>
      </rateConstant>
      */
      std::stringstream ss;
      ss << cinfo.rdbuf();
      string s = ss.str();
    }
  }

  // Test the  Chebyshev representation.
  void AnalyticalRepresentation::testRepresentation(
    const vector<vector<vector<double> > > &ChebyshevCoeff,
    const vector<vector<double> > &RCGrid,
    const vector<double> &Concentration,
    const vector<double> &Temperature,
    const vector<double> &CGrid,
    const vector<double> &TGrid ) const {

      for (size_t k=0; k < m_reactions.size() ; ++k ){
        ctest << "Comparison of fitted rate coefficients for reaction " << m_reactions[k] << endl;

        ostringstream concText ;	  
        string indent = padText(10) ;
        concText << indent << " | " ;
        for(size_t n = 0; n < m_NCpt ; ++n){
          concText << formatFloat(Concentration[n], 6, 22);
        }
        ctest << underlineText(concText.str()) ;

        for (size_t m(0), idx(0) ; m <m_NTpt ; ++m) {
          ctest << setw(10) << Temperature[m] << " | " ;
          for (size_t n(0); n <m_NCpt ; ++n, ++idx) {
            double ChebRate(0.);
            for (size_t i(0); i < m_ExpanSizeT; ++i) {    
              for(size_t j(0); j < m_ExpanSizeC ; ++j ){
                ChebRate += ChebyshevCoeff[i][j][k]*Cheb_poly(i, TGrid[m])*Cheb_poly(j, CGrid[n]);
              }
            }	
            ctest << setw(14) << pow(10.0,ChebRate) <<"/" << RCGrid[idx][k];
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

