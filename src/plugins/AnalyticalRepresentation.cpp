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
	  m_precision(DOUBLE),
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

    PersistPtr m_pp;
    const char* m_id;
    string m_format ;
    Precision m_precision ;
    bool m_AllBathGases;

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
    string m_RateUnits;

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
    m_pp = pp; //save to write XML
    // Determine the required format, default to Cantera if not supplied.

    const char* txt = pp->XmlReadValue("me:format", optional) ;
    m_format = (txt) ? string(txt) : string("cantera") ;

    // Determine the required precision, default to double if not supplied.

    txt = pp->XmlReadValue("me:precision", optional) ;
	if (txt) {
	  m_precision = txtToPrecision(txt) ;
	};

    // Units can be an attribute on either <me:chebMaxConc> or <me:chebMinConc>;
    // if on neither the units in defaults.xml are used.
    PersistPtr pPUnits = pp->XmlMoveTo("me:chebMaxConc") ;
    if (!pPUnits)
      throw(std::runtime_error("<me:chebMaxConc> must be present"));

    txt = pPUnits->XmlReadValue("units", optional) ;
    if(!txt)
    {
      pPUnits = pp->XmlMoveTo("me:chebMinConc") ;
      txt = pPUnits->XmlReadValue("units") ; //or from defaults.xml
    }
    m_PUnits = string(txt) ;
    if(m_PUnits=="molescm-3" || m_PUnits=="mol/cc" || m_PUnits=="moles/cc")
      m_RateUnits = "cm3mole-1s-1";
    else
      m_RateUnits = "cm3molecule-1s-1";

    //Use rate units specified as attribute on format element if present
    PersistPtr ppFormat = pp->XmlMoveTo("me:format");
    if(ppFormat)
    {
      const char* rUnits = ppFormat->XmlReadValue("rateUnits", optional);
      if(rUnits)
        m_RateUnits = rUnits;
    }

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
    
    // Generate a representation for each bath gas mentioned anywhere in data file
    m_AllBathGases = pp->XmlMoveTo("me:chebDoForAllBathGases");

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

    set<string> bathGases;
    if(m_AllBathGases)
      pSys->getAllBathGases(bathGases);
    else // just the bath gas specified in <conditions>
      bathGases.insert(pSys->getMoleculeManager()->get_BathGasName());

    for(set<string>::iterator iter = bathGases.begin(); iter!=bathGases.end(); ++iter)
    {
      cinfo << "\nBath gas is " << *iter << endl;

      // Create a grid of temperature and concentration (in ppcc) values.
      // For temperature we must account for the fact that 1/Tmax < 1/Tmin.
      vector<double> Temperature=Transform(TGrid, m_RpTMin, m_RpTMax);
      for (size_t i(0); i < m_NTpt; ++i) {
        Temperature[i] = 1/Temperature[i] ;
      }

      vector<double> Concentration=Transform(CGrid, m_lgCMax, m_lgCMin);
      for (size_t i(0); i < m_NCpt; ++i) {
        Concentration[i] = pow(10,(Concentration[i])) ;
      }

      // Get rate coefficients.
	  double fctr = (m_RateUnits == "cm3mole-1s-1") ? Constants::AvogadroC : 1.0 ;
      m_reactions.clear() ;
      vector<CTpoint> CTGrid ;
      bool flag(true) ; // First pass flag used to get reaction details.
      vector<vector<double> > RCGrid;
      for (size_t i(0); i < m_NTpt; ++i) {
        double Temp = Temperature[i] ;
        for (size_t j(0); j < m_NCpt; ++j) {
          double Conc = getConvertedP(m_PUnits, Concentration[j], Temp) ;
          CTGrid.push_back(CTpoint(TGrid[i],CGrid[j])) ;
          map<string, double> phenRates ;
          pSys->calculate(Temp, Conc, m_precision, phenRates, m_TMax, *iter);
          vector<double> rate ; 
          map<string, double>::const_iterator itr;
          for (itr = phenRates.begin() ; itr != phenRates.end(); ++itr) {
            // Expand the string in phenRates to include all the reactants and products.
            pair<string,Reaction*> presult= pSys->getReactionManager()->getCompleteReactantsAndProducts(itr->first);
            if (flag) {
              m_reactions.push_back(presult.first);
            }
            Reaction *r = presult.second;
            double concExcessReactant = r ? r->get_concExcessReactant() : 0.0 ;
			double rateCoefficient    = r ? r->normalizeRateCoefficient(itr->second) : itr->second ;
            rate.push_back(concExcessReactant>0 ? rateCoefficient*fctr : rateCoefficient) ;
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
      ctest << "\nBath gas is " << *iter << endl;
      testRepresentation(ChebyshevCoeff, RCGrid, Concentration, Temperature, CGrid, TGrid) ;
    }
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

  // The rate constants can be specified as cm3molecule-1s-1 or cm3mole-1s-1
  // in a rateUnits attribute on the <format> element. If this is not present
  // The rate constant output is in cm3molecule-1s-1 unless the concentrations
  // are specified in cm3moles-1 when the rate constants are in are in cm3mole-1s-1.
  // The units are written to a units line in Cantera or a REACTIONS line in Chemkin.

  // Write Chebyshev coefficients in Cantera format.
  // See http://cantera.github.io/docs/sphinx/html/cti/reactions.html.
  void AnalyticalRepresentation::writeCanteraFormat(const vector<vector<vector<double> > > &ChebyshevCoeff, System* pSys) const {
    
    ostringstream sinfo;
    sinfo << "units(length = 'cm', quantity = '" 
          << ((m_RateUnits=="cm3mole-1s-1") ? "mole')" : "molecule')") << endl;
    string header("chebyshev_reaction(") ;
    string coeffs("coeffs=[") ;
    sinfo << endl ;
    for (size_t k=0; k < m_reactions.size(); ++k ){
      string indent = padText(header.size()) ;
      sinfo << header << "'" << m_reactions[k] << "'," << endl ;
      sinfo << indent << "Tmin="  << setw(6) << m_TMin << ", Tmax=" << setw(6) << m_TMax << "," << endl  ;
      sinfo << indent << "Pmin=(" << setw(6) << m_CMin << ", '" << m_PUnits << "'), " 
        << "Pmax=(" << setw(6) << m_CMax << ", '" << m_PUnits << "'), "  << endl;
      sinfo << indent << coeffs << "[" ;
      indent += padText(coeffs.size()) ;
      for (size_t i(0); i < m_ExpanSizeT; ++i ) {
        for(size_t j(0); j < m_ExpanSizeC ; ++j ) {
          sinfo << formatFloat(ChebyshevCoeff[i][j][k], 6, 14) ;
          if (j < m_ExpanSizeC-1 ) 
            sinfo << "," ;
        }
        if (i < m_ExpanSizeT-1) {
          sinfo << "]," << endl << indent << "[";
        } else {  
          sinfo << "]])" << endl << endl ;
        }
      }
    }
    //Write both to log file and XML file.
    cinfo << sinfo.str() << endl;
    PersistPtr ppp(m_pp);
    ppp->XmlWriteValueElement("me:representation", sinfo.str(), true); //CDATA
  }

  // Write Chebyshev coefficients in Chemkin format.
  void AnalyticalRepresentation::writeChemkinFormat(const vector<vector<vector<double> > > &ChebyshevCoeff, System* pSys) const {
    ostringstream sinfo;
    sinfo << endl ;
    sinfo << "REACTIONS" << ((m_RateUnits=="cm3mole-1s-1") ? " MOLES" : " MOLECULES") << '\n' << endl;
    for (size_t k=0; k < m_reactions.size(); ++k ) {
      sinfo << m_reactions[k] << " (+M)  1.00  0.0  0.0" << endl ;
      sinfo << "! Data generated by MESMER " << MESMER_VERSION << " on " << __DATE__ << "." << endl ;
      sinfo << "    TCHEB/ " << formatFloat(m_TMin, 6, 14) << formatFloat(m_TMax, 6, 14) << "/" << endl ;
      sinfo << "    PCHEB/ " << formatFloat(m_CMin, 6, 14) << formatFloat(m_CMax, 6, 14) << "/" << endl ;
      sinfo << "    CHEB/"   << setw(6) << m_ExpanSizeT << setw(6) << m_ExpanSizeC << "/" << endl ;
      for (size_t i(0); i < m_ExpanSizeT; ++i ) {
        sinfo << "    CHEB/" ;
        for (size_t j(0); j < m_ExpanSizeC ; ++j ) {
          sinfo << formatFloat(ChebyshevCoeff[i][j][k], 6, 14) ;
        }
        sinfo << "/" << endl ;
      }
      sinfo << endl ;
    }
    //Write both to log file and XML file.
    cinfo << sinfo.str() << endl;
    PersistPtr ppp(m_pp);
    ppp->XmlWriteValueElement("me:representation", sinfo.str(), true); //CDATA
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

