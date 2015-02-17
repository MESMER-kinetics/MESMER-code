//-------------------------------------------------------------------------------------------
//
// SensitivityAnalysis.cpp
//
// Author: Struan Robertson
// Date:   25/Nov/2012
//
// This class implements sensitivity analysis algorithm used to identify this most important
// parameters from a fit. It is based on the Li et al, Chem. Eng. Sci. Vol. 57, 4445 (2002) and 
// guided by the matlab implementation of Tilo Ziehn, particularly the methods sub_alpha_1st.m,
// sub_beta_2nd.m and sub_sensitivity_indices.m.
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"
#include "FittingUtils.h"
#include "../Sobol.h"

namespace {

  // The first 15 shifted Legendre polynomials

  double ortho_poly_1  (const double &x) {return sqrt(3.0)*(2.0*x - 1.0); } 

  double ortho_poly_2  (const double &x) {return sqrt(5.0)*(6.0*x*x - 6.0*x + 1.0); }

  double ortho_poly_3  (const double &x) { 
    double x2 =  x*x ;
    double x3 = x2*x ;
    return sqrt(7.)*(20.0*x3 - 30.0*x2 + 12.0*x - 1.); } 

  double ortho_poly_4  (const double &x) {
    double x2 =  x*x;
    double x3 = x2*x ;
    double x4 = x3*x ;
    return 3.0*(70.*x4 - 140.*x3 + 90.0*x2 - 20.0*x + 1.); }

  double ortho_poly_5  (const double &x) {
    double x2 =  x*x ;
    double x3 = x2*x ;
    double x4 = x3*x ;
    double x5 = x4*x ;
    return 630.*sqrt(11.)*((2./5.)*x5-x4+(8./9.)*x3-(1./3.)*x2+(1./21.)*x-(1./630.)); }

  double ortho_poly_6  (const double &x) {
    double x2 =  x*x ;
    double x3 = x2*x ;
    double x4 = x3*x ;
    double x5 = x4*x ;
    double x6 = x5*x ;
    return 3150.*sqrt(13.)*((22./75.)*x6-(22./25.)*x5+x4-(8./15.)*x3+(2./15.)*x2-(1./75.)*x+(1./3150.)); }

  double ortho_poly_7  (const double &x) {
    double x2 =  x*x ;
    double x3 = x2*x ;
    double x4 = x3*x ;
    double x5 = x4*x ;
    double x6 = x5*x ;
    double x7 = x6*x ;
    return 16632.*sqrt(15.)*((13./63.)*x7-(13./18.)*x6+x5-(25./36.)*x4+(25./99.)*x3-(1./22.)*x2+(1./297.)*x-(1./16632.)); } 

  double ortho_poly_8  (const double &x) {
    double x2 =  x*x ;
    double x3 = x2*x ;
    double x4 = x3*x ;
    double x5 = x4*x ;
    double x6 = x5*x ;
    double x7 = x6*x ;
    double x8 = x7*x ;
    return 84084.*sqrt(17.)*((15./98.)*x8 - (30./49.)*x7 + x6 - (6./7.)*x5 + (75./182.)*x4 
      - (10./91.)*x3 + (15./1001.)*x2 - (6./7007.)*x + (1./84084.)); } 

  double ortho_poly_9  (const double &x) {
    double x2 =  x*x ;
    double x3 = x2*x ;
    double x4 = x3*x ;
    double x5 = x4*x ;
    double x6 = x5*x ;
    double x7 = x6*x ;
    double x8 = x7*x ;
    double x9 = x8*x ;
    return 420420.*sqrt(19.)*((17./147.)*x9 - (51./98.)*x8 + (48./49.)*x7 - x6 + (3./5.)*x5 - (3./14.)*x4 
      +(4./91.)*x3 - (3./637.)*x2 + (3./14014.)*x - (1./420420.)); } 

  double ortho_poly_10 (const double &x) {
    double x2  =  x*x ;
    double x3  = x2*x ;
    double x4  = x3*x ;
    double x5  = x4*x ;
    double x6  = x5*x ;
    double x7  = x6*x ;
    double x8  = x7*x ;
    double x9  = x8*x ;
    double x10 = x9*x ;
    return sqrt(21.)*(184756.*x10 - 923780.*x9 + 1969110.*x8 - 2333760.*x7 + 
      1681680.*x6 - 756756.*x5 + 210210.*x4 - 34320.*x3 + 2970.*x2 - 110.*x + 1.); } 

  double ortho_poly_11 (const double &x) {
    double x2  =   x*x ;
    double x3  =  x2*x ;
    double x4  =  x3*x ;
    double x5  =  x4*x ;
    double x6  =  x5*x ;
    double x7  =  x6*x ;
    double x8  =  x7*x ;
    double x9  =  x8*x ;
    double x10 =  x9*x ;
    double x11 = x10*x ;
    return sqrt(23.)*(705432.*x11 - 3879876.*x10 + 9237800.*x9 - 12471030.*x8 + 10501920.*x7 - 5717712.*x6 + 
      2018016.*x5 - 450450.*x4 + 60060.*x3 - 4290.*x2 + 132.*x - 1.); } 

  double ortho_poly_12 (const double &x) {
    double x2  =   x*x ;
    double x3  =  x2*x ;
    double x4  =  x3*x ;
    double x5  =  x4*x ;
    double x6  =  x5*x ;
    double x7  =  x6*x ;
    double x8  =  x7*x ;
    double x9  =  x8*x ;
    double x10 =  x9*x ;
    double x11 = x10*x ;
    double x12 = x11*x ;
    return 5.*(2704156.*x12 - 16224936.*x11 + 42678636.*x10 -64664600.*x9 + 62355150.*x8 - 39907296.*x7 + 
      17153136.*x6 - 4900896.*x5 + 900900.*x4 - 100100.*x3 + 6006.*x2 - 156.*x + 1.); }

  double ortho_poly_13 (const double &x) {
    double x2  =   x*x ;
    double x3  =  x2*x ;
    double x4  =  x3*x ;
    double x5  =  x4*x ;
    double x6  =  x5*x ;
    double x7  =  x6*x ;
    double x8  =  x7*x ;
    double x9  =  x8*x ;
    double x10 =  x9*x ;
    double x11 = x10*x ;
    double x12 = x11*x ;
    double x13 = x12*x ;
    return sqrt(27.)*(10400600.*x13 - 67603900.*x12 + 194699232.*x11 - 327202876.*x10 + 355655300.*x9 - 261891630.*x8 + 
      133024320.*x7 - 46558512.*x6 + 11027016.*x5 - 1701700.*x4 + 160160.*x3 - 8190.*x2 + 182.*x - 1.); } 

  double ortho_poly_14 (const double &x) {
    double x2  =   x*x ;
    double x3  =  x2*x ;
    double x4  =  x3*x ;
    double x5  =  x4*x ;
    double x6  =  x5*x ;
    double x7  =  x6*x ;
    double x8  =  x7*x ;
    double x9  =  x8*x ;
    double x10 =  x9*x ;
    double x11 = x10*x ;
    double x12 = x11*x ;
    double x13 = x12*x ;
    double x14 = x13*x ;
    return sqrt(29.)*(40116600.*x14 - 280816200.*x13 + 878850700.*x12 - 1622493600.*x11 + 1963217256.*x10 - 1636014380.*x9 +
      960269310.*x8 - 399072960.*x7 + 116396280.*x6 - 23279256.*x5 + 3063060.*x4 - 247520.*x3 + 10920.*x2 - 210.*x + 1.); }

  double ortho_poly_15 (const double &x) {
    double x2  =   x*x ;
    double x3  =  x2*x ;
    double x4  =  x3*x ;
    double x5  =  x4*x ;
    double x6  =  x5*x ;
    double x7  =  x6*x ;
    double x8  =  x7*x ;
    double x9  =  x8*x ;
    double x10 =  x9*x ;
    double x11 = x10*x ;
    double x12 = x11*x ;
    double x13 = x12*x ;
    double x14 = x13*x ;
    double x15 = x14*x ;
    return sqrt(32.)*(155117520.*x15 - 1163381400.*x14 + 3931426800.*x13 - 7909656300.*x12 + 10546208400.*x11 - 9816086280.*x10 +
      6544057520.*x9 - 3155170590.*x8 + 1097450640.*x7 - 271591320.*x6 + 46558512.*x5 - 5290740.*x4 + 371280.*x3 - 14280.*x2 + 240.*x - 1.); }
}

namespace mesmer
{
  class SensitivityAnalysis : public CalcMethod, private FittingUtils
  {
  public:

    SensitivityAnalysis(const char* id) : m_id(id), 
      m_nVar(0), 
      m_nOut(0), 
      m_maxIterations(0), 
      m_bGenerateData(true), 
      m_delta(0),
      m_order(1) { 
        Register();

        m_slpMap[1]  = ortho_poly_1  ;
        m_slpMap[2]  = ortho_poly_2  ;
        m_slpMap[3]  = ortho_poly_3  ;
        m_slpMap[4]  = ortho_poly_4  ;
        m_slpMap[5]  = ortho_poly_5  ;
        m_slpMap[6]  = ortho_poly_6  ;
        m_slpMap[7]  = ortho_poly_7  ;
        m_slpMap[8]  = ortho_poly_8  ;
        m_slpMap[9]  = ortho_poly_9  ;
        m_slpMap[10] = ortho_poly_10 ;
        m_slpMap[11] = ortho_poly_11 ;
        m_slpMap[12] = ortho_poly_12 ;
        m_slpMap[13] = ortho_poly_13 ;
        m_slpMap[14] = ortho_poly_14 ;
        m_slpMap[15] = ortho_poly_15 ;
    }

    virtual ~SensitivityAnalysis() {}
    virtual const char* getID()  { return m_id; }
    virtual bool ParseData(PersistPtr pp);

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:

    // This method does a complete sensitivity analysis.
    bool DoCalculationNew(System* pSys); 

    // This method generates data for external analysis.
    bool DoCalculationOld(System* pSys); 

    // This methods writes out the results of a sensitivity analysis.
    bool WriteOutAnalysis(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &alpha, double Temperature, double Concentration) ;

    // This method calculates sensitivity indicies.
    bool sensitivityIndicies(const vector<double> &f0, const vector<double> &alpha, const vector<double> &beta); 

    // Methods for generating values of shifted Legendre polynomials.
    typedef double (*slp)(const double &x) ;
    map<size_t, slp> m_slpMap ;
    double ShiftedLegendre(size_t order, const double x) {double f = m_slpMap[order](x) ; 
    return f ; } ;

    // This method generates the column header for the output tables.
    string columnHeader() const ;

    const char* m_id;

    size_t m_nVar ;                // Dimension of analysis - number of inputs.
    size_t m_nOut ;                // Dimension of analysis - number of outputs.

    size_t m_maxIterations ;
    bool   m_bGenerateData ;

    vector<double> m_delta ;

    size_t m_order ;               // The order of the HDMR analysis to use.

    vector<vector<double> > m_Di ;  // First order sensitivities.
    vector<vector<double> > m_Dij ; // Second order sensitivities.

  };

  ////////////////////////////////////////////////
  //Global instance
  SensitivityAnalysis theSensitivityAnalysis("SensitivityAnalysis");
  ///////////////////////////////////////////////

  bool SensitivityAnalysis::ParseData(PersistPtr pp)
  {
    // Read in sensitivity analysis parameters, or use values from defaults.xml.
    m_maxIterations = pp->XmlReadInteger("me:SensitivityAnalysisIterations");
    m_bGenerateData = pp->XmlReadBoolean("me:SensitivityGenerateData");
    m_order         = pp->XmlReadInteger("me:SensitivityAnalysisOrder");

    return true;
  }

  bool SensitivityAnalysis::DoCalculation(System* pSys) {
    return (m_bGenerateData) ? DoCalculationOld(pSys) :	DoCalculationNew(pSys) ;
  }

  bool SensitivityAnalysis::DoCalculationOld(System* pSys)
  {

    m_nVar = Rdouble::withRange().size() ;

    if (m_nVar < 1) { 
      cerr << "Sensitivity analysis requries at least one range variable to be set." << endl;
      return false ;
    }

    //Read variable uncertainties from range
    for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
      Rdouble var = *Rdouble::withRange()[iVar] ;
      double lower = var.get_lower();
      double upper = var.get_upper();
      double middle = upper-((upper - lower) / 2.0);
      m_delta.push_back(abs((middle - lower) / middle));
    }

    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are calculated).
    pSys->m_Flags.useTheSameCellNumber = true;

    // Uncomment to enable ctest output during fitting. Or use -w5 option in command.
    //ChangeErrorLevel e(obDebug); 

    //Default is to disable ctest during fitting. Restored when leaving this function.
    StopCTestOutput stop(true) ;

    //
    // Begin by finding the starting point chi-squared value.
    //

    vector<double> currentLocation(m_nVar,0.0) ; 
    vector<double> newLocation(m_nVar,0.0) ; 

    GetLocation(currentLocation) ;

    // Invoke SetLocation to catch any constrained parameters.
    SetLocation(currentLocation) ;

    vector<double> Temperature ;
    vector<double> Concentration ;
    pSys->getConditionsManager()->getConditions (Temperature, Concentration) ;

    // Instantiate a random vector generator.
    Sobol sobol ;

    // Loop over condiditons. 
    size_t nConditions = Temperature.size() ;
    for (size_t nCnd(0) ; nCnd < nConditions ; nCnd++) {

      // String stream to hold results. 
      stringstream sensitivityTable ;

      // Write out table header.

      sensitivityTable << endl ;
      sensitivityTable << "Sensitivity Table" << endl ;
      sensitivityTable << "  Temperature:   " << formatFloat(Temperature[nCnd],   5, 15) << " K"    << endl ;
      sensitivityTable << "  Concentration: " << formatFloat(Concentration[nCnd], 5, 15) << " ppcc" << endl ;
      sensitivityTable << endl ;
      for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
        Rdouble var = *Rdouble::withRange()[iVar] ;
        sensitivityTable << setw(15) << var.get_varname() ; 
      }
      sensitivityTable << endl ;

      // Loop over perturbed parameter values.

      long long seed(0) ;
      for (size_t itr(1) ; itr <= m_maxIterations ; itr++) {
        vector<double> rndmd(m_nVar,0.0) ;
        sobol.sobol(rndmd.size(), &seed, rndmd) ;

        // Use random vector generated by sobol method to perturb parameter values.

        for (size_t j(0) ; j < rndmd.size() ; j++) {
          newLocation[j] = currentLocation[j]*(1.0 + (m_delta[j])*(rndmd[j] - 0.5)) ;
        }

        // Set perturbed parameters and calculate new quantities.

        SetLocation(newLocation) ;
        vector<double> Quantities ;
        pSys->calculate(nCnd, Quantities, false) ;

        // Write perturbed parameters and values calculated from them.
        // SHR: Note the loop over the calculated values is two so as
        // to miss out the expt. values.

        for (size_t j(0) ; j < newLocation.size() ; j++) {
          sensitivityTable << formatFloat(rndmd[j], 5, 15) ;
        }
        for (size_t j(1) ; j < Quantities.size() ; j += 2) {
          sensitivityTable << formatFloat(Quantities[j], 5, 15) ;
        }
        sensitivityTable << endl ;

      }

      cinfo << sensitivityTable.str() << endl ;
    }

    return true;

  }

  bool SensitivityAnalysis::DoCalculationNew(System* pSys)
  {

    m_nVar = Rdouble::withRange().size() ;

    if (m_nVar < 2) { 
      cerr << "Sensitivity analysis requries at least two range variables to be set." << endl;
      return false ;
    }

    //Read variable uncertainties from range
    for (size_t iVar(0) ; iVar < m_nVar ; iVar++) {
      Rdouble var = *Rdouble::withRange()[iVar] ;
      double lower = var.get_lower();
      double upper = var.get_upper();
      m_delta.push_back(abs((upper - lower)/2.0 ));
    }

    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are calculated).
    pSys->m_Flags.useTheSameCellNumber = true;

    // Uncomment to enable ctest output during fitting. Or use -w5 option in command.
    //ChangeErrorLevel e(obDebug); 

    //Default is to disable ctest during fitting. Restored when leaving this function.
    StopCTestOutput stop(true) ;

    //
    // Begin by finding the starting point chi-squared value.
    //

    vector<double> originalLocation(m_nVar,0.0) ; 
    vector<double> newLocation(m_nVar,0.0) ; 

    GetLocation(originalLocation) ;

    // Invoke SetLocation to catch any constrained parameters.
    SetLocation(originalLocation) ;

    vector<double> Temperature ;
    vector<double> Concentration ;
    pSys->getConditionsManager()->getConditions (Temperature, Concentration) ;

    // Instantiate a random vector generator.
    Sobol sobol ;
    long long seed(0) ;

    // Loop over conditions. 
    size_t nConditions = Temperature.size() ;
    for (size_t nCnd(0) ; nCnd < nConditions ; nCnd++) {

      m_Di.clear() ;
      m_Dij.clear() ;

      vector<double> f0, f02, varf ;
      vector<double> alpha, beta ;
      vector<string> rxnId ;
      bool f0Initialized(false) ;

      // Loop over perturbed parameter values.

      double nrmlFctr(1.0/double(m_maxIterations)) ;
      for (size_t itr(1) ; itr <= m_maxIterations ; itr++) {
        vector<double> rndmd(m_nVar,0.0) ;
        sobol.sobol(m_nVar, &seed, rndmd) ;

        // Use random vector generated by sobol method to perturb parameter values.

        for (size_t j(0) ; j < rndmd.size() ; j++) {
          newLocation[j] = originalLocation[j] + m_delta[j]*(rndmd[j] - 0.5) ;
        }

        // Set perturbed parameters and calculate new quantities.

        SetLocation(newLocation) ;
        map<string, double> phenRates ;
        pSys->calculate(nCnd, phenRates);

        if (!f0Initialized) {
          m_nOut = phenRates.size() ;
          f0.resize(m_nOut, 0.0) ;
          f02.resize(m_nOut, 0.0) ;
          varf.resize(m_nOut, 0.0) ;
          alpha.resize(m_nOut*m_nVar*m_order, 0.0) ;
          beta.resize(m_nOut*m_nVar*m_nVar*m_order*m_order, 0.0) ;
          map<string, double>::const_iterator irxn = phenRates.begin();
          for (; irxn != phenRates.end(); irxn++) {
            rxnId.push_back(irxn->first) ;
          }
          f0Initialized = true ;
        }

        // Calculate alpha values.		
        map<string, double>::const_iterator irxn = phenRates.begin();
        for (size_t nOut(0), ida(0), idb(0) ; irxn != phenRates.end() ; irxn++, nOut++) {

          double output = irxn->second ;
          f0[nOut]  += output ;
          f02[nOut] += output*output ;

          // Calculate all alpha and beta coefficients for the up to m_order polynomials.

          for (size_t i(0) ; i < rndmd.size() ; i++) {

            // alpha coefficients.

            double input_i = rndmd[i] ;
            for (size_t k(1) ; k <= m_order ; k++, ida++) {
              alpha[ida] += output * ShiftedLegendre(k, input_i) ; ;
            }

            // beta coefficients.

            for (size_t j(0) ; j < rndmd.size() ; j++) {
              double input_j = rndmd[j] ;
              for (size_t k(1) ; k <= m_order ; k++) {
                for (size_t l(1) ; l <= m_order ; l++, idb++) {
                  beta[idb] += output * ShiftedLegendre(k, input_i) * ShiftedLegendre(l, input_j); 
                }
              }
            }

          }

        } // End of Ouputs loop.

      } // End of Iteration loop.

      for (size_t i(0) ; i < alpha.size() ; i++) 
        alpha[i] *= nrmlFctr ;

      for (size_t i(0) ; i < beta.size() ; i++) 
        beta[i] *= nrmlFctr ;

      // Calculate the variance of the outputs.
      for (size_t i(0) ; i < f0.size() ; i++) {
        f0[i]  *= nrmlFctr ;
        f02[i] *= nrmlFctr ;
        varf[i] = f02[i] - f0[i]*f0[i] ;
      }

      sensitivityIndicies(varf, alpha, beta) ;

      WriteOutAnalysis(rxnId, f0, alpha, Temperature[nCnd], Concentration[nCnd]) ;

    } // End of conditions loop. 

    return true;
  }

  // This method calculates sensitivity indicies.
  bool SensitivityAnalysis::sensitivityIndicies(const vector<double> &varf, const vector<double> &alpha, const vector<double> &beta) {

    for (size_t nOut(0), ida(0), idb(0) ; nOut < m_nOut ; nOut++) {
      double rvar(1.0/varf[nOut]) ;
      vector<double> tma ;
      for (size_t i(0) ; i < m_nVar ; i++) {

        // First order indicies.
        double sma(0.0) ;
        for (size_t k(1) ; k <= m_order ; k++, ida++) {
          sma += alpha[ida]*alpha[ida] ;
        }
        tma.push_back(sma*rvar) ;

        // Second order indicies.
        vector<double> tmb ;
        for (size_t j(0) ; j < m_nVar ; j++) {
          double smb(0.0) ;
          for (size_t k(1) ; k <= m_order ; k++) {
            for (size_t l(1) ; l <= m_order ; l++, idb++) {
              smb += beta[idb]*beta[idb];
            }
          }
          tmb.push_back(smb*rvar) ;
        }
        m_Dij.push_back(tmb) ;
      }
      m_Di.push_back(tma) ;
    }

    return true ;
  }

  // This method generates data for external analysis.
  bool SensitivityAnalysis::WriteOutAnalysis(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &alpha, double Temperature, double Concentration) {

    cinfo << endl << "Sensitivity Analysis: Temperature = " << setw(5) << setprecision(4) << Temperature << ", Concentration = " << Concentration << endl << endl;

    cinfo << "Input variable key:"  << endl << endl;
    for (size_t iVar(0), idx(1) ; iVar < m_nVar ; iVar++, idx++) {
      Rdouble var = *Rdouble::withRange()[iVar] ;
      cinfo << "  " << setw(30) << left << var.get_varname() << ": (" << idx << ")" << endl ;
    }
    cinfo << endl;

    for (size_t i(0), idx(0) ; i < m_nOut ; i++) {

      // First order sensistivity indices.

      vector<double> &tmp = m_Di[i] ;
      cinfo << " First order indices for " << rxnId[i] << ": " << endl ;
      cinfo << columnHeader() << endl;
      cinfo << "    " ;
      for (size_t j(0) ; j < m_nVar ; j++) {
        cinfo << formatFloat(tmp[j], 5, 15) ;
      }
      cinfo << endl << endl;

      // Second order sensistivity indicies.

      cinfo << " Second order indices for "  << rxnId[i] << ": "  << endl ;
      cinfo << columnHeader() << endl;
      for (size_t j(0) ; j < m_nVar ; j++, idx++) {
        stringstream ss ;
        ss << "(" << j+1 << ")" ;
        cinfo << setw(4)  << left << ss.str() ;
        vector<double> &tmp = m_Dij[idx] ;
        for (size_t k(0) ; k < j ; k++) {
          cinfo << formatFloat(tmp[k], 5, 15) ;
        }
        for (size_t k(j) ; k < m_nVar ; k++) {
          cinfo << setw(15) << right << "--" ;
        }
        cinfo << endl;
      }
      cinfo << endl;
    }
    cinfo << endl;

    return true ;
  }

  // This method generates the column header for the output tables.
  string SensitivityAnalysis::columnHeader() const {
    stringstream ss ;
    ss << "    " ;
    for (size_t j(1) ; j <= m_nVar ; j++) {
      stringstream tmp ;
      tmp <<  "(" << j << ")" ;
      ss << setw(15) << tmp.str() ;
    }	  
    return ss.str() ;
  }

} //namespace

