//-------------------------------------------------------------------------------------------
//
// SensitivityAnalysis.cpp
//
// Author: Struan Robertson
// Date:   25/Nov/2012
//
// This class implements sensitivity analysis algorithm used to identify this most important
// parameters from a fit. It is based on the Li et al, Chem. Eng. Sci. Vol. 57, 4445 (2002)
// and guided by the matlab implementation of Tilo Ziehn.
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"
#include "FittingUtils.h"
#include "../Sobol.h"
#include <string>
#include <cctype>

namespace {

  void ToUpper(std::string &str) {
    for (size_t i(0); i < str.size(); i++) {
      char c = std::toupper(str[i]);
      str[i] = c;
    }
  }

  // The first 15 shifted Legendre polynomials

  double Legendre_1(const double &x) { return sqrt(3.0)*(2.0*x - 1.0); }

  double Legendre_2(const double &x) { return sqrt(5.0)*(6.0*x*x - 6.0*x + 1.0); }

  double Legendre_3(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    return sqrt(7.)*(20.0*x3 - 30.0*x2 + 12.0*x - 1.);
  }

  double Legendre_4(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    return 3.0*(70.*x4 - 140.*x3 + 90.0*x2 - 20.0*x + 1.);
  }

  double Legendre_5(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    return 630.*sqrt(11.)*((2. / 5.)*x5 - x4 + (8. / 9.)*x3 - (1. / 3.)*x2 + (1. / 21.)*x - (1. / 630.));
  }

  double Legendre_6(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    return 3150.*sqrt(13.)*((22. / 75.)*x6 - (22. / 25.)*x5 + x4 - (8. / 15.)*x3 + (2. / 15.)*x2 - (1. / 75.)*x + (1. / 3150.));
  }

  double Legendre_7(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    return 16632.*sqrt(15.)*((13. / 63.)*x7 - (13. / 18.)*x6 + x5 - (25. / 36.)*x4 + (25. / 99.)*x3 - (1. / 22.)*x2 + (1. / 297.)*x - (1. / 16632.));
  }

  double Legendre_8(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    return 84084.*sqrt(17.)*((15. / 98.)*x8 - (30. / 49.)*x7 + x6 - (6. / 7.)*x5 + (75. / 182.)*x4
      - (10. / 91.)*x3 + (15. / 1001.)*x2 - (6. / 7007.)*x + (1. / 84084.));
  }

  double Legendre_9(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    return 420420.*sqrt(19.)*((17. / 147.)*x9 - (51. / 98.)*x8 + (48. / 49.)*x7 - x6 + (3. / 5.)*x5 - (3. / 14.)*x4
      + (4. / 91.)*x3 - (3. / 637.)*x2 + (3. / 14014.)*x - (1. / 420420.));
  }

  double Legendre_10(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    double x10 = x9*x;
    return sqrt(21.)*(184756.*x10 - 923780.*x9 + 1969110.*x8 - 2333760.*x7 +
      1681680.*x6 - 756756.*x5 + 210210.*x4 - 34320.*x3 + 2970.*x2 - 110.*x + 1.);
  }

  double Legendre_11(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    double x10 = x9*x;
    double x11 = x10*x;
    return sqrt(23.)*(705432.*x11 - 3879876.*x10 + 9237800.*x9 - 12471030.*x8 + 10501920.*x7 - 5717712.*x6 +
      2018016.*x5 - 450450.*x4 + 60060.*x3 - 4290.*x2 + 132.*x - 1.);
  }

  double Legendre_12(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    double x10 = x9*x;
    double x11 = x10*x;
    double x12 = x11*x;
    return 5.*(2704156.*x12 - 16224936.*x11 + 42678636.*x10 - 64664600.*x9 + 62355150.*x8 - 39907296.*x7 +
      17153136.*x6 - 4900896.*x5 + 900900.*x4 - 100100.*x3 + 6006.*x2 - 156.*x + 1.);
  }

  double Legendre_13(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    double x10 = x9*x;
    double x11 = x10*x;
    double x12 = x11*x;
    double x13 = x12*x;
    return sqrt(27.)*(10400600.*x13 - 67603900.*x12 + 194699232.*x11 - 327202876.*x10 + 355655300.*x9 - 261891630.*x8 +
      133024320.*x7 - 46558512.*x6 + 11027016.*x5 - 1701700.*x4 + 160160.*x3 - 8190.*x2 + 182.*x - 1.);
  }

  double Legendre_14(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    double x10 = x9*x;
    double x11 = x10*x;
    double x12 = x11*x;
    double x13 = x12*x;
    double x14 = x13*x;
    return sqrt(29.)*(40116600.*x14 - 280816200.*x13 + 878850700.*x12 - 1622493600.*x11 + 1963217256.*x10 - 1636014380.*x9 +
      960269310.*x8 - 399072960.*x7 + 116396280.*x6 - 23279256.*x5 + 3063060.*x4 - 247520.*x3 + 10920.*x2 - 210.*x + 1.);
  }

  double Legendre_15(const double &x) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double x6 = x5*x;
    double x7 = x6*x;
    double x8 = x7*x;
    double x9 = x8*x;
    double x10 = x9*x;
    double x11 = x10*x;
    double x12 = x11*x;
    double x13 = x12*x;
    double x14 = x13*x;
    double x15 = x14*x;
    return sqrt(32.)*(155117520.*x15 - 1163381400.*x14 + 3931426800.*x13 - 7909656300.*x12 + 10546208400.*x11 - 9816086280.*x10 +
      6544057520.*x9 - 3155170590.*x8 + 1097450640.*x7 - 271591320.*x6 + 46558512.*x5 - 5290740.*x4 + 371280.*x3 - 14280.*x2 + 240.*x - 1.);
  }

  // Some utility classes used to generate probability densities and calculate related 
  // orthogonal functions.

  // Abstract base class.

  class probDensity {
  public:
    typedef double(*ofn)(const double &x);
    probDensity() {
      m_slpMap[1] = Legendre_1;
      m_slpMap[2] = Legendre_2;
      m_slpMap[3] = Legendre_3;
      m_slpMap[4] = Legendre_4;
      m_slpMap[5] = Legendre_5;
      m_slpMap[6] = Legendre_6;
      m_slpMap[7] = Legendre_7;
      m_slpMap[8] = Legendre_8;
      m_slpMap[9] = Legendre_9;
      m_slpMap[10] = Legendre_10;
      m_slpMap[11] = Legendre_11;
      m_slpMap[12] = Legendre_12;
      m_slpMap[13] = Legendre_13;
      m_slpMap[14] = Legendre_14;
      m_slpMap[15] = Legendre_15;
    };
    virtual ~probDensity() {};

    // Read data need to initialise distribution.
    virtual bool initializeDist(System* pSys, string &err) = 0;

    // Determine new location from random variate.
    virtual void rndLocation(const vector<double> &rndmd, const vector<double> &currentLoc, vector<double> &newLoc) = 0;

    // Calculate orthogonal function of specifed order.
    virtual double orthFn(size_t order, const double x) {
      double f = m_slpMap[order](x);
      return f;
    };

  private:
    map<size_t, ofn> m_slpMap;
  };

  // Class for calculating a shifted uniform desnity and shifted Legendre polynomials.

  class ShiftedLegendre : public probDensity {
  public:
    ShiftedLegendre() : probDensity() { };
    ~ShiftedLegendre() { };

    bool initializeDist(System* pSys, string &err) {
      // Read variable uncertainties from range.
      size_t m_nVar = Rdouble::withRange().size();
      for (size_t iVar(0); iVar < m_nVar; iVar++) {
        Rdouble var = *Rdouble::withRange()[iVar];
        double lower = var.get_lower();
        double upper = var.get_upper();
        m_delta.push_back(abs((upper - lower) / 2.0));
      }
      return true;
    }

    virtual void rndLocation(const vector<double> &rndmd, const vector<double> &currentLoc, vector<double> &newLoc) {
      for (size_t j(0); j < rndmd.size(); j++) {
        newLoc[j] = currentLoc[j] + m_delta[j] * (rndmd[j] - 0.5);
      }
    }

  private:
    vector<double>   m_delta;
  };

  // Class for generating a Gaussian density and Legendre functions.

  class CorrelatedLegendre : public probDensity {
  public:
    CorrelatedLegendre() : probDensity() { };
    ~CorrelatedLegendre() { };

    bool initializeDist(System* pSys, string &err) {
      PersistPtr ppa = pSys->getPersistPtr()->XmlMoveTo("me:analysis");
      PersistPtr pp = ppa->XmlMoveTo("me:covariance");
      if (!pp)
        pp = ppa->XmlMoveTo("me:hessian");
      pp = pp->XmlMoveTo("matrix");
      if ((m_CorrelMtx = ReadMatrix<double>(pp))) {
        m_CorrelMtx->cholesky();
      }
      else {
        err = "Correlation Matrix not found";
        return false;
      }
      return true;
    };

    void rndLocation(const vector<double> &rndmd, const vector<double> &currentLoc, vector<double> &newLoc) {
      size_t nVar = rndmd.size();
      // Take inverse cumulative distribution of each sobol element.
      for (size_t j(0); j < nVar; j++) {
        newLoc[j] = NormalCDFInverse(rndmd[j]);
      }
      // Multiply the InvNorm with the cholesky decompostion.
      newLoc *= (*m_CorrelMtx);
      for (size_t j(0); j < nVar; j++) {
        newLoc[j] += currentLoc[j];
      }
    }

  private:
    dMatrix *m_CorrelMtx;
  };

}

namespace mesmer
{
  class SensitivityAnalysis : public CalcMethod, private FittingUtils
  {
  public:

    SensitivityAnalysis(const char* id) : m_probDensity(NULL),
      m_id(id),
      m_pSA(),
      m_VarRedMthd(RATIOCONTROL),
      m_nVar(0),
      m_nOut(0),
      m_nVred(5),
      m_nSample(0),
      m_bGenerateData(true),
      m_bCorrelatedData(false),
      m_order(1),
      m_Di(),
      m_Dij(),
      m_RSquared(),
      m_pParallelManager(NULL),
      m_rank(0),
      m_nRanks(1) {
      Register();
    }

    virtual ~SensitivityAnalysis() { delete m_probDensity; }
    virtual const char* getID() { return m_id; }
    virtual SensitivityAnalysis* Clone() { return new SensitivityAnalysis(*this); }
    virtual bool ParseData(PersistPtr pp);

    // Function to do the work.
    virtual bool DoCalculation(System* pSys);

  private:

    typedef vector<vector<double> > Table;
    typedef map<string, double>::const_iterator RxnItr;

    enum VarRedMthd {
      RATIOCONTROL,
      ADDITIVECONTROL,
      UNDEFINED
    };

    // This method does a complete sensitivity analysis.
    bool DoCalculationNew(System* pSys);

    // This method generates data for external analysis.
    bool DoCalculationOld(System* pSys);

    // This method writes out the results of a sensitivity analysis to test file.
    bool WriteOutAnalysisToTest(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &varf, double Temperature, double Concentration);

    // This method writes out the results of a sensitivity analysis.
    bool WriteOutAnalysis(const vector<string> &rxnId, const vector<double> &varf, double Temperature, double Concentration);

    // This method calculates the legendre expansion coefficients by MC integration and
    // applies the ratio control variate variance reduction method to improve convergence.
    bool CalcExpansionCoeffsRatioCV(const Table &locationData, const Table &outputData, Table &alphas, Table &betas);

    // This method calculates the legendre expansion coefficients by MC integration and
    // applies the additive control variate variance reduction method to improve convergence.
    bool CalcExpansionCoeffsAddCV(const Table &locationData, const Table &outputData, Table &alphas, Table &betas);

    // This method calculates sensitivity indicies.
    bool sensitivityIndicies(const vector<double> &f0, const Table &alphas, const Table &betas);

    // This method calculates R-squared values (measure of the quality of the fit).
    bool RSquared(const vector<double> &f0, vector<double> varf, Table &alphas, const Table &betas, const Table &locationData, const Table &outputData);

    // This method returns the legendre expansion estimate of the output value.
    double orthFnExpn(vector<double> &x, const double &f0, const vector<double> &alpha, const vector<double> &beta);

    // This method generates the column header for the output tables.
    string columnHeader() const;

    // This method writes the input variable key.
    void WriteInputVariableKey(stringstream &Key) const;

    // During parallel execution this method accumulates data across ranks.
    void AccumulateDataAcrossRanks(Table& table);

    // Set of orthogonal functions to be used in calculation.
    probDensity *m_probDensity;

    const char* m_id;

    PersistPtr m_pSA;

    VarRedMthd m_VarRedMthd;

    size_t m_nVar;             // Dimension of analysis - number of inputs.
    size_t m_nOut;             // Dimension of analysis - number of outputs.
    size_t m_nVred;            // Number of variance reduction iterations.

    size_t m_nSample;
    bool   m_bGenerateData;
    bool   m_bCorrelatedData;

    size_t m_order;            // The order of the HDMR analysis to use.

    Table m_Di;                // First order sensitivities.
    Table m_Dij;               // Second order sensitivities.

    vector<double> m_RSquared; // R^2 statistic.

    ParallelManager *m_pParallelManager; // Location of parallel manager.
    int m_rank;                          // MPI rank
    int m_nRanks;                        // MPI number of ranks

  };

  ////////////////////////////////////////////////
  //Global instance
  SensitivityAnalysis theSensitivityAnalysis("SensitivityAnalysis");
  ///////////////////////////////////////////////

  bool SensitivityAnalysis::ParseData(PersistPtr pp)
  {
    // Read in sensitivity analysis parameters, or use values from defaults.xml.
    m_nSample = pp->XmlReadInteger("me:sensitivityAnalysisSamples");
    m_bGenerateData = pp->XmlReadBoolean("me:sensitivityGenerateData");
    m_bCorrelatedData = pp->XmlReadBoolean("me:sensitivityCorrelatedData");
    size_t order = pp->XmlReadInteger("me:sensitivityAnalysisOrder", optional);
    if (!IsNan(order))
      m_order = order;
    size_t nVred = pp->XmlReadInteger("me:sensitivityNumVarRedIters", optional);
    if (!IsNan(nVred))
      m_nVred = nVred;
    const char* txt = pp->XmlReadValue("me:sensitivityVarRedMethod", optional); //hard-wired default is RATIOCONTROL
    if (txt) {
      string str(txt);
      ToUpper(str);
      m_VarRedMthd = (str == "RATIOCONTROL") ? RATIOCONTROL :
        (str == "ADDITIVECONTROL") ? ADDITIVECONTROL : UNDEFINED;

      if (m_bCorrelatedData) {
        m_probDensity = new CorrelatedLegendre();
      }
      else {
        m_probDensity = new ShiftedLegendre();
      }
    }

    // Store pointer for output.
    m_pSA = pp;

    return true;
  }

  bool SensitivityAnalysis::DoCalculation(System* pSys) {
      m_pParallelManager = pSys->getParallelManager();
      m_rank = m_pParallelManager->rank();
      m_nRanks = m_pParallelManager->size();
    return (m_bGenerateData) ? DoCalculationOld(pSys) : DoCalculationNew(pSys);
  }

  bool SensitivityAnalysis::DoCalculationOld(System* pSys)
  {
    m_nVar = Rdouble::withRange().size();

    if (m_nVar < 2) {
      cerr << "Sensitivity analysis requries at least two range variables to be set." << endl;
      return false;
    }

    string err("");
    if (!m_probDensity->initializeDist(pSys, err)) {
      cerr << err << endl;
      return false;
    }

    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are calculated).
    pSys->m_Flags.useTheSameCellNumber = true;

    //Default is to disable ctest during analysis. Restored when leaving this function.
    StopCTestOutput stop(true);

    vector<double> originalLocation(m_nVar, 0.0);
    vector<double> newLocation(m_nVar, 0.0);

    GetLocation(originalLocation);

    // Invoke SetLocation to catch any constrained parameters.
    SetLocation(originalLocation);

    vector<double> Temperature;
    vector<double> Concentration;
    pSys->getConditionsManager()->getConditions(Temperature, Concentration);

    // Instantiate a random vector generator.
    Sobol sobol;

    // Loop over condiditons. 
    size_t nConditions = Temperature.size();
    for (size_t nCnd(0); nCnd < nConditions; nCnd++) {

      // Test calculation to determine number of outputs.
      SetLocation(originalLocation);
      map<string, double> phenRates;
      pSys->calculate(nCnd, phenRates);

      m_nOut = phenRates.size();
      stringstream outputKey;
      outputKey << "Output variable key:" << endl << endl;
      size_t idx(1);
      for (RxnItr irxn = phenRates.begin(); irxn != phenRates.end(); irxn++, idx++) {
        outputKey << "  " << setw(30) << left << irxn->first << ": k(" << idx << ")" << endl;
      }
      stringstream inputKey;
      WriteInputVariableKey(inputKey);
      cinfo << inputKey.str();
      cinfo << outputKey.str();

      // String stream to hold results. 
      stringstream sensitivityTable;

      // Write out table header.

      sensitivityTable << endl;
      sensitivityTable << "Sensitivity Table" << endl;
      sensitivityTable << "  Temperature:   " << formatFloat(Temperature[nCnd], 5, 15) << " K" << endl;
      sensitivityTable << "  Concentration: " << formatFloat(Concentration[nCnd], 5, 15) << " ppcc" << endl;
      sensitivityTable << endl;
      for (size_t iVar(0), idx(1); iVar < m_nVar; iVar++, idx++) {
        stringstream strIdx;
        strIdx << "(" << idx << ")";
        sensitivityTable << setw(15) << strIdx.str();
      }
      for (size_t iOut(0), idx(1); iOut < m_nOut; iOut++, idx++) {
        stringstream strIdx;
        strIdx << "k(" << idx << ")";
        sensitivityTable << setw(15) << strIdx.str();
      }
      sensitivityTable << endl;

      // Loop over perturbed parameter values.

      long long seed(m_rank*m_nSample / m_nRanks);
      for (size_t itr(m_rank + 1); itr <= m_nSample; itr += m_nRanks) {
        vector<double> rndmd(m_nVar, 0.0);
        sobol.sobol(rndmd.size(), &seed, rndmd);

        // Use random vector generated by sobol method to perturb parameter values.

        m_probDensity->rndLocation(rndmd, originalLocation, newLocation);

        // Set perturbed parameters and calculate new quantities.

        SetLocation(newLocation);

        map<string, double> phenRates;

        // As some perturbations will produce unrealistic configurations the
    // diagonalizer will fail, in which case this configuration is skipped.
        try {
          pSys->calculate(nCnd, phenRates);
        }
        catch (...) {
          continue;
        }

        for (size_t j(0); j < newLocation.size(); j++) {
          sensitivityTable << formatFloat(rndmd[j], 5, 15);
        }
        for (RxnItr irxn = phenRates.begin(); irxn != phenRates.end(); irxn++) {
          sensitivityTable << formatFloat(irxn->second, 5, 15);
        }

        sensitivityTable << endl;

      }

      // Have each thread print out in order.
      for (size_t i(0); i < size_t(m_nRanks); i++) {
				if (size_t(m_nRanks) == i)
          cinfo << sensitivityTable.str() << endl;
        pSys->getParallelManager()->barrier();
      }

    }

    return true;
  }

  bool SensitivityAnalysis::DoCalculationNew(System* pSys)
  {
    m_nVar = Rdouble::withRange().size();

    if (m_nVar < 2) {
      cerr << "Sensitivity analysis requries at least two range variables to be set." << endl;
      return false;
    }

    // Prepare output.
    m_pSA = m_pSA->XmlWriteMainElement("me:sensitivityAnalysisTables", "", true); // Will replace an existing element.

    string err("");
    if (!m_probDensity->initializeDist(pSys, err)) {
      cerr << err << endl;
      return false;
    }

    // Do not output all the intermediate results to XML.
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are calculated).
    pSys->m_Flags.useTheSameCellNumber = true;

    vector<double> originalLocation(m_nVar, 0.0);
    vector<double> newLocation(m_nVar, 0.0);

    GetLocation(originalLocation);

    // Invoke SetLocation to catch any constrained parameters.
    SetLocation(originalLocation);

    vector<double> Temperature;
    vector<double> Concentration;
    pSys->getConditionsManager()->getConditions(Temperature, Concentration);

    // Instantiate a random vector generator.
    Sobol sobol;

    long long seed(m_rank*m_nSample / m_nRanks + 10);

    // Loop over conditions. 
    size_t nConditions = Temperature.size();
    for (size_t nCnd(0); nCnd < nConditions; nCnd++) {

      // Suppress output to .test.
      ctest.clear(std::ios::failbit);

      m_Di.clear();
      m_Dij.clear();

      // Test calculation to determine number of outputs.
      SetLocation(originalLocation);
      map<string, double> phenRates;
      pSys->calculate(nCnd, phenRates);

      m_nOut = phenRates.size();
      vector<double> f0(m_nOut, 0.0), varf(m_nOut, 0.0);
      vector<string> rxnId;
      for (RxnItr irxn = phenRates.begin(); irxn != phenRates.end(); irxn++) {
        rxnId.push_back(irxn->first);
      }
      Table locationData;
      Table outputData;

      // Loop over perturbed parameter values.

      size_t cnt(0);
      for (size_t itr(m_rank); itr < m_nSample; itr += m_nRanks) {

        vector<double> rndmd(m_nVar, 0.0);
        sobol.sobol(m_nVar, &seed, rndmd);

        // Use random vector generated by sobol method to perturb parameter values.

        m_probDensity->rndLocation(rndmd, originalLocation, newLocation);

        // Set perturbed parameters and calculate new quantities.

        SetLocation(newLocation);
        map<string, double> phenRates;

        // As some perturbations will produce unrealistic configurations the
        // diagonalizer will fail, in which case this configuration is skipped.
        try {
          pSys->calculate(nCnd, phenRates);
        }
        catch (...) {
          continue;
        }

				vector<double> trial;
        cnt++;
        RxnItr irxn = phenRates.begin();
        for (size_t nOut(0); irxn != phenRates.end(); irxn++, nOut++) {
          double output = irxn->second;
          trial.push_back(output);
          f0[nOut] += output;
          varf[nOut] += output*output;
        }
				locationData.push_back(rndmd);
				outputData.push_back(trial);
			}

      // Calculate mean and variance across ranks.
      double dCnt(cnt);
      m_pParallelManager->sumDouble(&dCnt, 1);
      m_pParallelManager->sumDouble(&f0[0], f0.size());
      m_pParallelManager->sumDouble(&varf[0], varf.size());
      for (size_t nOut(0); nOut < varf.size(); nOut++) {
        f0[nOut] /= dCnt;
        varf[nOut] -= f0[nOut] * f0[nOut];
        varf[nOut] /= (dCnt - 1);
      }
			AccumulateDataAcrossRanks(locationData);
			AccumulateDataAcrossRanks(outputData);

      // Calculate all alpha and beta coefficients for the up to m_order polynomials.

      Table alphas(m_nOut, vector<double>(m_nVar*m_order, 0.0));
      Table betas(m_nOut, vector<double>(m_nVar*(m_nVar - 1)*m_order*m_order / 2, 0.0));

      switch (m_VarRedMthd) {
      case RATIOCONTROL:
        CalcExpansionCoeffsRatioCV(locationData, outputData, alphas, betas);
        break;
      case ADDITIVECONTROL:
        CalcExpansionCoeffsAddCV(locationData, outputData, alphas, betas);
        break;
      default:
        break;
      }

      // Calculate the R^2 statistic. 
      RSquared(f0, varf, alphas, betas, locationData, outputData);

      // Calculate sensitivity indicies. 
      sensitivityIndicies(varf, alphas, betas);

      ctest.clear();

      // Write out results. 
      WriteOutAnalysisToTest(rxnId, f0, varf, Temperature[nCnd], Concentration[nCnd]);

      WriteOutAnalysis(rxnId, varf, Temperature[nCnd], Concentration[nCnd]);

    } // End of conditions loop. 

    return true;
  }

  // This method calculates the legendre expansion coefficients by MC integration and
  // applies the ratio control variate variance reduction method to improve convergence.
  bool SensitivityAnalysis::CalcExpansionCoeffsRatioCV(const Table &locationData, const Table &outputData, Table &alphas, Table &betas) {

    for (size_t nOut(0); nOut < m_nOut; nOut++) {

      vector<double> &alpha = alphas[nOut];
      vector<double> &beta = betas[nOut];
      double f0(0.0);

      for (size_t itr(0); itr < locationData.size(); itr++) {
        vector<double> rndmd = locationData[itr];
        double output = outputData[itr][nOut];

        f0 += output;

        for (size_t i(0), ida(0), idb(0); i < rndmd.size(); i++) {

          // alpha coefficients.

          double input_i = rndmd[i];
          for (size_t k(1); k <= m_order; k++, ida++) {
            alpha[ida] += output * m_probDensity->orthFn(k, input_i);
          }

          // beta coefficients.

          for (size_t j(0); j < i; j++) {
            double input_j = rndmd[j];
            for (size_t k(1); k <= m_order; k++) {
              for (size_t l(1); l <= m_order; l++, idb++) {
                beta[idb] += output * m_probDensity->orthFn(k, input_i)
                  * m_probDensity->orthFn(l, input_j);
              }
            }
          }

        }

      } // End of Iteration loop.

      // Normalize expansion coefficients.

      double nrmlFctr(1.0 / double(locationData.size()));

      f0 *= nrmlFctr;

      for (size_t j(0); j < alpha.size(); j++) {
        alpha[j] *= nrmlFctr;
      }

      for (size_t j(0); j < beta.size(); j++) {
        beta[j] *= nrmlFctr;
      }

      // Now iterate control variate ratio.

      vector<double> alpha0(alpha);
      vector<double> beta0(beta);
      for (size_t iVred(0); iVred < m_nVred; iVred++) {

        vector<double> alpha1(alpha.size(), 0.0);
        vector<double> beta1(beta.size(), 0.0);
        for (size_t itr(0); itr < locationData.size(); itr++) {
          vector<double> rndmd = locationData[itr];
          double h0 = orthFnExpn(rndmd, f0, alpha0, beta0);

          for (size_t i(0), ida(0), idb(0); i < rndmd.size(); i++) {

            // alpha coefficients.

            double input_i = rndmd[i];
            for (size_t k(1); k <= m_order; k++, ida++) {
              alpha1[ida] += h0 * m_probDensity->orthFn(k, input_i);
            }

            // beta coefficients.

            for (size_t j(0); j < i; j++) {
              double input_j = rndmd[j];
              for (size_t k(1); k <= m_order; k++) {
                for (size_t l(1); l <= m_order; l++, idb++) {
                  beta1[idb] += h0 * m_probDensity->orthFn(k, input_i)
                    * m_probDensity->orthFn(l, input_j);
                }
              }
            }
          }

        } // End of sample loop.

        for (size_t j(0); j < alpha.size(); j++) {
          alpha0[j] *= alpha[j] / (alpha1[j] * nrmlFctr);
        }

        for (size_t j(0); j < beta.size(); j++) {
          beta0[j] *= beta[j] / (beta1[j] * nrmlFctr);
        }
      } // End of iteration loop.

      alpha = alpha0;
      beta = beta0;

    } // End of Outputs loop.

    return true;
  }

  // This method calculates the legendre expansion coefficients by MC integration and
  // applies the additive control variate variance reduction method to improve convergence.
  bool SensitivityAnalysis::CalcExpansionCoeffsAddCV(const Table &locationData, const Table &outputData, Table &alphas, Table &betas) {

    for (size_t nOut(0); nOut < m_nOut; nOut++) {

      vector<double> &alpha = alphas[nOut];
      vector<double> &beta = betas[nOut];
      double f0(0.0);

      for (size_t iVred(0); iVred < m_nVred; iVred++) {

        vector<double> alpha0(alpha.size(), 0.0);
        vector<double> beta0(beta.size(), 0.0);
        double sum(0.0);
        for (size_t itr(0); itr < locationData.size(); itr++) {
          vector<double> rndmd = locationData[itr];

          double output = outputData[itr][nOut];
          sum += output;
          output -= orthFnExpn(rndmd, f0, alpha, beta);

          for (size_t i(0), ida(0), idb(0); i < rndmd.size(); i++) {

            // alpha coefficients.

            double input_i = rndmd[i];
            for (size_t k(1); k <= m_order; k++, ida++) {
              alpha0[ida] += output * m_probDensity->orthFn(k, input_i);
            }

            // beta coefficients.

            for (size_t j(0); j < i; j++) {
              double input_j = rndmd[j];
              for (size_t k(1); k <= m_order; k++) {
                for (size_t l(1); l <= m_order; l++, idb++) {
                  beta0[idb] += output * m_probDensity->orthFn(k, input_i)
                    * m_probDensity->orthFn(l, input_j);
                }
              }
            }

          }

        } // End of Iteration loop.

        // Normalize expansion coefficients.

        double nrmlFctr(1.0 / double(locationData.size()));

        f0 = nrmlFctr*sum;

        for (size_t j(0); j < alpha0.size(); j++) {
          alpha[j] += nrmlFctr*alpha0[j];
        }

        for (size_t j(0); j < beta0.size(); j++) {
          beta[j] += nrmlFctr*beta0[j];
        }
      }

    } // End of Outputs loop.

    return true;
  }

  // This method returns the legendre expansion estimate of the output value.
  double SensitivityAnalysis::orthFnExpn(vector<double> &x, const double &f0, const vector<double> &alpha, const vector<double> &beta) {

    // Initialize sum with the mean of the function.
    double sum(f0);

    for (size_t i(0), ida(0), idb(0); i < m_nVar; i++) {

      // alpha coefficients.

      double input_i = x[i];
      for (size_t k(1); k <= m_order; k++, ida++) {
        sum += alpha[ida] * m_probDensity->orthFn(k, input_i);
      }

      // beta coefficients.

      for (size_t j(0); j < i; j++) {
        double input_j = x[j];
        for (size_t k(1); k <= m_order; k++) {
          for (size_t l(1); l <= m_order; l++, idb++) {
            sum += beta[idb] * m_probDensity->orthFn(k, input_i)
              * m_probDensity->orthFn(l, input_j);
          }
        }
      }
    }

    return sum;
  }

  // This method calculates R-squared values (measure of the quality of the fit).
  bool SensitivityAnalysis::RSquared(const vector<double> &f0, vector<double> varf, Table &alphas, const Table &betas, const Table &locationData, const Table &outputData) {

    m_RSquared.resize(m_nOut, 0.0);
    vector<double> sumOfRes(m_nOut, 0.0);
    for (size_t nOut(0); nOut < f0.size(); nOut++) {

      for (size_t itr(0); itr < locationData.size(); itr++) {
        vector<double> rndmd = locationData[itr];
        double output = outputData[itr][nOut];

        // Calculate fitted value.

        double h0 = orthFnExpn(rndmd, f0[nOut], alphas[nOut], betas[nOut]);

        double residue = output - h0;
        sumOfRes[nOut] += residue*residue;
      }
    }

    for (size_t i(0); i < f0.size(); i++) {
      sumOfRes[i] /= double(locationData.size());
      m_RSquared[i] = 1.0 - (sumOfRes[i] / varf[i]);
    }

    return true;
  }

  // This method calculates sensitivity indicies.
  bool SensitivityAnalysis::sensitivityIndicies(const vector<double> &varf, const Table &alphas, const Table &betas) {

    for (size_t nOut(0); nOut < m_nOut; nOut++) {
      double rvar(1.0 / varf[nOut]);
      const vector<double> &alpha = alphas[nOut];
      const vector<double> &beta = betas[nOut];
      vector<double> tma;
      for (size_t i(0), ida(0), idb(0); i < m_nVar; i++) {

        // First order indicies.
        double sma(0.0);
        for (size_t k(1); k <= m_order; k++, ida++) {
          sma += alpha[ida] * alpha[ida];
        }
        tma.push_back(sma*rvar);

        // Second order indicies.
        vector<double> tmb;
        for (size_t j(0); j < i; j++) {
          double smb(0.0);
          for (size_t k(1); k <= m_order; k++) {
            for (size_t l(1); l <= m_order; l++, idb++) {
              smb += beta[idb] * beta[idb];
            }
          }
          tmb.push_back(smb*rvar);
        }
        m_Dij.push_back(tmb);
      }
      m_Di.push_back(tma);
    }

    return true;
  }

  // This method writes out the results of a sensitivity analysis.
  bool SensitivityAnalysis::WriteOutAnalysis(const vector<string> &rxnId, const vector<double> &varf, double Temperature, double Concentration) {

    // Begin table.

    PersistPtr pp = m_pSA->XmlWriteElement("me:sensitivityAnalysisTable");

    // Write out conditions

    pp->XmlWriteAttribute("temperature", Temperature, 2, true);
    pp->XmlWriteAttribute("concentration", Concentration, 2, false);

    for (size_t i(0), idx(0); i < m_nOut; i++) {

      PersistPtr ppSensInd = pp->XmlWriteElement("me:sensitivityIndices");
      ppSensInd->XmlWriteAttribute("reaction", rxnId[i]);

      // Write out first order sensitivities.

      vector<double> &FrtOdr = m_Di[i];
      for (size_t j(0); j < FrtOdr.size(); j++) {
        stringstream ss;
        ss << FrtOdr[j];
        PersistPtr pp1stOrdVal = ppSensInd->XmlWriteValueElement("me:firstOrderIndex", ss.str());
        Rdouble var = *Rdouble::withRange()[j];
        pp1stOrdVal->XmlWriteAttribute("key", var.get_varname());
      }

      // Write out second order sensitivities.

      for (size_t j(0); j < m_nVar; j++, idx++) {
        Rdouble var1 = *Rdouble::withRange()[j];
        vector<double> &SndOdr = m_Dij[idx];
        for (size_t k(0); k < j; k++) {
          stringstream ss;
          ss << SndOdr[k];
          PersistPtr pp1stOrdVal = ppSensInd->XmlWriteValueElement("me:secondOrderIndex", ss.str());
          Rdouble var2 = *Rdouble::withRange()[k];
          pp1stOrdVal->XmlWriteAttribute("key1", var1.get_varname());
          pp1stOrdVal->XmlWriteAttribute("key2", var2.get_varname());
        }
      }

      // Write out R^2 statistic and Standard deviation.

      stringstream ss;
      ss << m_RSquared[i];
      ppSensInd->XmlWriteValueElement("me:R_Squared", ss.str());
      ss.str("");
      ss << sqrt(varf[i]);
      ppSensInd->XmlWriteValueElement("me:standardDeviation", ss.str());
    }

    return true;
  }

  // This method generates data for external analysis.
  bool SensitivityAnalysis::WriteOutAnalysisToTest(const vector<string> &rxnId, const vector<double> &f0, const vector<double> &varf, double Temperature, double Concentration) {

    ctest << endl << "Sensitivity Analysis: Temperature = " << setw(5) << setprecision(4) << Temperature << ", Concentration = " << Concentration << endl << endl;

    stringstream Key;
    WriteInputVariableKey(Key);
    ctest << Key.str() << endl;

    for (size_t i(0), idx(0); i < m_nOut; i++) {

      ctest << " Standard deviation for " << rxnId[i] << ": " << sqrt(varf[i]) << endl << endl;

      // First order sensistivity indices.

      vector<double> &tmp = m_Di[i];
      ctest << " First order indices for " << rxnId[i] << ": " << endl;
      ctest << columnHeader() << endl;
      ctest << "    ";
      for (size_t j(0); j < m_nVar; j++) {
        ctest << formatFloat(tmp[j], 5, 15);
      }
      ctest << endl << endl;

      // Second order sensistivity indicies.

      ctest << " Second order indices for " << rxnId[i] << ": " << endl;
      ctest << columnHeader() << endl;
      for (size_t j(0); j < m_nVar; j++, idx++) {
        stringstream ss;
        ss << "(" << j + 1 << ")";
        ctest << setw(4) << left << ss.str();
        vector<double> &tmp = m_Dij[idx];
        for (size_t k(0); k < j; k++) {
          ctest << formatFloat(tmp[k], 5, 15);
        }
        for (size_t k(j); k < m_nVar; k++) {
          ctest << setw(15) << right << "--";
        }
        ctest << endl;
      }
      ctest << endl;
      ctest << " R-Squared statistic: " << m_RSquared[i] << endl << endl;
    }

    return true;
  }

  // This method generates the column header for the output tables.
  string SensitivityAnalysis::columnHeader() const {
    stringstream ss;
    ss << "    ";
    for (size_t j(1); j <= m_nVar; j++) {
      stringstream tmp;
      tmp << "(" << j << ")";
      ss << setw(15) << tmp.str();
    }
    return ss.str();
  }

  // This method writes the input variable key.
  void SensitivityAnalysis::WriteInputVariableKey(stringstream &Key) const {
    Key << "Input variable key:" << endl << endl;
    for (size_t iVar(0), idx(1); iVar < m_nVar; iVar++, idx++) {
      Rdouble var = *Rdouble::withRange()[iVar];
      Key << "  " << setw(30) << left << var.get_varname() << ": (" << idx << ")" << endl;
    }
    Key << endl;
  }

  // During parallel execution this method accumulates data across ranks.
  void SensitivityAnalysis::AccumulateDataAcrossRanks(Table& table) {
    // If there is only one rank return as there is no work to do.
    if (m_nRanks == 1)
      return;

		m_pParallelManager->barrier();

		Table wrk; // Workspace for the accumulated table.
		vector<double> tmp = table[0];
    for (size_t i(0); i < size_t(m_nRanks); i++) {
      int tableSize = int(table.size()) ;
			m_pParallelManager->broadcastInteger(&tableSize, 1, i);
			for (size_t j(0); j < size_t(tableSize); j++) {
				if (size_t(m_nRanks) == i)
          tmp = table[j];
        m_pParallelManager->broadcastDouble(&tmp[0], tmp.size(), i);
        wrk.push_back(tmp);
      }
    }
    table = wrk;
  }
} //namespace
