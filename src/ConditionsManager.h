#ifndef GUARD_ConditionsManager_h
#define GUARD_ConditionsManager_h

#include <vector>
#include "System.h"
#include "ParallelManager.h"
#include "AnalysisData.h"

namespace mesmer
{

  struct conditionSet
  {
  public:
    conditionSet(std::string ref1, std::string ref2, std::string refReaction, double value, double error) :
      m_ref1(ref1), m_ref2(ref2), m_refReaction(refReaction), m_value(value), m_error(error), m_calcValue(0.0)
    {}

    void get_conditionSet(std::string& ref1, std::string& ref2, std::string& refReaction, double& value, double& error) {
      ref1 = m_ref1;
      ref2 = m_ref2;
      refReaction = m_refReaction;
      value = m_value;
      error = m_error;
    }

    std::string m_ref1;
    std::string m_ref2;
    std::string m_refReaction;
    double m_value;
    double m_error;
    double m_calcValue;
  };

  struct RawDataSet
  {
    RawDataSet() :
      data(),
      m_calcTrace(),
      m_traceSize(0),
      m_ref1(""),
      m_Name(NULL),
      m_weight(0.0),
      m_pPersistPtr(NULL) {}
  public:
    std::vector<std::pair<double, double> > data; // Time series data.
    std::vector<double> m_traceErrors;            // Trace errors. 
    std::vector<double> m_calcTrace;              // Calculated trace. 
    size_t m_traceSize;                           // Number of points in trace (including baseline noise).
    std::string m_ref1;                           // Monitored species.
    const char* m_Name;                           // Name of trace.
    mutable double m_weight;                      // Weight data should be given - as determined by the experimenter.
    PersistPtr m_pPersistPtr;                     // Location of data.
  };

  // Functor to compare times in RawDataSet.
  class Before
  {
  private:
    double startT;
  public:
    Before(double start) :startT(start) {}
    bool operator()(std::pair<double, double> val) { return val.first < startT; }
  };

  // To make sure if there is a concentration or pressure definition, there is a temperature definition.
  struct CandTpair {

    double      m_concentration; // particles per cubic centimeter
    double      m_temperature;   // Kelvin
    double      m_radiationTemp; // Radiation temperature in Kelvin.
    Precision   m_precision;
    std::string m_pBathGasName;
    std::map<Reaction*, double> m_excessConcs; // Reaction, excess species conc in ppcc.
    std::map<Reaction*, bool>   m_percentExcessConc; // Reaction, excess species conc is a % of bath gas.

    std::vector<conditionSet> m_rates;
    std::vector<conditionSet> m_yields;
    std::vector<conditionSet> m_eigenvalues;
    std::vector<RawDataSet>   m_rawDataSets;
    std::vector<PersistPtr>   m_expDataPtrs;

    std::string m_group; // Used when comparing data from different Labs. 

    AnalysisData m_analysisData;

    CandTpair(double cp_, double t_, Precision _pre, const char* _bathGas,
      const map<Reaction*, double>& _excessConcs, const char* _group, double _radiationTemp)
      : m_concentration(cp_), m_temperature(t_), m_precision(_pre),
      m_pBathGasName(_bathGas), m_excessConcs(_excessConcs), m_group(_group), m_radiationTemp(_radiationTemp)  {}

    void set_experimentalRates(PersistPtr ppData, std::string ref1, std::string ref2, std::string refReaction, double value, double error) {
      if (ref1.size() > 0 && ref2.size() > 0) {
        m_rates.push_back(conditionSet(ref1, ref2, refReaction, value, error));
        m_expDataPtrs.push_back(ppData);
      }
    }

    void set_experimentalYields(PersistPtr ppData, std::string ref, std::string yieldTime, double value, double error) {
      if (ref.size() > 0 && yieldTime.size() > 0) {
        m_yields.push_back(conditionSet(ref, yieldTime, std::string(""), value, error));
        m_expDataPtrs.push_back(ppData);
      }
    }

    void set_experimentalEigenvalues(PersistPtr ppData, std::string eigenvalueID, double value, double error) {
      m_eigenvalues.push_back(conditionSet(std::string(""), std::string(""), eigenvalueID, value, error));
      m_expDataPtrs.push_back(ppData);
    }

    void set_calcRates(std::vector<double>& calcRates) {
      for (size_t i(0); i < calcRates.size(); ++i) {
        m_rates[i].m_calcValue = calcRates[i];
      }
    };

    void set_calcYields(std::vector<double>& calcYields) {
      for (size_t i(0); i < calcYields.size(); ++i) {
        m_yields[i].m_calcValue = calcYields[i];
      }
    };

    void set_calcEigenvalues(std::vector<double>& calcEigenvalues) {
      for (size_t i(0); i < calcEigenvalues.size(); ++i) {
        m_eigenvalues[i].m_calcValue = calcEigenvalues[i];
      }
    };

  };


  //****************************************************************************************************
  class ConditionsManager
  {
  public:

    ConditionsManager(System* pSys);

    // Reads the general bathgas from <me:conditions>
    bool ParseBathGas(PersistPtr ppConditions);
    // Reads the rest of  <me:conditions>
    bool ParseConditions();

    size_t      getNumPTPoints() const { return PandTs.size(); }
    size_t      getTotalNumPoints() const;
    double      PTPointTemp(int index) const { return PandTs[index].m_temperature; }
    double      PTPointRadT(int index) const { return PandTs[index].m_radiationTemp; }
    double      PTPointConc(int index) const { return PandTs[index].m_concentration; }
    const char* PTPointBathGas(int index) { return PandTs[index].m_pBathGasName.c_str(); }
    Precision   PTPointPrecision(int index) { return PandTs[index].m_precision; }
    map<Reaction*, double> PTPointExcessConcs(int index) { return PandTs[index].m_excessConcs; }

    // Calculate the effective reciprocal temperature,
    // taking into account a radiation field if there is one.
    double getBeta(int index) const {
      double collTemp = PandTs[index].m_temperature ;
      double radTemp  = PandTs[index].m_radiationTemp;
      double conc     = PandTs[index].m_concentration;
      double effTemp = collTemp;
      if (radTemp > 0.0) {
        effTemp = radTemp + (collTemp - radTemp) * conc / (1.e14 + conc);
      }
      return 1.0 / (boltzmann_RCpK * effTemp);
    }

    // An accessor method to get conditions and related properties for
    // use with plugins classes etc.
    bool getConditions(std::vector<double>& Temperature, std::vector<double>& Concentration);

    double getMaxTemperature();

    void get_experimentalRates(unsigned index, std::vector<conditionSet>& rates) const
    {
      rates = PandTs[index].m_rates;
    }

    void set_calculatedRates(unsigned index, std::vector<double>& calcRates)
    {
      PandTs[index].set_calcRates(calcRates);
    }

    void get_experimentalYields(unsigned index, std::vector<conditionSet>& yields) const
    {
      yields = PandTs[index].m_yields;
    }

    void set_calculatedYields(unsigned index, std::vector<double>& calcYields)
    {
      PandTs[index].set_calcYields(calcYields);
    }

    void get_experimentalEigenvalues(unsigned index, std::vector<conditionSet>& eigenvalues) const
    {
      eigenvalues = PandTs[index].m_eigenvalues;
    }

    void set_calculatedEigenvalues(unsigned index, std::vector<double>& calcEigenvalues)
    {
      PandTs[index].set_calcEigenvalues(calcEigenvalues);
    }

    std::vector<RawDataSet>& get_experimentalrawDataSets(unsigned index)
    {
      return PandTs[index].m_rawDataSets;
    }

    PersistPtr get_experimentalDataPtr(unsigned index, size_t i) const
    {
      return PandTs[index].m_expDataPtrs[i];
    }

    AnalysisData* get_analysisData(unsigned index)
    {
      return &(PandTs[index].m_analysisData);
    }

    GeneralAnalysisData* get_generalAnalysisData()
    {
      return &generalAnalysisData;
    }

    // Collect bath gas names from PandTs.
    void getAllBathGases(std::set<std::string>& bathGases);

    // Reconcile table across processs.
    void reconcileTable();

    // Calculate ChiSquared.
    void calculateChiSquared(double& chiSquared, vector<double>& residuals) const;

    // Write data table.
    void WriteDataTable() const;

    // Write calculated data to output.
    void AddCalcValToXml() const;

    // Write analysis data to <me:analysis> section of the XML output.
    void WriteAnalysisXML(PersistPtr m_ppIOPtr);

    // Write the general analysis data to <me:analysis> section of the XML output.
    void WriteXMLandClear();

    // Set the the index to the data set being analysed.
    void set_currentData(int i) { currentSet = i; };

    // Set the the index to the data set being analysed.
    int get_currentData() { return currentSet; };

  private:

    bool readPTs();

    bool ReadRange(const std::string& name,
      std::vector<double>& vals,
      PersistPtr            ppbase,
      bool                  MustBeThere = true);

    // Write calculated data to XML.
    void WriteDataToXml(PersistPtr pp, const vector<conditionSet>& data) const;

    // Normalize Experimental weights.
    void NormalizeExptWeights();

    System* m_pSys; //parent System

    PersistPtr m_ppConditions;

    PersistPtr m_ppAnalysis;

    // The excess reactant concentrations as specified in <Reaction>.
    // Provides unchanging base data when excess concs are individually
    // specified in <PTPair>
    std::map<Reaction*, double> baseExcessConcs; // Reaction, conc in ppcc

        // Paired concentration and pressure points.
    std::vector<CandTpair> PandTs;

    // Location of the parallel mananger.
    ParallelManager* m_pParallelManager;

    // General analysis data.
    GeneralAnalysisData generalAnalysisData;

    // Index of the current set being analysed.
    int currentSet;

  };

}//namespace
#endif
