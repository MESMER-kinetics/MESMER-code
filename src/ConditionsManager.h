#ifndef GUARD_ConditionsManager_h
#define GUARD_ConditionsManager_h

#include <vector>
#include "System.h"

namespace mesmer
{

struct conditionSet
{
  public:
  conditionSet(std::string ref1, std::string ref2, std::string refReaction, double value, double error):
    m_ref1(ref1), m_ref2(ref2), m_refReaction(refReaction), m_value(value), m_error(error)
    {}

    void get_conditionSet(std::string& ref1, std::string& ref2, std::string& refReaction, double& value, double& error) {
      ref1        = m_ref1 ;
      ref2        = m_ref2 ;
      refReaction = m_refReaction ;
      value       = m_value ;
      error       = m_error ;
    }

  private:
  std::string m_ref1;
  std::string m_ref2;
  std::string m_refReaction;
  double m_value;
  double m_error;
};

struct RawDataSet
{
  std::vector<std::pair<double, double>> data;
  const char* m_Name;
  double m_StartTime;           //sec
};

// To make sure if there is a concentration or pressure definition, there is a temperature definition.
struct CandTpair{

  double      m_concentration; // particles per cubic centimeter
  double      m_temperature; // Kelvin
  Precision   m_precision;
  std::string m_pBathGasName;
  std::map<Reaction*, double> m_excessConcs; // Reaction, conc in ppcc

  std::vector<conditionSet> m_rates;
  std::vector<conditionSet> m_yields;
  std::vector<conditionSet> m_eigenvalues;
  std::vector<PersistPtr> m_expDataPtrs;
  std::vector<RawDataSet> m_rawDataSets;

  CandTpair(double cp_, double t_, Precision _pre, const char* _bathGas,
            const map<Reaction*,double>& _excessConcs)
            : m_concentration(cp_), m_temperature(t_), m_precision(_pre),
              m_pBathGasName(_bathGas), m_excessConcs(_excessConcs) {}

  void set_experimentalRates(PersistPtr ppData, std::string ref1, std::string ref2,
                             std::string refReaction, double value, double error) {
    m_rates.push_back(conditionSet(ref1, ref2, refReaction, value, error));
    m_expDataPtrs.push_back(ppData);
  }

  void set_experimentalYields(PersistPtr ppData, std::string ref,  std::string yieldTime, double value, double error) {
    m_yields.push_back(conditionSet(ref, yieldTime, std::string(""), value, error)) ;
    m_expDataPtrs.push_back(ppData);
  }

  void set_experimentalEigenvalues(PersistPtr ppData, std::string eigenvalueID, double value, double error) {
    m_eigenvalues.push_back(conditionSet(std::string(""), std::string(""), eigenvalueID, value, error)) ;
    m_expDataPtrs.push_back(ppData);
  }
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

  size_t      getNumPTPoints() const            { return PandTs.size(); }
  double      PTPointTemp(int index)  { return PandTs[index].m_temperature; }
  double      PTPointConc(int index){ return PandTs[index].m_concentration;}
  const char* PTPointBathGas(int index)      {return PandTs[index].m_pBathGasName.c_str();}
  Precision   PTPointPrecision(int index)    {return PandTs[index].m_precision;}
  map<Reaction*,double> PTPointExcessConcs(int index)  {return PandTs[index].m_excessConcs; }

  // An accessor method to get conditions and related properties for
  // use with plugins classes etc.
  bool getConditions (std::vector<double> &Temperature, std::vector<double> &Concentration) ;

  double getMaxTemperature();

  void get_experimentalRates(unsigned index, std::vector<conditionSet>& rates) const
  { rates = PandTs[index].m_rates; }

  void get_experimentalYields(unsigned index, std::vector<conditionSet>& yields) const
  { yields =  PandTs[index].m_yields ; }

  void get_experimentalEigenvalues(unsigned index, std::vector<conditionSet>& eigenvalues) const
  { eigenvalues =  PandTs[index].m_eigenvalues ; }

  PersistPtr get_experimentalDataPtr(unsigned index, size_t i) const
  { return PandTs[index].m_expDataPtrs[i]; }

  //Collect bath gas names from PandTs
  void getAllBathGases(std::set<std::string>& bathGases);

  void AddRawDataSet(RawDataSet ds){ RawDataSets.push_back(ds); }
  RawDataSet GetRawDataSet(const unsigned index) const { return RawDataSets[index]; }

private:
  bool readPTs();

  bool ReadRange(const std::string&    name,
    std::vector<double>&  vals,
    PersistPtr            ppbase,
    bool                  MustBeThere=true);

  System* m_pSys; //parent System
  PersistPtr m_ppConditions;

  // The excess reactant concentrations as specified in <Reaction>.
  // Provides unchanging base data when excess concs are individually
  // specified in <PTPair>
  std::map<Reaction*, double> baseExcessConcs; // Reaction, conc in ppcc

  // Paired concentration and pressure points.
  std::vector<CandTpair> PandTs;

  // All the raw data sets in this conditions block
  std::vector < RawDataSet > RawDataSets;

};


}//namespace
#endif
