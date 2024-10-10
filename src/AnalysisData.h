#ifndef GUARD_AnalysisData_h
#define GUARD_AnalysisData_h

#include <string>
#include <sstream>
#include <vector>
#include "dMatrix.h"

namespace mesmer
{
  // This structure holds condition specific data.

  struct AnalysisData
  {
  public:
    AnalysisData() : m_number(0),
      m_selection(""),
      m_eigenvalues(),
      m_lossRef(),
      m_lossRateCoeff(),
      m_firstOrderReactionType(),
      m_firstOrderFromRef(),
      m_firstOrderToRef(),
      m_firstOrderRateCoeff(),
      m_secondOrderReactionType(),
      m_secondOrderFromRef(),
      m_secondOrderToRef(),
      m_secondOrderRateCoeff(),
      m_timePoints(),
      m_aveEnergyRef(),
      m_aveEnergy(),
      m_PopRef(),
      m_Pop(),
      m_grnTimes(),
      m_grnDists(),
      m_warning() {}

    ~AnalysisData() {};

    void clear() {
      m_eigenvalues.clear();
      m_lossRef.clear();
      m_lossRateCoeff.clear();
      m_firstOrderReactionType.clear();
      m_firstOrderFromRef.clear();
      m_firstOrderToRef.clear();
      m_firstOrderRateCoeff.clear();
      m_secondOrderReactionType.clear();
      m_secondOrderFromRef.clear();
      m_secondOrderToRef.clear();
      m_secondOrderRateCoeff.clear();
      m_timePoints.clear();
      m_aveEnergyRef.clear();
      m_aveEnergy.clear();
      m_PopRef.clear();
      m_Pop.clear();
      m_grnTimes.clear();
      m_grnDists.clear();
      m_warning.clear();
    };

    // Parameters
    std::string m_parameters;

    // Eigenvalues.
    int m_number;
    std::string m_selection;
    std::vector<double> m_eigenvalues;

    // Loss rate coefficients.
    std::vector<std::string> m_lossRef;
    std::vector<double> m_lossRateCoeff;

    // First order rate coefficients.
    std::vector<std::string> m_firstOrderReactionType;
    std::vector<std::string> m_firstOrderFromRef;
    std::vector<std::string> m_firstOrderToRef;
    std::vector<double> m_firstOrderRateCoeff;

    // Second order rate coefficients.
    std::vector<std::string> m_secondOrderReactionType;
    std::vector<std::string> m_secondOrderFromRef;
    std::vector<std::string> m_secondOrderToRef;
    std::vector<double> m_secondOrderRateCoeff;

    // Average Energy and population.
    std::vector<double> m_timePoints;
    std::vector<std::string> m_aveEnergyRef;
    std::vector<double> m_aveEnergy;
    std::vector<std::string> m_PopRef;
    std::vector<double> m_Pop;

    // Grain populations.
    std::vector<double> m_grnTimes;
    std::map<string, std::vector<std::vector<double> > >  m_grnDists;

    // Warning.
    std::string m_warning;
  };

  // This structure holds data that applied to all conditions.

  class GeneralAnalysisData
  {
  public:
    GeneralAnalysisData() : m_covariance(1), m_writeCovar(false) {}

    ~GeneralAnalysisData() { clear(); };

    void clear() {
      m_covariance.resize(1);
    };

    void setCovariance(qdMatrix &covariance) {
      m_covariance = covariance;
      m_writeCovar = true;
    };

    void writeCovariance(PersistPtr ppAnalysis) const {
      if (m_writeCovar) {
        PersistPtr ppCovariance = ppAnalysis->XmlWriteElement("me:covariance", "");
        m_covariance.WriteToXML(ppCovariance);
      }
    };

  private:

    // Covariance.
    qdMatrix m_covariance;
    bool m_writeCovar;

  };

}//namespace

#endif // GUARD_AnalysisData_h

