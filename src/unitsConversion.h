#ifndef GUARD_unitsConversion_h
#define GUARD_unitsConversion_h

//-------------------------------------------------------------------------------------------
//
// unitsConversion.h
//
// Defines units conversion of concentration and energy.
//-------------------------------------------------------------------------------------------
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include "Constants.h"
#include "Persistence.h"

using namespace std;

namespace mesmer
{
  //Global variable for the XML addresses of range variables. Definition in unitsConversion.cpp.
  extern std::vector<PersistPtr> RangeXmlPtrs;

  double atomMass(std::string symb);

  enum Precision {
    DOUBLE,
    DOUBLE_DOUBLE,
    QUAD_DOUBLE,
    UNDEFINED_PRECISION
  } ;

  struct conditionSet{
  public:

    conditionSet(string ref1, string ref2, string refReaction, double value, double error):
        m_ref1(ref1), m_ref2(ref2), m_refReaction(refReaction), m_value(value), m_error(error)
        {}

        void get_conditionSet(string& ref1, string& ref2, string& refReaction, double& value, double& error) {
          ref1        = m_ref1 ;
          ref2        = m_ref2 ;
          refReaction = m_refReaction ;
          value       = m_value ;
          error       = m_error ;
        }

  private:

    string m_ref1;
    string m_ref2;
    string m_refReaction;
    double m_value;
    double m_error;

  };

  // to make sure if there is a concentration or pressure definition, there is a temperature definition.
  class CandTpair{

  public:

    CandTpair(double cp_, double t_): m_concentration(cp_), m_temperature(t_), m_precision(DOUBLE){}
    CandTpair(double cp_, double t_, Precision _pre): m_concentration(cp_), m_temperature(t_), m_precision(_pre){}

    // Accessors
    const double    get_concentration(){ return m_concentration; }
    const double    get_temperature()  { return m_temperature;   }
    const Precision get_precision()    { return m_precision;     }

    void set_experimentalRates(string ref1, string ref2, string refReaction, double value, double error){
      m_rates.push_back(conditionSet(ref1, ref2, refReaction, value, error));
    }

    void get_experimentalRates(std::vector<conditionSet>& rates){ rates = m_rates; }

    void set_experimentalYields(string ref1, string ref2, double value, double error){
      m_yields.push_back(conditionSet(ref1, ref2, string(""), value, error)) ;
    }

    void get_experimentalYields(std::vector<conditionSet>& yields){ yields = m_yields ; }

    void set_experimentalEigenvalues(string ref1, string ref2, double value, double error){
      m_eigenvalues.push_back(conditionSet(ref1, ref2, string(""), value, error)) ;
    }

    void get_experimentalEigenvalues(std::vector<conditionSet>& eigenvalues){ eigenvalues = m_eigenvalues ; }

  private:

    double    m_concentration; // particles per cubic centimeter
    double    m_temperature; // Kelvin
    Precision m_precision;

    vector<conditionSet> m_rates;
    vector<conditionSet> m_yields;
    vector<conditionSet> m_eigenvalues;
  };

  // mapping the conversion of concentration, pressure
  static std::map<std::string, int> concentrationMap;

  // mapping the conversion of energy
  static std::map<std::string, double> energyMap;

  void initializeConversionMaps();
  double getConvertedP(const string& unitInput, const double concentrationInput, const double temperatureInp);
  double getConvertedEnergy(const string& unitInput, const double energyInput);
  double ConvertFromWavenumbers(const string& unitInput, const double energyInput);

}//namespace

#endif // GUARD_unitsConversion_h
