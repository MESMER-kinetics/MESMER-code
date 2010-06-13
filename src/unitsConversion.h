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

    CandTpair(double cp_, double t_): concentration(cp_), temperature(t_), precision(0){}
    CandTpair(double cp_, double t_, int _pre): concentration(cp_), temperature(t_), precision(_pre){}

    // Accessors
    const double get_concentration(){ return concentration; }
    const double get_temperature()  { return temperature;   }
    const int    get_precision()    { return precision;     }

    void set_experimentalRate(string ref1, string ref2, string refReaction, double value, double error){
      conditionSet tCS(ref1, ref2, refReaction, value, error);
      rates.push_back(tCS);
    }

    void get_experimentalRates(std::vector<conditionSet>& rates_){
      rates_ = rates;
    }

  private:

    double  concentration; // particles per cubic centimeter
    double  temperature; // Kelvin
    int     precision;

    vector<conditionSet> rates;
    vector<conditionSet> yields;
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
