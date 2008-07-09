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
#include "marray.h"

using namespace std;

namespace mesmer
{
  struct conditionSet{
  public:

    conditionSet(string ref1_, string ref2_, double value_, double error_):
      ref1(ref1_), ref2(ref2_), value(value_), error(error_)
    {}

    void get_conditionSet(string& ref1_, string& ref2_, double& value_, double& error_){
      ref1_ = ref1;
      ref2_ = ref2;
      value_ = value;
      error_ = error;
    }

  private:

    string ref1;
    string ref2;
    double value;
    double error;

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

    void set_condition(string ref1_, string ref2_, double value_, double error_){
      conditionSet tCS(ref1_, ref2_, value_, error_);
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

  };

  // defining the unit concentration units conversion rule
  typedef std::map<std::string, int> stringMapType;

  // mapping the conversion of concentration, pressure
  static stringMapType concentrationMap;

  // mapping the conversion of energy
  static stringMapType energyMap;

  void initializeConversionMaps();
  double getConvertedP(string unitInput, double concentrationInput, double temperatureInp);
  double getConvertedEnergy(string unitInput, double energyInput);

}//namespace

#endif // GUARD_unitsConversion_h
