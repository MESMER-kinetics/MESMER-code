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
#include "Constants.h"

using namespace std;

namespace mesmer
{
  // to make sure if there is a concentration or pressure definition, there is a temperature definition.
  struct CandTpair{
    CandTpair(double _cp, double _t): concentration(_cp), temperature(_t), precision(0){}
    CandTpair(double _cp, double _t, int _pre): concentration(_cp), temperature(_t), precision(_pre){}
    double concentration; // particles per cubic centimeter
    double temperature; // Kelvin
    int precision;
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
