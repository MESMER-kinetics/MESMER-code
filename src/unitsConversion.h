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
  double atomMass(std::string symb);

  enum Precision {
	DOUBLE,
	DOUBLE_DOUBLE,
	QUAD_DOUBLE,
	UNDEFINED_PRECISION
  } ;

  Precision txtToPrecision (const char *txt) ;

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
