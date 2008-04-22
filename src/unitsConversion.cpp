//-------------------------------------------------------------------------------------------
//
// unitsConversion.cpp
//
//-------------------------------------------------------------------------------------------
#include "unitsConversion.h"
#include "Constants.h"
using namespace std;
using namespace Constants;

namespace mesmer
{ 
  void initializeConversionMaps(){
    // give a name to return a number for the switch functions
    // define rule of conversion for concentration 
    concentrationMap["particles per cubic centimeter"] = 0;
    concentrationMap["PPCC"] = 0;
    concentrationMap["number density"] = 0;
    concentrationMap["Torr"] = 1;
    concentrationMap["mmHg"] = 1;
    concentrationMap["mbar"] = 2;
    concentrationMap["atm"] = 3;
    concentrationMap["psi"] = 4;

    // define rule of conversion for energy
    energyMap["wavenumber"] = 0;
    energyMap["cm-1"] = 0;
    energyMap["kJ/mol"] = 1;
    energyMap["kJ per mol"] = 1;
    energyMap["kcal/mol"] = 2;
    energyMap["kcal per mol"] = 2;
    energyMap["Hartree"] = 3;
    energyMap["au"] = 3;
  }

  // Returns particles per cubic centimeter no matter what unit the user has provided.
  double getConvertedCP(string unitInput, double concentrationInput, double temperatureInput)
  {
    
    // switch
    switch (concentrationMap[unitInput])
    {
      case 0: // number density
        return concentrationInput;
      case 1: // Torr
        return ((concentrationInput / Atm_in_Torr) * atm_in_pascal * AvogadroC / (idealGasC * temperatureInput * 1.0e6));
      case 2: // mbar
        return (concentrationInput * mbar_in_pascal * AvogadroC / (idealGasC * temperatureInput * 1.0e6));
      case 3: // atm
        return (concentrationInput * atm_in_pascal * AvogadroC / (idealGasC * temperatureInput * 1.0e6));
      case 4: // psi
        return (concentrationInput * psi_in_pascal * AvogadroC / (idealGasC * temperatureInput * 1.0e6));
    }
    return 0.;
  }

  // Returns particles per cubic centimeter no matter what unit the user has provided.
  double getConvertedEnergy(string unitInput, double energyInput)
  {
    
    // switch
    switch (energyMap[unitInput])
    {
      case 0: // wavenumber
        return energyInput;
      case 1: // kJ/mol
        return energyInput * kJPerMol_in_RC;
      case 2: // kcal/mol
        return energyInput * kCalPerMol_in_RC;
      case 3: // Hartree
        return energyInput * Hartree_In_kJperMol * kJPerMol_in_RC;
    }
    return 0.;
  }
}
