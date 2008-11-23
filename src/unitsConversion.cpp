//-------------------------------------------------------------------------------------------
//
// unitsConversion.cpp
//
//-------------------------------------------------------------------------------------------
#include "oberror.h"
#include "unitsConversion.h"
using namespace std;
using namespace Constants;

namespace mesmer
{ 
  void initializeConversionMaps(){
    // define rule of conversion for concentration 
    concentrationMap["particles per cubic centimeter"] = 1.0; //internal concentration units
    concentrationMap["PPCC"] = 1.0;
    concentrationMap["number density"] = 1.0;
    concentrationMap["Torr"] = atm_in_pascal * AvogadroC / (idealGasC * Atm_in_Torr * 1.0e6);
    concentrationMap["mmHg"] = atm_in_pascal * AvogadroC / (idealGasC * Atm_in_Torr * 1.0e6);
    concentrationMap["mbar"] = mbar_in_pascal * AvogadroC / (idealGasC * 1.0e6);
    concentrationMap["atm"]  = atm_in_pascal * AvogadroC / (idealGasC * 1.0e6);
    concentrationMap["psi"]  = psi_in_pascal * AvogadroC / (idealGasC * 1.0e6);

    // define rule of conversion for energy
    energyMap["wavenumber"]   = 1.0; //internal energy units
    energyMap["cm-1"]         = 1.0;
    energyMap["kJ/mol"]       = kJPerMol_in_RC;
    energyMap["kJ per mol"]   = kJPerMol_in_RC;
    energyMap["kcal/mol"]     = kCalPerMol_in_RC;
    energyMap["kcal per mol"] = kCalPerMol_in_RC;
    energyMap["Hartree"]      = Hartree_In_kJperMol * kJPerMol_in_RC;
    energyMap["au"]           = Hartree_In_kJperMol * kJPerMol_in_RC;
  }

  // Returns particles per cubic centimeter no matter what unit the user has provided.
  double getConvertedP(const string& unitInput, const double concentrationInput, const double temperatureInput)
  {
    if(concentrationMap.count(unitInput)==0) {
      cerr << "Unrecognized concentration unit: " + unitInput << endl;
      return 0.;
    }
    return concentrationInput * concentrationMap[unitInput] /temperatureInput;
    
/*     
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
*/
  }

  // Returns particles cm-1 no matter what unit the user has provided.
  double getConvertedEnergy(const string& unitInput, const double energyInput)
  {
    if(energyMap.count(unitInput)==0) {
      cerr << "Unrecognized energy unit: " + unitInput << endl;
      return 0.;
    }
    return energyInput * energyMap[unitInput];
    
/*    // switch
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
*/
  }
  //Given energyInput in cm-1 returns value in specified units (no check on units validity)
  double ConvertFromWavenumbers(const string& unitInput, const double energyInput)
  {
    return energyInput / energyMap[unitInput];
  }    

}
