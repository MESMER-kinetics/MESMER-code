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
  //Lookup table for atomic masses
  double atomMass(std::string symb)
  {
    static std::map<std::string, double> mass;
    if(mass.empty())
    {
      //Masses of most common isotope
      // reference page http://pntpm.ulb.ac.be/private/divers.htm
      mass["H"]  =  1.007825032;
      mass["D"]  =  2.014101778;
      mass["2H"] =  2.014101778;
      mass["3H"] =  3.016049268;
      mass["C"]  = 12.000000000;
      mass["13C"]= 13.003354838;
      mass["14C"]= 14.003241991;
      mass["13N"]= 13.005738584;
      mass["N"]  = 14.003074007;
      mass["15O"]= 15.003065460;
      mass["O"]  = 15.994914620;
      mass["17O"]= 16.999131501;
      mass["F"]  = 18.998403220;
      mass["S"]  = 31.972071000;
      mass["Cl"] = 34.968852680;
    }
    return mass[symb];
  }

  void initializeConversionMaps(){
    // define rule of conversion for concentration 
    concentrationMap["particles per cubic centimeter"] = 0;
    concentrationMap["PPCC"] = 0;
    concentrationMap["number density"] = 0;
    concentrationMap["Torr"] = 1;
    concentrationMap["mmHg"] = 1;
    concentrationMap["mbar"] = 2;
    concentrationMap["atm"] = 3;
    concentrationMap["psi"] = 4;
    concentrationMap["moles/cc"] =5;

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
      case 5: //moles/cc
        return concentrationInput * AvogadroC;
    }
    return 0.;

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
