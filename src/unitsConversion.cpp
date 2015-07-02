//-------------------------------------------------------------------------------------------
//
// unitsConversion.cpp
//
//-------------------------------------------------------------------------------------------
#include "error.h"
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
      mass["H"]   =  1.007825032;
      mass["D"]   =  2.014101778;
      mass["2H"]  =  2.014101778;
      mass["3H"]  =  3.016049268;
      mass["C"]   = 12.000000000;
      mass["13C"] = 13.003354838;
      mass["14C"] = 14.003241991;
      mass["13N"] = 13.005738584;
      mass["N"]   = 14.003074007;
      mass["15O"] = 15.003065460;
      mass["O"]   = 15.994914620;
      mass["17O"] = 16.999131501;
      mass["F"]   = 18.998403220;
      mass["S"]   = 31.972071000;
      mass["Cl"]  = 34.968852680;
      mass["B"]   = 10.81;
      mass["Na"]  = 22.990;
      mass["Mg"]  = 24.305;
      mass["Al"]  = 26.982;
      mass["Si"]  = 28.085;
      mass["P"]   = 30.974;
      mass["K"]   = 39.098;
      mass["Ca"]  = 40.078;
      mass["Sc"]  = 44.956;
      mass["Ti"]  = 47.867;
      mass["V"]   = 50.942;
      mass["Cr"]  = 51.996;
      mass["Mn"]  = 54.938;
      mass["Fe"]  = 55.845;
      mass["Co"]  = 58.933;
      mass["Ni"]  = 58.693;
      mass["Cu"]  = 63.546;
      mass["Zn"]  = 65.38;
      mass["Ga"]  = 69.723;
      mass["Ge"]  = 72.63;
      mass["As"]  = 74.922;
      mass["Se"]  = 78.96;
      mass["Br"]  = 79.904;
      mass["Rb"]  = 85.468;
      mass["Sr"]  = 87.62;
      mass["Y"]   = 88.906;
      mass["Zr"]  = 91.224;
      mass["Nb"]  = 92.906;
      mass["Mo"]  = 95.96;
      mass["Tc"]  = 97.91;
      mass["Ru"]  = 101.07;
      mass["Rh"]  = 102.91;
      mass["Pd"]  = 106.43;
      mass["Ag"]  = 107.87;
      mass["Cd"]  = 112.41;
      mass["Ln"]  = 114.82;
      mass["Sn"]  = 118.71;
      mass["Sb"]  = 121.76;
      mass["Te"]  = 127.60;
      mass["I"]   = 126.9;

    }
    return mass[symb];
  }

  void initializeConversionMaps(){
    // define rule of conversion for concentration 
    concentrationMap["particles per cubic centimeter"] = 0;
    concentrationMap["PPCC"] = 0;
    concentrationMap["number density"] = 0;
    concentrationMap["Torr"] = 1;
    concentrationMap["torr"] = 1;
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

    timeMap["sec"]      = 1.0; //internal time units
    timeMap["msec"]     = 1.0e-3;
    timeMap["millisec"] = 1.0e-3;
    timeMap["usec"]     = 1.0e-6;
    timeMap["microsec"] = 1.0e-6;
    timeMap["nsec"]     = 1.0e-9;
    timeMap["nanosec"]  = 1.0e-9;
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
    if(unitInput.empty())
    {
      cinfo << "No energy units specified" << endl;
      return energyInput;
    }
    else if(energyMap.count(unitInput)==0) 
    {
      cerr << "Unrecognized energy unit: " + unitInput << endl;
      return 0.;
    }
    return energyInput * energyMap[unitInput];
  }

  double getConvertedTime(const string& unitInput, const double timeInput)
  {
    if (unitInput.empty())
    {
      cinfo << "No time units specified" << endl;
      return timeInput;
    }
    else if (timeMap.count(unitInput) == 0)
    {
      cerr << "Unrecognized time unit: " + unitInput << endl;
      return 0.;
    }
    return timeInput * timeMap[unitInput];
  }

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
  
  //Given energyInput in cm-1 returns value in specified units (no check on units validity)
  double ConvertFromWavenumbers(const string& unitInput, const double energyInput)
  {
    return energyInput / energyMap[unitInput];
  }

  //Converts the units of energyInput from unitInput to unitOutput
  double ConvertEnergy(const string& unitInput, const string& unitOutput, const double energyInput)
  {
    return energyInput * energyMap[unitInput]/ energyMap[unitOutput];
  }

  Precision txtToPrecision (const char *txt) {
	string strPrcsn(txt) ;
	Precision precision;
	if      (strPrcsn == "1d" || strPrcsn == "d"  || strPrcsn == "double")  
	  precision = DOUBLE ;
	else if (strPrcsn == "2d" || strPrcsn == "dd" || strPrcsn == "double-double") 
	  precision = DOUBLE_DOUBLE ;
	else if (strPrcsn == "4d" || strPrcsn == "qd" || strPrcsn == "quad-double")
	  precision = QUAD_DOUBLE ;
	else {
	  cerr << "Unknown precision. Calculation will be run using double precision" << endl ;
	  precision = DOUBLE ;
	}
	return precision ;
  }

}
