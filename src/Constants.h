#ifndef GUARD_Constants_h
#define GUARD_Constants_h

//-------------------------------------------------------------------------------------------
//
// Constants.h
//
// Author: Struan Robertson
// Date:   10/Mar/2003
//
// This header file contains the definitions of common physical constants.
//
//-------------------------------------------------------------------------------------------
#include <cmath>
#include<limits>

#ifndef M_PI
#define M_PI acos(-1.)
#endif //M_PI

namespace Constants {
  //More meaningful synonym used in parsing
  static const bool   optional = false;
  // default values
  static const double sigmaDefault = 5.0; // Angstroms.
  static const double epsilonDefault = 50.0; // K.

  //
  // Basic constants
  //
  static const double boltzmann_C = 1.3806503e-23;    // Boltzmann constant (J*K-1).
  static const double SpeedOfLight_in_cm = 2.99792458e+10;  // speed of light in centimeter/s
  static const double PlancksConstant_in_JouleSecond = 6.6260689633e-34; // Planck's constant in Joule times second
  static const double AvogadroC = 6.0221367e+23;
  static const double Calorie_in_Joule = 4.184;
  static const double Atm_in_Torr = 760.0 / 1.01325;   // The divison by 1.01325 is to account for the change in definition of 1 atom to be 1 bar.
  static const double mbar_in_pascal = 100.;
  static const double bar_in_pascal = 1.e+05;
  static const double atm_in_pascal = 1.e+05;           //Pascal: N*m^-2
  static const double Hartree_In_kJperMol = 2625.478;
  static const double psi_in_pascal = 6894.76;
  static const double bohr_in_angstrom = 0.529177;

  //
  // Derived constants (do not alter contents below)
  //

  static const double boltzmann_RCpK = boltzmann_C / (SpeedOfLight_in_cm * PlancksConstant_in_JouleSecond); //0.695035612 ;
                                                                // Boltzmann constant (cm*K)-1. (Reciprocal Centimeter per Kelvin)
                                                                // boltzmann_C /(SpeedOfLight_in_cm * PlancksConstant_in_JouleSecond);
  static const double kJPerMol_in_RC = 1.0e3 / (AvogadroC * PlancksConstant_in_JouleSecond * SpeedOfLight_in_cm);
  // kilo Joule per mol to reciprocal centimeter
  // 1.0e3 / (AvogadroC * PlancksConstant_in_JouleSecond * SpeedOfLight_in_cm);
  static const double kCalPerMol_in_RC = kJPerMol_in_RC * Calorie_in_Joule;       // 349.757 ;
                                                                // kilo Calorie per mol in reciprocal centimeter

  static const double Hartree_in_RC = kJPerMol_in_RC * Hartree_In_kJperMol;    // 2.1947e5 ;
                                                                // kilo Calorie per mol in reciprocal centimeter


  static const double tp_C = pow((2.0 * M_PI * SpeedOfLight_in_cm) / (1.0e3 * PlancksConstant_in_JouleSecond * AvogadroC * 1.0e4), 1.5);    //3.24331e+20 cm-3;
  static const double amu = 1.0 / (AvogadroC * 1.0e3); // Atomic mass unit in kg
  static const double idealGasC = boltzmann_C * AvogadroC;

  // Convertion factor needed to obtain rotational constant (cm-1) from moment of Inertia (amu Ang^2).
  // The definitiion used here is h/(8 pi^2c), but note, authors differ on this definition.

  // static const double conMntInt2RotCnt = 16.85917 ;
  static const double conMntInt2RotCnt = PlancksConstant_in_JouleSecond * 1.0e+20 / (8.0 * M_PI * M_PI * SpeedOfLight_in_cm * amu);

  // Convertion factor needed to obtain frequencies (cm-1) from mass weighted Hessain (kJ/mol/amu/Ang^2).

  static const double conHess2Freq = 1.e13 / SpeedOfLight_in_cm;
}
#endif // GUARD_Constants_h
