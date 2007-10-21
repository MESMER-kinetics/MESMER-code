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

#ifndef M_PI
#define M_PI acos(-1.)
#endif //M_PI

namespace Constants {
  //
  // Constants in terms of wavenumbers.
  //
  static const double boltzmann_RCpK  = 0.695035612 ;   // Boltzmann constant (cm*K)-1. (Reciprocal Centimeter per Kelvin)
  static const double KcalPerMolToRC  = 349.757 ;       // kilo Calorie per mol to reciprocal centimeter
  
  //static const double OneOverSpeedOfLightInCentimeter = 1.0/2.998e+10 ; // In wavenumber units this comes out as
  //                                                                      // the reciprocal of the speed of light.
  static const double SpeedOfLight_cm = 2.99792458e+10 ; // speed of light in centimeter

  static const double kJPerMolToRC    = 83.593461;      // kilo Joule per mol to reciprocal centimeter
  static const double boltzmannC      = 1.3806503e-23;  // m2 kg s-2 K-1
  static const double AvogadroC       = 6.0221367e+23;
  static const double pascalPerAtm    = 1.01325e+05;
  static const double idealGasC       = 8.31451l;
}

#endif // GUARD_Constants_h
