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
  // default values
  static const double sigmaDefault            = 5.0;
  static const double epsilonDefault          = 50.0;

  //temporary flags for debug purposes
  static const bool collisionOperatorCheck    = false;
  static const bool debugFlag                 = false;

  //
  // Basic constants
  //
  static const double boltzmann_C             = 1.3806503e-23;  // Boltzmann constant (J*K-1).
  static const double SpeedOfLight_cm         = 2.99792458e+10 ;// speed of light in centimeter
  static const double PlancksConstant         = 6.6260689633e-34;   // Planck's constant in Joule times second
  static const double AvogadroC               = 6.0221367e+23;
  static const double CalorieInJoule          = 4.184;
  static const double AtmInMmHg               = 760.0;
  static const double pascalPerAtm            = 1.01325e+05;    //Pascal: N*m^-2

  //
  // Derived constants (do not alter contents below)
  //

  static const double boltzmann_RCpK          = boltzmann_C /(SpeedOfLight_cm * PlancksConstant); //0.695035612 ;
                                                                // Boltzmann constant (cm*K)-1. (Reciprocal Centimeter per Kelvin)
                                                                // boltzmann_C /(SpeedOfLight_cm * PlancksConstant);
  static const double kJPerMolInRC            = 1.0e3 / (AvogadroC * PlancksConstant * SpeedOfLight_cm);
                                                                // kilo Joule per mol to reciprocal centimeter
                                                                // 1.0e3 / (AvogadroC * PlancksConstant * SpeedOfLight_cm);
  //static const double KcalPerMolInRC          = kJPerMolInRC * CalorieInJoule;       // 349.757 ;
                                                                // kilo Calorie per mol in reciprocal centimeter

  static const double tp_C                    = pow((2.0 * M_PI * SpeedOfLight_cm) / (1.0e3 * PlancksConstant * AvogadroC * 1.0e4),1.5);    //3.24331e+20;
  static const double amu                     = 1.0 / (AvogadroC * 1.0e3);
  static const double idealGasC               = boltzmann_C * AvogadroC;

}
#endif // GUARD_Constants_h
