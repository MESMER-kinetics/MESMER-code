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

namespace Constants {

    //
    // Constants in terms of wavenumbers.
    //

    static const double boltzmann    = 0.695029 ;

    static const double KCMLTOPCM    = 349.757 ;

    static const double plancksConst = 1.0/2.998e+10 ; // In wavenumber units this comes out as
	                                                   // the reciprocal of the speed of light.
}

#endif // GUARD_Constants_h
