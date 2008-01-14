#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h
#include "oberror.h"

namespace mesmer
{
  struct MesmerEnv
  {
    MesmerEnv();

    double beta;
    double conc;

    double GrainSize;  //Grain size in cm-1
    int    MaxGrn;     //The number of grains
    int    MaxCell;    //The number of cells
    double MaxT;       //Maximum temperature for the purposes of setting the energy range
    
    double EMin, EMax; // The absolute lowest and highest energies in the system, cm-1

    bool   microRateEnabled;  // Whether to output microcanonical rate coefficients
    bool   grainDOSEnabled;   // Enabled printing grain DOS
    bool   cellDOSEnabled;   // Enabled printing cell DOS
    bool   collisionOCSEnabled; // Enabled printing collision operator column Sums
    bool   kECellsEnabled;    // Enabled printing k(E) cells
    double EAboveWell; //Max energy above the highest well [was 10kT]
    int printEigenValuesNum; // number of eigen values to be printed: -1 for all of them, otherwise specified.
  };
}//namespace


#endif // GUARD_MesmerEnv_h

