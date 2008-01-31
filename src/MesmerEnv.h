#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h
#include "oberror.h"

namespace mesmer
{
  struct MesmerEnv
  {
    MesmerEnv();

    // environmental values
    double beta;
    double conc;

    // granularity of the system
    double GrainSize;  //Grain size in cm-1
    int    MaxGrn;     //The number of grains
    int    MaxCell;    //The number of cells
    double MaxT;       //Maximum temperature for the purposes of setting the energy range
    double EMin, EMax; // The absolute lowest and highest energies in the system, cm-1

    // decide what to report
    bool   testDOSEnabled;    // Whether to output test of DOS to mesmer.test
    bool   microRateEnabled;  // Whether to output microcanonical rate coefficients
    bool   grainBoltzmannEnabled;   // Enabled printing grain boltzmann distribution
    bool   grainDOSEnabled;   // Enabled printing grain DOS
    bool   cellDOSEnabled;    // Enabled printing cell DOS
    bool   collisionOCSEnabled; // Enabled printing collision operator column Sums
    bool   kfECellsEnabled;     // Enabled printing k_f(E) cells
    bool   kfEGrainsEnabled;    // Enabled printing k_f(E) grains
    bool   kbECellsEnabled;     // Enabled printing k_b(E) cells
    bool   kbEGrainsEnabled;    // Enabled printing k_b(E) grains
    double EAboveWell; //Max energy above the highest well [was 10kT]
    int printEigenValuesNum; // number of eigen values to be printed: -1 for all of them, otherwise specified.
  };
}//namespace


#endif // GUARD_MesmerEnv_h

