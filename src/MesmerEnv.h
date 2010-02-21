#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h

namespace mesmer
{
  struct MesmerEnv
  {
    MesmerEnv();

    // environmental values
    double beta;
    double conc;

    // granularity of the system
    int    GrainSize;                     // Grain size in cm-1
    int    MaxGrn;                        // The number of grains
    int    MaxCell;                       // The number of cells
    double MaximumTemperature;            // Maximum temperature for the purposes of setting the energy range
    double EMin, EMax;                    // The absolute lowest and highest energies in the system, cm-1
    double EAboveHill;                    // Max energy above the highest Hill [in kT]
    bool   useBasisSetMethod;             // Use the contracted basis set method.
    size_t nBasisSet;                     // Number of basis set functions to use.
  };

}//namespace


#endif // GUARD_MesmerEnv_h

