#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h

namespace mesmer
{
  struct MesmerEnv
  {
    MesmerEnv();
    // whether do the fitting or grid search
    int    searchMethod; // 0 for single calculation, 1 for grid search, 2 for fitting, and more...?

    // environmental values
    double beta;
    double conc;

    // granularity of the system
    int    GrainSize;                     //Grain size in cm-1
    int    MaxGrn;                        //The number of grains
    int    MaxCell;                       //The number of cells
    double MaximumTemperature;            //Maximum temperature for the purposes of setting the energy range
    double EMin, EMax;                    // The absolute lowest and highest energies in the system, cm-1

    // decide what to report
    bool   testDOSEnabled;                // Whether to output test of DOS to mesmer.test
    bool   testRateConstantEnabled;       // Option to output canonical rate constant
    bool   microRateEnabled;              // Whether to output microcanonical rate coefficients
    bool   grainBoltzmannEnabled;         // Enabled printing grain boltzmann distribution
    bool   grainDOSEnabled;               // Enabled printing grain DOS
    bool   cellDOSEnabled;                // Enabled printing cell DOS
    bool   collisionOCSEnabled;           // Enabled printing collision operator column Sums
    bool   kfEGrainsEnabled;              // Enabled printing k_f(E) grains
    bool   kbEGrainsEnabled;              // Enabled printing k_b(E) grains
    bool   TunnellingCoeffEnabled;         // Enabled printing Tunneling coefficients
    bool   cellTSFluxEnabled;             // Enabled printing transition state flux
    bool   grainTSFluxEnabled;            // Enabled printing transition state flux
    bool   rateCoefficientsOnly;          // Calculate rate coefficients only without doing collision operators
    bool   useTheSameCellNumber;          // Option to use the same cell number or not in various conditions
    bool   grainedProfileEnabled;         // Option to print out grained species profile (before summation to individual species)
    bool   speciesProfileEnabled;         // Option to print species profile
    double EAboveHill;                    // Max energy above the highest Hill [in kT]
    double maxEvolutionTime;              // Maximum time of evolution for the species profile
    int printEigenValuesNum;              // Number of eigen values to be printed: -1 for all of them, otherwise specified.
  };
}//namespace


#endif // GUARD_MesmerEnv_h

