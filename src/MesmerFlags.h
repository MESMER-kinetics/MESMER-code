#ifndef GUARD_MesmerFlags_h
#define GUARD_MesmerFlags_h

namespace mesmer
{
  struct MesmerFlags
  {
    MesmerFlags();
    // whether do the fitting or grid search
    int    searchMethod; // 0 for single calculation, 1 for grid search, 2 for fitting, and more...?

    // decide what to report
    bool   testDOSEnabled;                // Whether to output test of DOS to mesmer.test
    bool   testRateConstantEnabled;       // Option to output canonical rate constant
    bool   microRateEnabled;              // Whether to output microcanonical rate coefficients
    bool   grainBoltzmannEnabled;         // Enabled printing grain boltzmann distribution
    bool   grainDOSEnabled;               // Enabled printing grain DOS
    bool   cyclePrintGrainDOS;            // Controls the print-out of grain DOS in each cycle (This is only for source term)
    bool   cellDOSEnabled;                // Enabled printing cell DOS
    bool   cyclePrintCellDOS;             // Controls the print-out of cell DOS in each cycle (This is only for source term)
    bool   reactionOCSEnabled;           // Enabled printing collision operator column Sums
    bool   kfEGrainsEnabled;              // Enabled printing k_f(E) grains
    bool   kbEGrainsEnabled;              // Enabled printing k_b(E) grains
    bool   TunnellingCoeffEnabled;        // Enabled printing Tunneling coefficients
    bool   cellFluxEnabled;             // Enabled printing transition state flux
    bool   grainFluxEnabled;            // Enabled printing transition state flux
    bool   rateCoefficientsOnly;          // Calculate rate coefficients only without doing collision operators
    bool   useTheSameCellNumber;          // Option to use the same cell number or not in various conditions
    bool   grainedProfileEnabled;         // Option to print out grained species profile (before summation to individual species)
    bool   speciesProfileEnabled;         // Option to print species profile
    bool   viewEvents;                    // Print events timestamps
    double maxEvolutionTime;              // Maximum time of evolution for the species profile
    int printEigenValuesNum;              // Number of eigen values to be printed: -1 for all of them, otherwise specified.
    int printReactionOperatorNum;         // Size of printed reaction operator before and after diagonalization: -1
                                          // for all of them, -2 for 1/2 of them, -3 for 1/3 of them, otherwise specified
                                          // by positive integers.
    bool allowSmallerDEDown;              // decide whether allows <delta E>d to be smaller than grain size.
    bool useFFTIntegration;               // integrates numbers by FFT.
    bool print_TabbedMatrices;            // print tabbed instead of fixed-widthed matrices.
  };
}//namespace


#endif // GUARD_MesmerFlags_h

