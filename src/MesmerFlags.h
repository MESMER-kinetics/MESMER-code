#ifndef GUARD_MesmerFlags_h
#define GUARD_MesmerFlags_h

#include <string>
#include <sstream>

namespace mesmer
{

  struct MesmerFlags
  {
    MesmerFlags();

    // decide what to report
    bool   testDOSEnabled;                 // Whether to output test of DOS to mesmer.test
    bool   testRateConstantEnabled;        // Option to output canonical rate constant
    bool   microRateEnabled;               // Whether to output microcanonical rate coefficients
    bool   grainBoltzmannEnabled;          // Enabled printing grain boltzmann distribution
    bool   grainDOSEnabled;                // Enabled printing grain DOS
    bool   grainTSsosEnabled;			   // enables printing of TS sum of states
    bool   cyclePrintGrainDOS;             // Controls the print-out of grain DOS in each cycle (This is only for source term)
    bool   cellDOSEnabled;                 // Enabled printing cell DOS
    bool   cyclePrintCellDOS;              // Controls the print-out of cell DOS in each cycle (This is only for source term)
    bool   reactionOCSEnabled;             // Enabled printing reaction operator column Sums
    bool   kfEGrainsEnabled;               // Enabled printing k_f(E) grains
    bool   kbEGrainsEnabled;               // Enabled printing k_b(E) grains
    bool   TunnellingCoeffEnabled;         // Enabled printing Tunneling coefficients
    bool   CrossingCoeffEnabled;           // Enabled printing Crossing coefficients
    bool   cellFluxEnabled;                // Enabled printing transition state flux
    bool   grainFluxEnabled;               // Enabled printing transition state flux
    bool   rateCoefficientsOnly;           // Calculate rate coefficients only without doing collision operators
    bool   useTheSameCellNumber;           // Option to use the same cell number or not in various conditions
    bool   grainedProfileEnabled;          // Option to print out grained species profile (before summation to individual species)
    bool   speciesProfileEnabled;          // Option to print species profile.
	bool   printPhenomenologicalEvolution; // Option to print species profile calculated from phenomenological equations.
    bool   InitialDistEnabled;             // Option to print initial distribution.
    bool   viewEvents;                     // Print events timestamps
    double shortestTimeOfInterest;         // Shortest time of interest
    double maxEvolutionTime;               // Maximum time of evolution for the species profile
    int    printEigenValuesNum;            // Number of eigen values to be printed: -1 for all of them, otherwise specified.
	bool   printEigenVectors;              // If true and eigenvalues are printed, then print the associated eigenvectors.
    int    printReactionOperatorNum;       // Size of printed reaction operator before and after diagonalization: -1
                                           // for all of them, -2 for 1/2 of them, -3 for 1/3 of them, otherwise specified
                                           // by positive integers.
	bool   printSinkFluxes;                // Prints the fluxes into each sink from each isomer or pseudoisomer.
    bool   allowSmallerDEDown;             // decide whether allows <delta E>d to be smaller than grain size.
    bool   print_TabbedMatrices;           // print tabbed instead of fixed-widthed matrices.
    int    showCollisionOperator;          // Show collision operator before and after normalization for each well.
                                           // 1: after normalization, 2: after deducting with I
                                           // 0 or anything else: before normalization
    bool   useDOSweightedDT;               // Use number of states to weigh the downward transition in collisionOperator()
    std::string punchSymbols;              // a string holds the symbols of rates.
    std::string punchNumbers;
    std::string punchFileName;
    bool   overwriteXmlAnalysis;          // Set when fitting so that all the intermediate results are not output to XML
    bool   autoSetMaxEne;                 // Check equilibrium population of all species is less than popThreshold at energy cut-off. 
	double popThreshold;                  // Equilibrium population criteria used in cutt-off check.
    bool   bForceMacroDetailedBalance;    // Impose detailed balance at the macroscopic level for non-conservative systems.
  };
}//namespace


#endif // GUARD_MesmerFlags_h

