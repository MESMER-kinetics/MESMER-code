#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h

#include <stddef.h>
#include <string>

namespace mesmer
{
  struct MesmerEnv
  {
    MesmerEnv();

    // environmental values
    double beta;
    double conc;
    std::string bathGasName;

    // granularity of the system
    size_t GrainSize;           // Grain size in cm-1
	int    stochTrials;         // Number of Stochatic simulations.
	double stochEndTime;        // End time for stochastic trajectory.
	double stochStartTime;      // Start time for stochastic trajectory.
	int    stochThermThresh;    //Threshold for thermalisation in stochastic simulation
	int    stochEquilThresh;    //Threshold for species to species equillibration in stochastic simulation
	int    stochAXDlimit;       // Distance below lowest barrier for top AXD boundary. 
    double CellSize;            // Cell size in cm-1
    size_t MaxGrn;              // The number of grains
    size_t MaxCell;             // The number of cells
    double MaximumTemperature;  // Maximum temperature for the purposes of setting the energy range
    double EMin, EMax;          // The absolute lowest and highest energies in the system, cm-1
    double EAboveHill;          // Max energy above the highest Hill [in kT]
    bool   useBasisSetMethod;   // Use the contracted basis set method.
    size_t nBasisSet;           // Number of basis set functions to use.



public:

    const size_t cellPerGrain() const { return size_t(double(GrainSize)/CellSize) ; } ;
  };

}//namespace


#endif // GUARD_MesmerEnv_h

