#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h

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
  };
}//namespace


#endif // GUARD_MesmerEnv_h

