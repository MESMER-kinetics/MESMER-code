#ifndef GUARD_MesmerEnv_h
#define GUARD_MesmerEnv_h

namespace mesmer
{
  class MesmerEnv
  {
  public:
    //
    // Declare System as a friend class. Not sure that this is the best or
    // most OOP way to go as it clearly defeats the point of encapsulation
    // but the System needs to know a lot about the reactions it is
    // combining. Review latter. SHR 2/Apr/2003.
    //
    friend class System ;
   
    MesmerEnv();

    double temp;
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

