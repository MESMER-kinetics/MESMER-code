#include "MesmerEnv.h"

namespace mesmer{
  MesmerEnv::MesmerEnv()
    :beta(1.0/200.), // Initialized to 300 K expressed in cm-1.
    conc(0.0),
    bathGasName(),
    GrainSize(0),
    CellSize(1.0),
    MaxGrn(0),
    MaxCell(10000L), // Initialized to allow parsing of <me:Hf298>. Overridden for real calculation.
    MaximumTemperature(0.0),
    EMin(0.0),
    EMax(0.0),
    EAboveHill(20.),
    useBasisSetMethod(false),
    nBasisSet(2){}
}//namespace

