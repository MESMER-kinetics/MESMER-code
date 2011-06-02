#include "MesmerEnv.h"

namespace mesmer{
  MesmerEnv::MesmerEnv():
    beta(200.), // Initialized to 300 K expressed in cm-1.    conc(.0),
    GrainSize(100),
    MaxGrn(0),
    MaxCell(10000), //initialized to allow parsing of <me:Hf298>. Overridden for real calculation.
    MaximumTemperature(.0),
    EMin(.0),
    EMax(.0),
    EAboveHill(20.),
    useBasisSetMethod(false),
    nBasisSet(2){}
}//namespace

