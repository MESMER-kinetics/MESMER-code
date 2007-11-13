#include "MesmerEnv.h"

namespace mesmer{
  MesmerEnv::MesmerEnv():
    beta(1.), // not sure what value to initialize. CHL
    conc(.0),
    GrainSize(100.),
    MaxGrn(0),
    MaxCell(0),
    MaxT(.0),
    EMin(.0),
    EMax(.0),
    microRateEnabled(false){}
}//namespace

