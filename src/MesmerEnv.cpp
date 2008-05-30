#include "MesmerEnv.h"

namespace mesmer{
  MesmerEnv::MesmerEnv():
    beta(1.), // not sure what value to initialize. CHL
    conc(.0),
    GrainSize(100),
    MaxGrn(0),
    MaxCell(0),
    MaximumTemperature(.0),
    EMin(.0),
    EMax(.0),
    testDOSEnabled(false),
    microRateEnabled(false),
    grainBoltzmannEnabled(false),
    grainDOSEnabled(false),
    cellDOSEnabled(false),
    collisionOCSEnabled(false),
    kfEGrainsEnabled(false),
    kbEGrainsEnabled(false),
    TunnellingCoeffEnabled(false),
    cellTSFluxEnabled(false),
    grainTSFluxEnabled(false),
    rateCoefficientsOnly(false),
    useTheSameCellNumber(false),
    grainedProfileEnabled(false),
    EAboveHill(20.),
    maxEvolutionTime(0.),
    printEigenValuesNum(-1)
    {}
}//namespace

