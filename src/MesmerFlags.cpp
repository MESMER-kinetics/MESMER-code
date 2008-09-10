#include "MesmerFlags.h"

namespace mesmer{
  MesmerFlags::MesmerFlags():
    searchMethod(0),
    testDOSEnabled(false),
    testRateConstantEnabled(false),
    microRateEnabled(false),
    grainBoltzmannEnabled(false),
    grainDOSEnabled(false),
    cyclePrintGrainDOS(false),
    cellDOSEnabled(false),
    cyclePrintCellDOS(false),
    reactionOCSEnabled(false),
    kfEGrainsEnabled(false),
    kbEGrainsEnabled(false),
    TunnellingCoeffEnabled(false),
    cellTSFluxEnabled(false),
    grainTSFluxEnabled(false),
    rateCoefficientsOnly(false),
    useTheSameCellNumber(false),
    grainedProfileEnabled(false),
    speciesProfileEnabled(false),
    viewEvents(false),
    maxEvolutionTime(0.),
    printEigenValuesNum(-1),
    printReactionOperatorNum(0)
    {}
}//namespace

