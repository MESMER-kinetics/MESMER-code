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
    cellFluxEnabled(false),
    grainFluxEnabled(false),
    rateCoefficientsOnly(false),
    useTheSameCellNumber(false),
    grainedProfileEnabled(false),
    speciesProfileEnabled(false),
    viewEvents(false),
    maxEvolutionTime(0.),
    printEigenValuesNum(-1),
    printReactionOperatorNum(0),
    allowSmallerDEDown(false),
    useFFTIntegration(true),
    print_TabbedMatrices(false)
    {}
}//namespace

