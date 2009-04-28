#include "MesmerFlags.h"

namespace mesmer{
  MesmerFlags::MesmerFlags():
    searchMethod(UNDEFINED_METHOD),
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
    print_TabbedMatrices(true),
    showCollisionOperator(0),
    useDOSweighedDT(false),
    doBasisSetMethod(false),
    punchSymbols(),
    punchNumbers(),
    punchFileName()
    {}
}//namespace

