//MesmerTools.h
#ifndef GUARD_MesmerTools_h
#define GUARD_MesmerTools_h

#include <vector>
#include <stdexcept>
#include "MesmerMath.h"
#include "TimeCounter.h"
#include "unitsConversion.h"
#include "dMatrix.h"
#include "marray.h"
using namespace Constants;
using namespace std;

namespace mesmer
{
  void Beyer_Swinehart(const vector<double>& VibFreq, vector<double>& cellDOS);
  
  void reverseBeyer_Swinehart(double VibFreq, vector<double>& cellDOS);

  // translation contribution for the partition function of two molecules
  double translationalContribution(const double m1, const double m2, const double beta);

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta);

  // shift cell DOS and energy vectors according to cellOffset
  void shiftCells(int MaximumCell, int cellOffset, const vector<double>& cellDOS, const vector<double>& cellEne,
    std::vector<double>& shiftedCellDOS, std::vector<double>& shiftedCellEne);

  // Calculate the average grain energy and then number of states per grain.
  void calcGrainAverages(const int MaximumGrain, const int GrainSize, const std::vector<double>& shiftedCellDOS,
    const std::vector<double>& shiftedCellEne, vector<double>& grainDOS, vector<double>& grainEne) ;

}

#endif // GUARD_MesmerTools_h
