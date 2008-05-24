//MesmerTools.h
#ifndef GUARD_MesmerTools_h
#define GUARD_MesmerTools_h

#include <vector>
#include "MesmerMath.h"
#include "TimeCounter.h"
#include "unitsConversion.h"
#include "dMatrix.h"
using namespace Constants;
using namespace std;

namespace mesmer
{
  void Beyer_Swinehart(const vector<double>& VibFreq, vector<double>& cellDOS);

  // translation contribution for the partition function of two molecules
  long double translationalContribution(const long double m1, const long double m2, const long double beta);

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta);

  //// convolutes v2 into v1
  //void convolution(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& sum);
}

#endif // GUARD_MesmerTools_h
