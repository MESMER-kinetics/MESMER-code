//MesmerTools.cpp
#include "MesmerTools.h"

using namespace Constants;

namespace mesmer
{

  // translation contribution for the partition function of two molecules
  double translationalContribution(const double m1, const double m2, const double beta) {
    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tp_C = 2.0593e19 * pow(2. * M_PI ,1.5);

    return (tp_C * pow(m1 * m2 / ((m1 + m2) * beta), 1.5));
  }

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta) {

    double CanPrtnFn(0.0);
    for (size_t i(0), j(DOS.size() - 1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0)
        CanPrtnFn += exp(log(DOS[j]) - beta * Ene[j]);
    }
    return CanPrtnFn;

  }

  // Calculate the mean and variance of a boltzmann distribution.
  void canonicalMeanEnergy(const vector<double>& DOS, const vector<double>& Ene, const double beta, double& meanEnergy, double& varEnergy) {

    double meanEnergySqr(0.0), CanPrtnFn(0.0);
    meanEnergy = 0.0;
    for (size_t i(0), j(DOS.size() - 1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0) {
        double ene = Ene[j];
        double tmp = exp(log(DOS[j]) - beta * ene);
        CanPrtnFn += tmp;
        meanEnergy += ene * tmp;
        meanEnergySqr += ene * ene*tmp;
      }
    }

    meanEnergy /= CanPrtnFn;
    varEnergy = meanEnergySqr / CanPrtnFn - meanEnergy * meanEnergy;

    return;
  }

  //
  // Calculate the average grain energy and then number of states per grain.
  //
  void calcGrainAverages(const size_t &MaximumGrain, const size_t &cellPerGrain, const size_t &cellOffset, const vector<double>& CellDOS, const vector<double>& CellEne, vector<double>& grainDOS, vector<double>& grainEne, string name)
  {
    grainEne.clear();
    grainDOS.clear();
    grainEne.resize(MaximumGrain, 0.);
    grainDOS.resize(MaximumGrain, 0.);

    // Check that there are enough cells.
    //if (GrainSize < 1) {
    //  throw (std::runtime_error("The number of Cells is insufficient to produce requested number of Grains.")); 
    //}

    size_t idx1 = 0;
    size_t idx2 = 0;
    for (size_t i(0); i < MaximumGrain; ++i) {

      size_t idx3(idx1);

      // Account for the cell off set against PES grid by altering
      // the range of first grain average.

      const size_t cellRange = (i == 0) ? cellPerGrain - cellOffset : cellPerGrain;

      // Calculate the number of states in a grain.
      double gNOS = 0.0;
      for (size_t j(0); j < cellRange; ++j, ++idx1) {
        gNOS += CellDOS[idx1];
      }

      // Calculate average energy of the grain if it contains sum states.
      if (gNOS > 0.0) {
        double gSE = 0.0; // grain sum of state energy
        for (size_t j(0); j < cellRange; ++j, ++idx3) {
          gSE += CellEne[idx3] * CellDOS[idx3];
        }
        grainDOS[idx2] = gNOS;
        grainEne[idx2] = gSE / gNOS;
        idx2++;
      }

    }

    // Issue warning if number of grains produced is less that requested.

    if (idx2 != MaximumGrain) {
      cinfo << "Number of grains produced is not equal to that requested for: " << name << once << endl
        << "Number of grains requested: " << MaximumGrain << once << endl
        << "Number of grains produced : " << idx2 << once << endl;
    }
  }
}
