//MesmerTools.h
#ifndef GUARD_MesmerTools_h
#define GUARD_MesmerTools_h

#include <vector>
#include <stdexcept>
#include "MesmerMath.h"
#include "unitsConversion.h"

template<typename T>
std::string toString(T t)
{
	std::ostringstream s; s << t; return s.str();
}

namespace mesmer
{
  // translation contribution for the partition function of two molecules
  double translationalContribution(const double m1, const double m2, const double beta);

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta);

  // Calculate the mean and variance of a boltzmann distribution.
  void canonicalMeanEnergy(const vector<double>& DOS, const vector<double>& Ene, const double beta, double& meanEnergy, double& varEnergy) ;

  // Calculate the average grain energy and then number of states per grain.
  void calcGrainAverages(const size_t &MaximumCell, const size_t &cellPerGrain, const size_t &cellOffset, const std::vector<double>& shiftedCellDOS,
    const std::vector<double>& shiftedCellEne, vector<double>& grainDOS, vector<double>& grainEne, std::string name = "Unknown") ;

  template<class T>
  void quickSort(vector<T> &v, vector<size_t> &index, size_t lo0, size_t hi0)
  {
    size_t lo = lo0;
    size_t hi = hi0;
    if (lo >= hi) return;
    T mid = v[(lo + hi) / 2];

    while (lo < hi)
    {
      while (lo<hi && v[lo] > mid) lo++;
      while (lo<hi && v[hi] < mid) hi--;
      if (lo < hi) {
        swap(v[lo],v[hi]) ;
        swap(index[lo],index[hi]) ;
      }
    }

    if (hi < lo){
      swap(v[lo],v[hi]) ;
      swap(index[lo],index[hi]) ;
    }

    quickSort(v, index, lo0, lo);
    quickSort(v, index, (lo == lo0) ? lo+1 : lo, hi0);
  }

}

#endif // GUARD_MesmerTools_h
