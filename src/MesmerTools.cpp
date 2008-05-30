//MesmerTools.cpp
#include "MesmerTools.h"


namespace mesmer
{

  void Beyer_Swinehart(const vector<double>& VibFreq, vector<double>& cellDOS){
    // Implementation of the Beyer-Swinehart algorithm.
    const int MaximumCell = int(cellDOS.size());
    for ( vector<double>::size_type j = 0 ; j < VibFreq.size() ; ++j ) {
      int freq = static_cast<int>(VibFreq[j]) ;
      for (int i = 0 ; i < MaximumCell - freq ; ++i ){
        cellDOS[i + freq] += cellDOS[i] ;
      }
    }
  }

  // translation contribution for the partition function of two molecules
  long double translationalContribution(const long double m1, const long double m2, const long double beta){
    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tp_C = 2.0593e19 * pow(2. * M_PI ,1.5);

    return (tp_C * pow(m1 * m2 / ((m1 + m2) * beta), 1.5l));
  }

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta){
    double CanPrtnFn(0.);
    int vsize = static_cast<int>(DOS.size());
    for (int i = 0; i < vsize; ++i) {
      if (DOS[i] > 0.0)
        CanPrtnFn += exp( log(DOS[i]) - beta*Ene[i] ) ;
    }
    return CanPrtnFn;
  }

  //// convolutes v2 into v1
  //void convolution(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& sum){
  //  if (sum.size())
  //    sum.clear();
  //
  //  int MaximumCell = int(cellDOS.size());
  //  for(int i = 0; i < MaximumCell; ++i) {
  //    sum.push_back(0.0);
  //    for(int j = 0; j <= i; ++j){
  //      sum[i] += v1[i - j] * v2[j];
  //    }
  //  }
  //}

}
