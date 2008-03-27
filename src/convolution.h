//convolution.h
#include <vector>

void Beyer_Swinehart(const std::vector<double>& VibFreq, std::vector<double>& cellDOS){
  // Implementation of the Beyer-Swinehart algorithm.
  int MaxCell = int(cellDOS.size());
  for ( vector<double>::size_type j = 0 ; j < VibFreq.size() ; ++j ) {
    int freq = static_cast<int>(VibFreq[j]) ;
    for (int i = 0 ; i < MaxCell - freq ; ++i ){
      cellDOS[i + freq] += cellDOS[i] ;
    }
  }
}

//// convolutes v2 into v1
//void convolution(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& sum){
//  if (sum.size())
//    sum.clear();
//
//  int MaxCell = int(cellDOS.size());
//  for(int i = 0; i < MaxCell; ++i) {
//    sum.push_back(0.0);
//    for(int j = 0; j <= i; ++j){
//      sum[i] += v1[i - j] * v2[j];
//    }
//  }
//}
