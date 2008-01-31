//Bayer_Swinehart.h
#include <vector>

void Bayer_Swinehart(const std::vector<double>& VibFreq, std::vector<double>& cellDOS, const int MaxCell){
  // Implementation of the Bayer-Swinehart algorithm.
  for ( vector<double>::size_type j = 0 ; j < VibFreq.size() ; ++j ) {
    int iFreq = static_cast<int>(VibFreq[j]) ;
    for (int i = 0 ; i < MaxCell - iFreq ; ++i ){
      cellDOS[i + iFreq] += cellDOS[i] ;
    }
  }
}
