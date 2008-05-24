#include "MesmerMath.h"

//convolutes DOSs
void DOSconvolution(const std::vector<double> &f1,
                    const std::vector<double> &f2,
                    std::vector<double> &conv)
{
  std::vector<double>::size_type vSize = f1.size();
  for (std::vector<double>::size_type i = 0; i < vSize; ++i){
    conv[i] = 0.;
    for (std::vector<double>::size_type j = 0; j <= i; ++j){
      conv[i] += f1[j] * f2[i - j];
    }
  }
}
