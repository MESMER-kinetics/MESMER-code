#include "MesmerMath.h"


//convolutes DOSs
void DOSconvolution(const std::vector<double> &f1,
                    const std::vector<double> &f2,
                    std::vector<double> &conv)
{
  for (unsigned int i = 0; i < f1.size(); ++i){
    conv[i] = 0.;
    for (unsigned int j = 0; j <= i; ++j)
      conv[i] += f1[j] * f2[i - j];
  }
}
