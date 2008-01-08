#include "MesmerMath.h"
#include "TimeCounter.h"

using namespace mesmer;
//convolutes DOSs
void DOSconvolution(const std::vector<double> &f1,
                    const std::vector<double> &f2,
                    std::vector<double> &conv)
{
  int vSize = f1.size();
  for (int i = 0; i < vSize; ++i){
    conv[i] = 0.;
    for (int j = 0; j <= i; ++j){
      conv[i] += f1[j] * f2[i - j];
    }
  }
}
