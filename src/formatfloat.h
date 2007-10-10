#ifndef GUARD_formatfloat_h
#define GUARD_formatfloat_h

#include <iostream>
#include <cmath>

namespace mesmer
{
  void formatFloat(std::ostream& out, const double& datum, const int precision, const int width);
}
#endif //GUARD_formatfloat_h
