#ifndef GUARD_MesmerPrecision_h
#define GUARD_MesmerPrecision_h

#include "MesmerConfig.h"
//#include <qd/dd_real.h>
//#include <qd/qd_real.h>

//#include <qd/dd.h>  // for chmlin14
//#include <qd/qd.h>  // for chmlin14
//#include <qd/fpu.h> // for chmlin14

// This is a dummy conversion for double to double itself (just to escape the compiler error). 
// QD has its own to_double() function converting QD to double or DD to double.
// See include files such as qd_inline.h in qd include folder.
inline double to_double(const double &a) {
  return a;
}

#endif // GUARD_MesmerPrecision_h

