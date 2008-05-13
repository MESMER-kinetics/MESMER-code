#ifndef GUARD_MesmerPrecision_h
#define GUARD_MesmerPrecision_h

#include "MesmerConfig.h"
//#include <qd/dd_real.h>
//#include <qd/qd_real.h>


//-------------------   Precision   ----------------
//uncomment one of these to choose precision package
#define USE_DOUBLE
//#define USE_DD
//#define USE_QD
//--------------------------------------------------

//---------------------
//variable definition
#define varUseDouble              1000
#define varUseDoubleDouble        2000
#define varUseQuadDouble          3000
//---------------------


#if defined USE_DD
#include <qd/dd_real.h>
#define MesmerPrecisionMethod dd_real
#define precisionTag varUseDoubleDouble

#elif defined USE_QD
#include <qd/qd_real.h>
#define MesmerPrecisionMethod qd_real
#define precisionTag varUseQuadDouble

#else
#define MesmerPrecisionMethod double
#define precisionTag varUseDouble
#endif

// ------------- Mesmer High precision ----------------
#ifndef MesmerHP
#define MesmerHP MesmerPrecisionMethod
#endif //MesmerHP
// ----------------------------------------------------

// This is a dummy conversion for double to double itself (just to escape the compiler error). 
// QD has its own to_double() function converting QD to double or DD to double.
// See include files such as qd_inline.h in qd include folder.
inline double to_double(const double &a) {
  return a;
}

#endif // GUARD_MesmerPrecision_h

