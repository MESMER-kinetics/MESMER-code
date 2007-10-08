#ifndef GUARD_MesmerPrecision_h
#define GUARD_MesmerPrecision_h

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


inline double to_double(const double &a) {
  return a;
}


#endif // GUARD_MesmerPrecision_h

