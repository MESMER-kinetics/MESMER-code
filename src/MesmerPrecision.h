#ifndef GUARD_MesmerPrecision_h
#define GUARD_MesmerPrecision_h

//-------------------   Precision   ----------------
//uncomment one of these to choose precision package
#define USE_DOUBLE
//#define USE_DD
//#define USE_QD
//#define USE_CXSC
//--------------------------------------------------

#if defined USE_DD
#include <qd/dd_real.h>
#define MesmerPrecisionMethod dd_real
#define precisionTag 1000

#elif defined USE_QD
#include <qd/qd_real.h>
#define MesmerPrecisionMethod qd_real
#define precisionTag 2000


#elif defined USE_CXSC
#include <l_real.h>
#define MesmerPrecisionMethod l_real
using namespace cxsc;
#define precisionTag 3000

#else
#define MesmerPrecisionMethod double
#define precisionTag 4000
#endif

// ------------- Mesmer High precision ----------------
#ifndef MesmerHP
#define MesmerHP MesmerPrecisionMethod
#endif //MesmerHP
// ----------------------------------------------------


#endif // GUARD_MesmerPrecision_h

