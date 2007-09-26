#ifndef GUARD_MesmerConfig_h
#define GUARD_MesmerConfig_h

// -------------------   Compiler specific configuration
#ifdef __CYGWIN__
#define IsNan isnan
#else

#ifdef __UNIX__
#define IsNan isnan
#else

#ifdef __LINUX__
#define IsNan isnan
#else

#ifdef __HPUX__
#define IsNan isnan
#else

#ifdef __WIN32__
#define IsNan _isnan
#include <conio.h>
#else //suppose it is LINUX but not defined

#define IsNan isnan

#endif //__WIN32__
#endif //__HPUX__
#endif //__LINUX__
#endif //__UNIX__
#endif //__CYGWIN__

// -------------------   Precision

// Mesmer High precision
#ifndef MesmerHP

#define MesmerHP long double
#endif //MesmerHP

#endif // GUARD_MesmerConfig_h
