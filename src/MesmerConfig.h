#ifndef GUARD_MesmerConfig_h
#define GUARD_MesmerConfig_h

// -------------------   Compiler specific configuration
#if defined (CYGWIN) // Platform definition 
#define IsNan isnan
#elif defined (UNIX)
#define IsNan isnan
#elif defined (LINUX)
#define IsNan isnan
#elif defined (HPUX)
#define IsNan isnan
#elif defined (WIN32)
#define IsNan _isnan
#include <conio.h>
#else //suppose it is LINUX but not defined

#define IsNan isnan

#endif // Platform definition 

#endif // GUARD_MesmerConfig_h
