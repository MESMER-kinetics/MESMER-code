#ifndef GUARD_MesmerConfig_h
#define GUARD_MesmerConfig_h

// Compiler specific configuration - all in this file

#ifdef __CYGWIN__
#define IsNan isnan
#endif

#if defined(WIN32)
#define IsNan _isnan
#include <conio.h>

#endif

#endif // GUARD_MesmerConfig_h
