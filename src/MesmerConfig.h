#ifndef GUARD_MesmerConfig_h
#define GUARD_MesmerConfig_h

#include <float.h>
#include <string>

// -------------------   Compiler specific configuration


#if defined (CYGWIN) // Platform definition 
#define IsNan isnan
#elif defined (UNIX)
#define IsNan isnan
#elif defined (LINUX)
#define IsNan isnan
#elif defined (HPUX)
#define IsNan isnan

#elif defined (_MSC_VER)
#define IsNan _isnan
#include <conio.h>
//Ignore warning C4800: 'int' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning( disable : 4800 )
const char FileSeparatorChar = '\\';
const std::string TestSubFolder("baselines\\Win32\\mesmer.test"); //relative to the data file
#else //suppose it is LINUX but not defined

#define IsNan isnan
const char FileSeparatorChar = '/';

#endif // Platform definition 

#endif // GUARD_MesmerConfig_h
