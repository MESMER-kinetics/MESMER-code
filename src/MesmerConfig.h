#ifndef GUARD_MesmerConfig_h
#define GUARD_MesmerConfig_h

#include <float.h>
#include <string>


// -------------------   Compiler specific configuration
#if defined (CYGWIN)      // CYGWIN

#define IsNan isnan
static const double NaN = NAN;

#elif defined (UNIX)      // UNIX

#define IsNan isnan
static const double NaN = NAN;

#elif defined (LINUX)     // LINUX

#define IsNan isnan
static const double NaN = NAN;

#elif defined (HPUX)      // HPUX

#define IsNan isnan
static const double NaN = NAN;

#elif defined (_MSC_VER)  // MS Windows

#define IsNan _isnan
#include <conio.h>
//Ignore warning C4800: 'int' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning( disable : 4800 )
const char FileSeparatorChar = '\\';
static const double NaN = std::numeric_limits<double>::quiet_NaN();

#else                     //suppose it is LINUX but not defined

#define IsNan isnan
const char FileSeparatorChar = '/';
static const double NaN = NAN;

#endif
// -------------------   Compiler specific configuration



const std::string TestSubFolder("baselines\\Win32\\mesmer.test"); //relative to the data file

#endif // GUARD_MesmerConfig_h
