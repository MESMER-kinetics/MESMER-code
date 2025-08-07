#ifndef GUARD_MesmerConfig_h
#define GUARD_MesmerConfig_h

#include <float.h>
#include <string>
#include <limits>

// -------------------   Compiler specific configuration
#if defined (__CYGWIN__)      // CYGWIN

#define IsNan isnan
#define _isnan isnan

static const double NaN = std::numeric_limits<double>::quiet_NaN();


#elif defined (UNIX)      // UNIX

#define IsNan isnan
static const double NaN = NAN;

#elif defined (LINUX)     // LINUX

#include <math.h>

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
// #define isfinite _finite

#else                     //suppose it is LINUX but not defined

#define IsNan isnan
const char FileSeparatorChar = '/';
static const double NaN = NaN;

#endif
// -------------------   Compiler specific configuration


#if defined (_MSC_VER)  // MS Windows
const std::string TestSubFolder("baselines\\Win32\\mesmer.test"); //relative to the data file
const std::string TestFolder("baselines\\Win32\\"); //relative to the data file
//#elif defined (SOMETHING)
//const std::string TestSubFolder("baselines/Something/mesmer.test"); //relative to the data file
//const std::string TestFolder("baselines/Something/"); //relative to the data file
#else
const std::string TestSubFolder("baselines/Linux64/mesmer.test"); //relative to the data file
const std::string TestFolder("baselines/Linux64/"); //relative to the data file
#endif

#endif // GUARD_MesmerConfig_h
