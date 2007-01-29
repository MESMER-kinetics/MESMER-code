#include <sstream>
#include <iomanip>
#include <time.h>

#include "Persistence.h"

namespace mesmer
{
//
// Format floating point datum.
// This method is required as the stream manipulators do not appear to work
// for floating point variables, at least under linux (RH 7.0). SHR 16/Mar/2003.
//
void formatFloat( std::ostream&     out,       // Output Stream.
                                  const double datum,     // Floating point value.
                                  const int    precision, // Required precision.
                                  const int    width      // Required width.
){
  std::ostringstream sstrdatum ;

  sstrdatum << std::setprecision(precision) << datum ;

  out.setf(std::ios::right, std::ios::adjustfield) ;

    out << std::setw(width) << sstrdatum.str().c_str() ;
}

}//namespace mesmer