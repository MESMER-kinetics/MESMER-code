#include <sstream>
#include <iomanip>
#include "formatfloat.h"
using namespace std;
namespace mesmer
{
//
// Format floating point datum.
// This method is required as the stream manipulators do not appear to work
// for floating point variables, at least under linux (RH 7.0). SHR 16/Mar/2003.
//
void formatFloat( ostream&      out,       // Output Stream.
                           const double& datum,     // Floating point value.
                           const int     precision, // Required precision.
                           const int     width      // Required width.
){

  ostringstream sstrdatum ;

     sstrdatum << setprecision(precision) << datum ;

     out.setf(ios::right, ios::adjustfield) ;

     out << setw(width) << sstrdatum.str().c_str() ;

}

}
