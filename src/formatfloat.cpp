#include <sstream>
#include <iomanip>
#include "formatfloat.h"

namespace mesmer
{
  //
  // Format floating point datum.
  // This method is required as the stream manipulators do not appear to work
  // for floating point variables, at least under linux (RH 7.0). SHR 16/Mar/2003.
  //
//  void formatFloat( ostream&      out,       // Output Stream.
//                    const double& datum,     // Floating point value.
//                    const int     precision, // Required precision.
//                    const int     width      // Required width.
//  ){
//
//    ostringstream sstrdatum ;
//
//    sstrdatum << setprecision(precision) << datum ;
//
//    string thisString = sstrdatum.str();
//
//    sstrdatum.clear();
//
//    // insert a '0' if there is only one zero after the +/- sign
//    string::size_type position = thisString.find('e', 0);
//    if (position != string::npos){
//      string sub = thisString.substr(position + 2);
//      if (sub.length() == 2){
//        thisString.insert(position + 2, 1, '0');
//      }
//    }
//
//    ostringstream sstrdatum1 ;
//    sstrdatum1 << thisString;
//
//    out.setf(ios::right, ios::adjustfield) ;
//
//    out << setw(width) << sstrdatum1.str() ;
//
//  }

}
