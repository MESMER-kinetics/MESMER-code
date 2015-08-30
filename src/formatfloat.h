#ifndef GUARD_formatfloat_h
#define GUARD_formatfloat_h

#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
using namespace std;

namespace mesmer
{
  template<typename T>
  void formatFloat( ostream&      out,       // Output Stream.
                    const T&      datum,     // Floating point value.
                    const int     precision, // Required precision.
                    const int     width      // Required width.
  ){

    ostringstream sstrdatum ;

    sstrdatum << setprecision(precision) << datum ;

    string thisString = sstrdatum.str();

    sstrdatum.clear();

#ifdef USE_OLD_EXPONENT //(see in main())
    // insert a '0' if there is only one zero after the +/- sign
    string::size_type position = thisString.find('e', 0);
    if (position != string::npos){
      string sub = thisString.substr(position + 2);
      if (sub.length() == 2){
        thisString.insert(position + 2, 1, '0');
      }
    }
#endif

    ostringstream sstrdatum1 ;
    sstrdatum1 << thisString;

    out.setf(ios::right, ios::adjustfield) ;

    out << setw(width) << sstrdatum1.str() ;

  }

  template<typename T>
  string formatFloat(const T&      datum,     // Floating point value.
                     const int     precision, // Required precision.
                     const int     width      // Required width.
  ){

    ostringstream sstrdatum ;

    sstrdatum << setprecision(precision) << datum ;

    string thisString = sstrdatum.str();

    sstrdatum.clear();

#ifdef USE_OLD_EXPONENT
    // insert a '0' if there is only one zero after the +/- sign
    string::size_type position = thisString.find('e', 0);
    if (position != string::npos){
      string sub = thisString.substr(position + 2);
      if (sub.length() == 2){
        thisString.insert(position + 2, 1, '0');
      }
    }
#endif

    ostringstream sstrdatum1 ;
    sstrdatum1.setf(ios::right, ios::adjustfield) ;

    sstrdatum1 << setw(width) << thisString ;
    
    return sstrdatum1.str();

  }

}
#endif //GUARD_formatfloat_h
