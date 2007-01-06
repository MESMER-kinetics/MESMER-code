//-------------------------------------------------------------------------------------------
//
// Persistence.cpp
//
// Author: Struan Robertson 
// Date:   30/Dec/2006
//
// This header file contains the definition of the persistence interface utility methods.
//
//-------------------------------------------------------------------------------------------

#include <sstream>
#include <iomanip>
#include "Persistence.h"

using namespace std ;

namespace mesmer
{

/// Returns the effective content of an CML <property> element
/// Looks for child elements pnList of the form: 
/// <property dictRef="name">
///     <scalar> content </scalar>
/// </property>
/// The property can have <array>, <string>, or anything, in place of <scalar>
/// Returns NULL if the appropriate property is not found or if it has no content.

const char* IPersistObject::ReadProperty(TiXmlElement* pnList, const string& name, bool MustBeThere)
{
  TiXmlElement* pnProp = pnList->FirstChildElement("property");
  while(pnProp)
  {
    const char* pAtt = pnProp->Attribute("dictRef");
    if(pAtt && name==pAtt)
    {
      TiXmlElement* pnChild = pnProp->FirstChildElement(); //could be <array> or <scalar> or <string>
      if(pnChild)
        return pnChild->GetText();
    }
    pnProp = pnProp->NextSiblingElement();
  }
  if(MustBeThere)
    cerr << "Ill-formed " << name << " in the molecule " << endl;
  return NULL;
}

//
// Format floating point datum.
// This method is required as the stream manipulators do not appear to work
// for floating point variables, at least under linux (RH 7.0). SHR 16/Mar/2003.
//
void IPersistObject::formatFloat( ostream&     out,       // Output Stream.
                                  const double datum,     // Floating point value.
                                  const int    precision, // Required precision.
                                  const int    width      // Required width.
){
	ostringstream sstrdatum ;

    sstrdatum << setprecision(precision) << datum ;

    out.setf(ios::right, ios::adjustfield) ;

    out << setw(width) << sstrdatum.str().c_str() ;

}

}