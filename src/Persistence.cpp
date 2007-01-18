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
#include <time.h>
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

void IPersistObject::WriteValueElement( TiXmlElement* parentElement, const std::string name, 
                 const double datum, const int precision)
{
  ostringstream sstrdatum ;
  if(precision>=0)
    sstrdatum << setprecision(precision) << datum ;
  else
    sstrdatum << datum ;

  
  TiXmlElement* item = new TiXmlElement(name);
  parentElement->LinkEndChild(item);
  item->LinkEndChild(new TiXmlText(sstrdatum.str()));
}

TiXmlElement* IPersistObject::WriteElement(TiXmlElement* parentElement, const std::string name)
{
  TiXmlElement* item = new TiXmlElement( name );
  parentElement->LinkEndChild(item);
  return item;
}

TiXmlElement* IPersistObject::WriteMainElement(TiXmlElement* parentElement, 
                        const std::string name, const std::string comment, bool replaceExisting)
{
  // Delete any existing element of the same name, unless explicitly asked not to.
  if(replaceExisting)
  {
    TiXmlNode* pnExisting = parentElement->FirstChild(name);
    if(pnExisting)
      parentElement->RemoveChild(pnExisting);
  }

  // Add new element with specified name
  TiXmlElement* pnel = new TiXmlElement( name );
  parentElement->LinkEndChild(pnel);

  // Add attribute to show when data was calculated, removing trailing 0x0a
  time_t ltime;
  time( &ltime );
  string timestring(ctime(&ltime));
  pnel->SetAttribute("calculated", timestring.erase(timestring.size()-1));

  //Add explanatory comment in description element
  TiXmlElement* pncom = new TiXmlElement("me:description");
  pnel->LinkEndChild( pncom );
  pncom->LinkEndChild(new TiXmlText(comment));

  return pnel;
}

} //namespace