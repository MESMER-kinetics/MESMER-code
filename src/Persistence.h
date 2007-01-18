#ifndef GUARD_Persistence_h
#define GUARD_Persistence_h

//-------------------------------------------------------------------------------------------
//
// Persistence.h 
//
// Author: Struan Robertson 
// Date:   30/Dec/2006
//
// This header file contains the declaration of the persistence interface.
//
//-------------------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include "tinyxml.h"

namespace mesmer
{

class IPersistObject {

public:

	IPersistObject(){} ;
	virtual ~IPersistObject(){} ;

	virtual bool ReadFromXML(TiXmlElement* pnMol) = 0 ;

	// Utillity methods:

	// Returns the effective content of a CML property element
  static const char* ReadProperty(TiXmlElement* pnList, const std::string& name, bool MustBeThere=true);

  static void formatFloat(std::ostream& out, const double datum, const int precision, const int width);

  /// Inserts into XML document a new element /<name/>
  TiXmlElement* WriteElement(TiXmlElement* parentElement, const std::string name);

  /// Inserts into XML document a new element /<name/> containing a formatted number
  static void WriteValueElement( TiXmlElement* parentElement, const std::string name, 
                                 const double datum, const int precision=-1);

  //Insert into XML document a new element /<name/> and gives it a timestamp attribute and comment 
  static TiXmlElement* WriteMainElement(TiXmlElement* parentElement, const std::string name,
                                  const std::string comment, bool replaceExisting=true);
} ;

}

#endif // GUARD_Persistence_h

