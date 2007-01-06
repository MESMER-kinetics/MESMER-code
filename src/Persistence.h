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
   const char* ReadProperty(TiXmlElement* pnList, const std::string& name, bool MustBeThere=true);

   static void formatFloat(std::ostream& out, const double datum, const int precision, const int width);

} ;

}

#endif // GUARD_Persistence_h

