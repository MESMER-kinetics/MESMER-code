#ifndef GUARD_ParseForPlugin_h
#define GUARD_ParseForPlugin_h

#include "plugin.h"

// This function is intended to be used to set up all types of plugin, both
// parsing the XML input file for a specified plugin type, making an instance
// of the requested plugin from its name and parsing for its data using the
// plugin's virtual function. Handles both old and new formats.
// Returns a pointer to the plugin or NULL if a valid plugin or if required
// data is not found.
// The parameter pParent is Molecule*, Reaction* or System* depending on the
// plugin.

template <class T,class U>
T* ParseForPlugin(U* pParent, const char* elTypeName, PersistPtr ppTop, bool MustBeThere=true)
{
  // Parse XML input for a specified plugin (with me: prefix), either in a text element
  // which is a child of ppTop e.g.<me:MCRCMethod>SimpleRRKM</me:MCRCMethod>
  // OR in a name or a xsi:type attribute e.g <me:MCRCMethod xsi:type ="me:MesmerILT">
  // If the xsi:type form is found, ppTop points to the element name on return.
  // Call this function with an explicit template parameter 
  const char* ptxt = NULL;
  PersistPtr pp = ppTop->XmlMoveTo(elTypeName);
  if(pp)
    ptxt = pp->XmlReadValue("xsi:type",optional);
  if(!ptxt)
  {
    if(pp)
      ptxt = pp->XmlReadValue("name",optional);
    if(!ptxt)
    {
      //Plugin name in text element or use value from defaults.xml
      ptxt = ppTop->XmlReadValue(elTypeName, MustBeThere) ; //older style
      if(ptxt && !*ptxt) //newer style
        ptxt = (ppTop->XmlMoveTo(elTypeName))->XmlReadValue("name",optional);
    }
    // With older formats, the data may be in sibling elements   
    pp = ppTop;
  }
  if(!ptxt)
    return NULL; //Element is not present or unlikely failure
  if(ptxt[2]==':')
    ptxt+=3; //Remove prefix "me:"

  //Look up plugin name
  T* pPlugin = T::Find(ptxt);

  //Read data if present
  bool dataOK;
  if(pPlugin)
  {
    pPlugin->setParent(pParent);
    dataOK = pPlugin->ParseData(pp);//use plugin's virtual function
    if(!dataOK)
    {
      cerr << "Data not found for plugin " << pPlugin->getID() << endl;
      return NULL;
    }
  }
  return pPlugin;
}

#endif

