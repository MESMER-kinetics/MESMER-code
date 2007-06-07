#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include "XMLPersist.h"

using namespace std;

namespace mesmer
{
PersistPtr XMLPersist::Create(const std::string& inputfilename, const std::string& title)
{
  TiXmlDocument* pdoc = new TiXmlDocument( inputfilename.c_str() );//Deleted in destructor
  if( !pdoc->LoadFile() )
  {
    cerr << "Could not load file " << inputfilename 
         << "\nIt may not exist or it may not consist of well-formed XML." << endl;
    delete pdoc;
    return PersistPtr(NULL);
  }

  TiXmlElement* root = pdoc->RootElement();
  if(!title.empty() && root->ValueStr()!=title)
  {
    cerr << inputfilename << " does not have a root element or title named " << title << endl;
    delete pdoc;
    return PersistPtr(NULL);
  }
  return PersistPtr(new XMLPersist(root, pdoc));//this XMLPersist object has a document pointer
}

XMLPersist::~XMLPersist()
{
  delete pDocument; //doesn't matter that pDocument is usually NULL
}

PersistPtr XMLPersist::MoveTo(const std::string& name) const
{
    TiXmlElement* pnEl = pnNode->FirstChildElement(name);
    if(!pnEl)
      pnEl = pnNode->NextSiblingElement(name);
    return PersistPtr(new XMLPersist(pnEl));
}

///Look first to see if there is a child element of this name.
///Look second for an attribute of this name.
///If either found, return its value. 
///Otherwise return NULL. If MustBeThere is true(the default) also give an error message.
const char* XMLPersist::ReadValue(const std::string& name, bool MustBeThere) const
{
  const char* ptext=NULL;
  //Look first to see if there is a child element of this name and, if so, return its vale
  TiXmlElement* pnEl = pnNode->FirstChildElement(name);
  if(pnEl)
    ptext = pnEl->GetText();
  else
    ptext = pnNode->Attribute(name.c_str());

  if(!ptext && MustBeThere)
    cerr << "The " << name << " element or attribute was missing or empty." << endl;

  return ptext;
}

/// Returns the effective content of an CML <property> element
/// Looks for child elements pnList of the form: 
/// <property dictRef="name">
///     <scalar> content </scalar>
/// </property>
/// The property can have <array>, <string>, or anything, in place of <scalar>
/// Returns NULL if the appropriate property is not found or if it has no content.

const char* XMLPersist::ReadProperty(const string& name, bool MustBeThere) const
{
  TiXmlElement* pnProp = pnNode->FirstChildElement("property");
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
    cerr << "Ill-formed " << name << " in the molecule "; //calling function to add name of molecule
  return NULL;
}

void XMLPersist::WriteValueElement(const std::string& name, 
                 const double datum, const int precision)
{
  ostringstream sstrdatum ;
  if(precision>=0)
    sstrdatum << setprecision(precision) << datum ;
  else
    sstrdatum << datum ;

  
  TiXmlElement* item = new TiXmlElement(name);
  pnNode->LinkEndChild(item);
  item->LinkEndChild(new TiXmlText(sstrdatum.str()));
}

PersistPtr XMLPersist::WriteElement(const std::string& name)
{
  TiXmlElement* item = new TiXmlElement( name );
  pnNode->LinkEndChild(item);
  return PersistPtr(new XMLPersist(item));
}

PersistPtr XMLPersist::WriteMainElement( 
                        const std::string& name, const std::string& comment, bool replaceExisting)
{
  // Delete any existing element of the same name, unless explicitly asked not to.
  if(replaceExisting)
  {
    TiXmlNode* pnExisting = pnNode->FirstChild(name);
    if(pnExisting)
      pnNode->RemoveChild(pnExisting);
  }

  // Add new element with specified name
  TiXmlElement* pnel = new TiXmlElement( name );
  pnNode->LinkEndChild(pnel);

  // Add attribute to show when data was calculated, removing trailing 0x0a
  time_t ltime;
  time( &ltime );
  string timestring(ctime(&ltime));
  pnel->SetAttribute("calculated", timestring.erase(timestring.size()-1));

  //Add explanatory comment in description element
  TiXmlElement* pncom = new TiXmlElement("me:description");
  pnel->LinkEndChild( pncom );
  pncom->LinkEndChild(new TiXmlText(comment));

  return PersistPtr(new XMLPersist(pnel));
}

bool XMLPersist::SaveFile(const std::string& outfilename)
  {
    return pnNode->GetDocument()->SaveFile(outfilename);
  }
}//namespacer mesmer