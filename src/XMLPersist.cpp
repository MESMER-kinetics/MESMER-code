#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include "XMLPersist.h"

using namespace std;

namespace mesmer
{
PersistPtr XMLPersist::XmlCreate(const std::string& inputfilename, const std::string& title)
{
  TiXmlDocument* pdoc = new TiXmlDocument( inputfilename.c_str() );//Deleted in destructor
  if( !pdoc->LoadFile() )
  {
    stringstream errorMsg;
    errorMsg << "Could not load file " << inputfilename
             << "\nIt may not exist or it may not consist of well-formed XML." << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);

    delete pdoc;
    return PersistPtr(NULL);
  }

  TiXmlElement* root = pdoc->RootElement();
  if(title.size() && root->ValueStr()!=title)
  {
    stringstream errorMsg;
    errorMsg << inputfilename << " does not have a root element or title named " << title << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    delete pdoc;
    return PersistPtr(NULL);
  }
  return PersistPtr(new XMLPersist(root, pdoc));//this XMLPersist object has a document pointer
}

XMLPersist::~XMLPersist()
{
  delete pDocument; //doesn't matter that pDocument is usually NULL
}

PersistPtr XMLPersist::XmlMoveTo(const std::string& name) const
{
    TiXmlElement* pnEl = pnNode->FirstChildElement(name);
    if(!pnEl)
      pnEl = pnNode->NextSiblingElement(name);
    return PersistPtr(new XMLPersist(pnEl));
}

const char* XMLPersist::XmlRead()const
{
  return pnNode->GetText();
}

///Look first to see if there is a child element of this name.
///Look second for an attribute of this name.
///If either found, return its value.
///Otherwise return NULL. If MustBeThere is true(the default) also give an error message.
const char* XMLPersist::XmlReadValue(const std::string& name, bool MustBeThere) const
{
  const char* ptext=NULL;
  //Look first to see if there is a child element of this name and, if so, return its value
  TiXmlElement* pnEl = pnNode->FirstChildElement(name);
  if(pnEl)
  {
    ptext = pnEl->GetText();
    if(!ptext)
      return ""; //element exists but is empty
  }
  else
    ptext = pnNode->Attribute(name.c_str());

  if(!ptext && MustBeThere){
    stringstream errorMsg;
    errorMsg << "The " << name << " element or attribute was missing or empty.";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
  }
  return ptext;
}

/// Returns the effective content of an CML <property> element
/// Looks for child elements pnList of the form:
/// <property dictRef="name">
///     <scalar> content </scalar>
/// </property>
/// The property can have <array>, <string>, or anything, in place of <scalar>
/// Returns NULL if the appropriate property is not found or if it has no content.

const char* XMLPersist::XmlReadProperty(const string& name, bool MustBeThere) const
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
  if(MustBeThere){
    stringstream errorMsg;
    errorMsg << "Ill-formed " << name << " in the molecule ";//calling function to add name of molecule
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
  }
  return NULL;
}

  /// Returns true if datatext associated with name is "1" or "true" or nothing;
  //  returns false if datatext is something else or if element is not found.
  bool XMLPersist::XmlReadBoolean( const std::string& name)const
  {
    const char* txt = XmlReadValue(name, false);
    if(txt)
    {
      string s(txt);
      return s.empty() || s=="1" || s=="yes";
    }
    return false;
  }


void XMLPersist::XmlWriteValueElement(const std::string& name,
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

PersistPtr XMLPersist::XmlWriteElement(const std::string& name)
{
  TiXmlElement* item = new TiXmlElement( name );
  pnNode->LinkEndChild(item);
  return PersistPtr(new XMLPersist(item));
}

void XMLPersist::XmlWriteAttribute(const std::string& name, const std::string& value)
{
  pnNode->SetAttribute(name, value);
}

PersistPtr XMLPersist::XmlWriteMainElement(
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

  //No timestamp or comment if comment is empty
  if(comment.size())
  {
    //----------------------------------------
    TimeCount events;
    std::string thisEvent, timeString;
    thisEvent = "Record calculated data";
    timeString = events.setTimeStamp(thisEvent);
    cout << thisEvent << " at " << timeString << endl;
    pnel->SetAttribute("calculated", timeString);
    //----------------------------------------

    //Add explanatory comment in description element
    TiXmlElement* pncom = new TiXmlElement("me:description");
    pnel->LinkEndChild( pncom );
    pncom->LinkEndChild(new TiXmlText(comment));
  }
  return PersistPtr(new XMLPersist(pnel));
}

bool XMLPersist::XmlSaveFile(const std::string& outfilename)
  {
    return pnNode->GetDocument()->SaveFile(outfilename);
  }

}//namespacer mesmer
