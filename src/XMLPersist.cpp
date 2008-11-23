#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "XMLPersist.h"

using namespace std;

namespace mesmer
{
PersistPtr XMLPersist::XmlLoad(const std::string& inputfilename, const std::string& title)
{
  TiXmlDocument* pdoc = new TiXmlDocument();//Deleted in destructor

  bool ret;
  if(inputfilename.empty())
    ret = pdoc->LoadFile(stdin);
  else
    ret = pdoc->LoadFile(inputfilename);
  if( !ret )
  {
    cerr << "Could not load file " << inputfilename
             << "\nIt may not exist or it may not consist of well-formed XML." << endl;

    delete pdoc;
    return PersistPtr(NULL);
  }

  TiXmlElement* root = pdoc->RootElement();
  if(title.size() && root->ValueStr()!=title)
  {
    cerr << inputfilename << " does not have a required root element named " << title << endl;
    delete pdoc;
    return PersistPtr(NULL);
  }
  return PersistPtr(new XMLPersist(root, pdoc));//this XMLPersist object has a document pointer
}

XMLPersist::~XMLPersist()
{
  delete pDocument; //doesn't matter that pDocument is usually NULL
  //shall we delete pnNode as well?? CHL
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
    cinfo << "The " << name << " element or attribute was missing or empty." << endl;
  }
  return ptext;
}

/** 
Returns the effective content of an CML <property> element
Looks for child elements pnList of the form:
 <property dictRef="name">
   <scalar> content </scalar>
 </property>
 or alternatively:
 <property title="name">
   <scalar> content </scalar>
 </property>
 In the second case, only  the part of name after any colon (the "localname") is used
 The property can have <array>, <string>, or anything, in place of <scalar>
 Returns NULL if the appropriate property is not found or if it has no content.
 **/

const char* XMLPersist::XmlReadProperty(const string& name, bool MustBeThere) const
{
  TiXmlElement* pnProp = pnNode->FirstChildElement("property");
  while(pnProp)
  {
    size_t pos=0;
    const char* pAtt = pnProp->Attribute("dictRef");
    if(!pAtt)
    {
      pAtt = pnProp->Attribute("title");
      //The position of the start of the localname(the bit after the colon)
      pos=name.find(':') + 1; //ok even if there is no colon (if npos=-1)
    }
    if(pAtt && name.compare(pos, name.size()-pos, pAtt)==0)
    {
      TiXmlElement* pnChild = pnProp->FirstChildElement(); //could be <array> or <scalar> or <string>
      if(pnChild)
        return pnChild->GetText();
    }
    pnProp = pnProp->NextSiblingElement();
  }
//  if(MustBeThere)
//    meErrorLog.ThrowError(__FUNCTION__, "The property " + name + " is missing or empty", obError);
  return NULL;
}

/// Returns the attName attribute of an CML <property> element
/// See XMLPersist::XmlReadProperty for details
const char* XMLPersist::XmlReadPropertyAttribute(const string& name, const string& attName, bool MustBeThere) const
{
  TiXmlElement* pnProp = pnNode->FirstChildElement("property");
  while(pnProp)
  {
    size_t pos=0;
    const char* pAtt = pnProp->Attribute("dictRef");
    if(!pAtt)
    {
      pAtt = pnProp->Attribute("title");
      //The position of the start of the localname(the bit after the colon)
      pos=name.find(':') + 1; //ok even if there is no colon (if npos=-1)
    }
    if(pAtt && name.compare(pos, name.size()-pos, pAtt)==0)
    {
      TiXmlElement* pnChild = pnProp->FirstChildElement(); //could be <array> or <scalar> or <string>
      if(pnChild){
        return pnChild->Attribute(attName.c_str());
      }
    }
    pnProp = pnProp->NextSiblingElement();
  }
//  if(MustBeThere)
//    meErrorLog.ThrowError(__FUNCTION__, "The property " + name + " is missing or empty", obError);
  return NULL;
}

  /// Returns true if datatext associated with name is "1" or "true" or "yes" or nothing;
  //  returns false if datatext is something else or if element is not found.
  bool XMLPersist::XmlReadBoolean( const std::string& name)const
  {
    const char* txt = XmlReadValue(name, false);
    if(txt)
    {
      string s(txt);
      return s.empty() || s=="1" || s=="yes" || s=="true";
    }
    return false;
  }


PersistPtr XMLPersist::XmlWriteValueElement(const std::string& name,
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
  return PersistPtr(new XMLPersist(item));
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

//e.g. <metadata name="dc:source" content="LibraryMols.xml" timestamp="20080705_104810" />
PersistPtr XMLPersist::XmlWriteMetadata(const std::string& name, const std::string& content)
{
  TiXmlElement* item = new TiXmlElement("metadata");
  pnNode->LinkEndChild(item);
  TimeCount events;
  std::string timeString;
  item->SetAttribute("name",name);
  item->SetAttribute("name",name);
  item->SetAttribute("content",content);
  item->SetAttribute("timestamp",events.setTimeStamp(timeString));
  return PersistPtr(new XMLPersist(item));
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
    pnel->SetAttribute("calculated", events.setTimeStamp(timeString));
    //----------------------------------------

    //Add explanatory comment in description element
    TiXmlElement* pncom = new TiXmlElement("me:description");
    pnel->LinkEndChild( pncom );
    pncom->LinkEndChild(new TiXmlText(comment));
  }
  return PersistPtr(new XMLPersist(pnel));
}

///Insert into XML document a new property element
/**If the paramaeter units is not empty a timestamp and a units attribute are added. Like:
  
  <property dictRef="me:ZPE">
    <scalar calculated="20081122_230151" units="kJ/mol">139.5</scalar>
  </property>

  Returns a PersistPtr to the <scalar> element (so that more attributes can abe added).
**/
PersistPtr XMLPersist::XmlWriteProperty( const std::string& name, 
                                        const std::string& content, const std::string& units)
{
  TiXmlElement* pnprop = new TiXmlElement( "property" );
  pnNode->LinkEndChild(pnprop);
  pnprop->SetAttribute("dictRef", name);
  TiXmlElement* pnscal = new TiXmlElement( "scalar" );
  pnprop->LinkEndChild(pnscal);
  pnscal->LinkEndChild(new TiXmlText(content));
  if(!units.empty())
  {
    TimeCount events;
    std::string thisEvent, timeString;
    pnscal->SetAttribute("calculated", events.setTimeStamp(timeString));
    pnscal->SetAttribute("units", units);
  }
  return PersistPtr(new XMLPersist(pnscal));
}

bool XMLPersist::XmlCopyElement(PersistPtr ppToBeCopied)
{
  TiXmlNode* pnBefore = pnNode->FirstChild();
  if(!pnBefore) return false;
  XMLPersist* pxP = dynamic_cast<XMLPersist*> (ppToBeCopied.get());
  return pnNode->InsertBeforeChild(pnBefore, *pxP->pnNode);
}

bool XMLPersist::XmlSaveFile(const std::string& outfilename)
{
  if(outfilename.empty())
    return pnNode->GetDocument()->SaveFile(stdout);
  else
    return pnNode->GetDocument()->SaveFile(outfilename);
}

}//namespacer mesmer
