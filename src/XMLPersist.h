#ifndef GUARD_XMLPersist_h
#define GUARD_XMLPersist_h

#include "Persistence.h"
#include "tinyxml.h"

namespace mesmer
{
class XMLPersist : public IPersist
{
protected:
  TiXmlElement* pnNode;
  TiXmlDocument* pDocument; //NULL except when XMLPersist made from Create()

private:
  //Constructor private so that reference counting cannot be subverted.
  //XMLPersist objects are made from outside only in Create().
  //They can also be made internally, e.g. in MoveTo().
  XMLPersist(TiXmlElement* pn=NULL, TiXmlDocument* pdoc=NULL) : pnNode(pn), pDocument(pdoc){}
public:
  ~XMLPersist();
  operator bool() const
  {
    return pnNode!=NULL;
  }

///Makes an instance of XMLPersist. Opens the file and reads the contents.
///If title is not empty gives an error message if the title or root element 
///does not have this name.
  static PersistPtr Create(const std::string& filename, const std::string& title="");


///Returns the first child element with this name, or if not found, the next sibling element with this name
  virtual PersistPtr MoveTo(const std::string& name) const;

  virtual const char* ReadValue(const std::string& name, bool MustBeThere=true) const;
  virtual const char* ReadProperty( const std::string& name, bool MustBeThere=true) const;

  /// Inserts into XML document a new element
  virtual PersistPtr WriteElement(const std::string& name);

  /// Inserts into XML document a new element  containing a formatted number
  virtual void WriteValueElement(const std::string& name, 
                                 const double datum, const int precision=-1);

  ///Insert into XML document a new element, name, and gives it a timestamp attribute and comment 
  virtual PersistPtr WriteMainElement( const std::string& name,
                                  const std::string& comment, bool replaceExisting=true);

  virtual bool SaveFile(const std::string& outfilename);

};

}//namespace mesmer
#endif //GUARD_XMLPersist_h
