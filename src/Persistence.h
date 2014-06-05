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
#include <sstream>
#include <time.h>


namespace mesmer
{
  // Although the methods have names starting with Xml, the Ipersist interface 
  // could potentially be used with non-XML files.
  // The derived class XMLPersist is XML specific.

  class PersistPtr;

  class IPersist
  {
  public:
    friend class PersistPtr;
    unsigned count_;

  public:
    IPersist() : count_(0){}
    virtual ~IPersist(){} ;
    virtual operator bool() const=0;

    ///Incorporates the content of the file with name filename into the main data file.
    ///Returns false if there is an error.
    virtual bool XmlInclude(const std::string& filename)=0;

    //Reading methods
    /// Returns a PersistPtr which can be used to read further down the input or output data
    //  (child or sibling element)
    virtual PersistPtr XmlMoveTo(const std::string& name)const =0;

    ///Returns next item from the input document. (Value of this element) or NULL
    virtual const char* XmlRead()const=0;

    /// Returns the value of the datatext associated with name (Value of child element or attribute) or NULL
    virtual const char* XmlReadValue(const std::string& name, bool MustBeThere=true, const char* pluginTypeID=NULL)=0;

    /// Returns the value of a specified child element or attribute.
    // If the element or attribute is not found and the default is either not requested
    // or is not found, then NaN is returned.
    virtual double XmlReadDouble(const std::string& name, bool MustBeThere=true)=0;

    /// Returns the value of a specified child element or attribute.
    virtual int XmlReadInteger(const std::string& name, bool MustBeThere=true)=0;

    /// Returns the data associated with name (CML property element) or NULL
    virtual const char* XmlReadProperty( const std::string& name, bool MustBeThere=true)=0;

    /// Returns a PersistPtr to the scalar or array  element of a CML property element.
    //  To call from withhin a previous property set NextProperty true.
    virtual PersistPtr XmlMoveToProperty(const std::string& name,bool ToNextProperty=0)=0;

    /// Returns the value associated with name (CML property element)
    // If the property is not found and the default is either not requested
    // or is not found, then NaN is returned.
    virtual double XmlReadPropertyDouble(const std::string& name, bool MustBeThere=true)=0;

    /// Returns the value associated with name (CML property element)
    virtual int  XmlReadPropertyInteger(const std::string& name, bool MustBeThere=true)=0;

    /// Returns the attribute associated with name and attName (CML property element) or NULL
    virtual const char* XmlReadPropertyAttribute(const std::string& name, const std::string& attName, bool MustBeThere=true)=0;

    /// Returns true if datatext associated with name is "1" or "true" or nothing;
    //  returns false if datatext is something else or if element is not found.
    virtual bool XmlReadBoolean( const std::string& name)=0;

    //Writing methods

    /// Writes a value to the element
    virtual void XmlWrite(const std::string& value)=0;

    /// Inserts into XML document a new element
    virtual PersistPtr XmlWriteElement(const std::string& name)=0;

    /// Adds an XML attribute (or equivalent)
    virtual void XmlWriteAttribute(const std::string& name, const std::string& value)=0;

    /// Adds an XML attribute containing a formatted number
    virtual void XmlWriteAttribute(const std::string& name,
                     const double datum, const int precision, const bool fixedOnly)=0;

    /// Inserts into XML document a new element containing a formatted number
    virtual PersistPtr XmlWriteValueElement(const std::string& name,
                     const double datum, const int precision=-1, const bool fixedOnly=false)=0;

    /// Inserts into XML document a new element  containing a string
    virtual PersistPtr XmlWriteValueElement(const std::string& name, const std::string& value,
                                            bool cdatamode=false)=0;

    /// Inserts into XML document meta data information
    ///like <metadata name="dc:source" content="LibraryMols.xml" timestamp="20080705_104810" />
    virtual PersistPtr XmlWriteMetadata(const std::string& name, const std::string& content)=0;

    ///Insert into XML document a new element, name, and gives it a timestamp attribute and comment (if comment not empty)
    virtual PersistPtr XmlWriteMainElement( const std::string& name,
                                   const std::string& comment, bool replaceExisting=true)=0;

    ///Insert into XML document a new property element, with timestamp and units attributes (if units is not empty)
    virtual PersistPtr XmlWriteProperty( const std::string& name, 
                                  const std::string& content, const std::string& units)=0;

    ///Replace or insert an element and return a pointer to the copy
    virtual PersistPtr XmlCopy(PersistPtr ppToBeCopied, PersistPtr ppToBeReplaced)=0;

    virtual bool XmlSaveFile(const std::string& outfilename)=0;

    ///Utility function to return a string with the current time and date.
    //static std::string TimeString()
    //{
    //  time_t ltime;
    //  time( &ltime );
    //  std::string timestring(ctime(&ltime));
    //  return timestring.erase(timestring.size()-1);
    //}

  } ;

  //Reference counting taken from:
  //http://www.parashift.com/c++-faq-lite/freestore-mgmt.html#faq-16.22
  //modified to handle p_==NULL, and with bool conversion
  class PersistPtr
  {
  public:
    operator bool()const
    {
      return p_ && *p_;
    }
    IPersist* operator-> () { return p_; }
    IPersist& operator* ()  { return *p_; }
    PersistPtr(IPersist* p) : p_(p)
    {
      if(p) ++p_->count_; //p can be NULL
    }
    ~PersistPtr()           { if (p_ && --p_->count_ == 0) delete p_; }
    PersistPtr(const PersistPtr& p) 
    {
      if(p)
      {
        p_ = p.p_;
        ++p_->count_;
      }
      else
        p_ = NULL;
    }
    PersistPtr() :p_(NULL){}
    PersistPtr& operator= (const PersistPtr& p)
    { // DO NOT CHANGE THE ORDER OF THESE STATEMENTS!
      // (This order properly handles self-assignment)
      // (This order also properly handles recursion, e.g., if a IPersist contains PersistPtrs)
      if(p)
      {
        IPersist* const old = p_;
        p_ = p.p_;
        ++p_->count_;
        if (old && --old->count_ == 0) delete old;
      }
      else
        p_ = NULL;
      return *this;
    }

    IPersist* get() {return p_;};//needed for XmlCopyElement()

  private:
    IPersist* p_;
  };

}

/**
The input/output code is designed to minimize dependencies between files,
and to make replacement of the XML library easier should this be necessary.

All I/O is done by calling virtual functions of the abstract class IPersist.
A derived class XMLPersist implements this interface for XML files using the
TinyXML library, but alternative implementations or even different I/O formats
could be used by replacing the code for XMLPersist in xmlpersist.h and
xmlpersist.cpp.

Alternatively, another class derived from IPersist could be written and the
line in main() that references XMLPersist modified. This is the only reference
to XMLPersist, xmlpersist.h is #included only in main.cpp. Also, tinyxml.h is
#included only in xmlpersist.h.

An XMLPersist object can be regarded as encapsulating a pointer to the XML
tree. Because it is convenient to have simultaneously pointers to several
different parts of the tree, there need to be several XMLPersist objects.
Their scopes need to be wider than a single function, so that they need to
be made on the heap with new. Deleting them when they are finished with, to
avoid memory leaks, is a challenge, but can be met by using a reference-counted
smart pointer. The class shared_ptr will be in the next C++ standard but in
the meantime a library like Boost needs to be used to provide a cross-platform
implementation. This is a bit heavyweight, so that a built-in reference-counting
mechanism has been used, based on an example in C++ FAQ Lite. All mentions of
XMLPersist (or any Ipersist derived class) are wrapped by this smart pointer,
PersistPtr. So this becomes the only I/O class directly referred to in the
main part of Mesmer.

A PersistPtr object referencing the root of the XML document is generated
in main() and passed as a parameter in System::parse(). Additional instances
of PersistPtr which point to other locations in the XML tree are generated
during the parsing. For example, one that points to a <molecule> element
is stored as a member variable in an instance of the Molecule class.
These additional PersistPtrs are mainly made by the XmlMoveTo() function,
used like:
\code
   PersistPtr newPtr = oldPtr->XmlMoveTo(new_element_name);
\endcode

You can copy and assign to PersistPtrs in the normal way. They are also
testable as a bool and 'NULL' pointers are returned from some I/O functions
generally to indicate non-existence of what was requested.

To allow possible use of other I/O formats and to simplify the interface,
the I/O functions are rather generalised and some of the details of the XML
structure are ignored. So, XmlMoveTo(new_element_name) will find both child and
sibling elements called new_element_name; XmlReadValue() will return the content
of either a child element or an attribute.

A 'top-level' XMLPersist object (wrapped as always by  a PersistPtr) can be
made by the static XMLPersist function XmlLoad(). It contains a pointer to the
TinyXml Document, and the lifetimes of the two objects are synchronised as is
required, the document being deleted in the XMLPersist destructor. Most
'subsidiary' XMLPersist objects do not do this and are made in normal XMLPersist
functions. The constructor of XMLPersist is private, so that external instances,
which may subvert the reference-counting mechanism, cannot be made.

**/
#endif // GUARD_Persistence_h

