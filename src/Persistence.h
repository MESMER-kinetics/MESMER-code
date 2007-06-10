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
#include "formatfloat.h"


namespace mesmer
{

  // Global Utility methods:

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

  //Reading methods
  /// Returns an PersistPtr which can be used to read further down the input or output data
  virtual PersistPtr MoveTo(const std::string& name)const =0;

  ///Returns next item from the input document. Not used with XML IO.
  virtual const char* Read(){ return NULL; }

  /// Returns the value of the datatext associated with name
  virtual const char* ReadValue(const std::string& name, bool MustBeThere=true) const=0;

	/// Returns the data associated with name
 virtual const char* ReadProperty( const std::string& name, bool MustBeThere=true)const=0;

  //Writing methods
  /// Inserts into XML document a new element
  virtual PersistPtr WriteElement(const std::string& name)=0;

    /// Inserts into XML document a new element  containing a formatted number
  virtual void WriteValueElement(const std::string& name, 
                                 const double datum, const int precision=-1)=0;

  ///Insert into XML document a new element, name, and gives it a timestamp attribute and comment 
  virtual PersistPtr WriteMainElement( const std::string& name,
                                  const std::string& comment, bool replaceExisting=true)=0;

  virtual bool SaveFile(const std::string& outfilename)=0;

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
  PersistPtr(const PersistPtr& p) : p_(p.p_) { ++p_->count_; }
  PersistPtr() :p_(NULL){}
  PersistPtr& operator= (const PersistPtr& p)
  { // DO NOT CHANGE THE ORDER OF THESE STATEMENTS!
    // (This order properly handles self-assignment)
    // (This order also properly handles recursion, e.g., if a IPersist contains PersistPtrs)
    IPersist* const old = p_;
    p_ = p.p_;
    ++p_->count_;
    if (old && --old->count_ == 0) delete old;
    return *this;
  }

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
These additional PersistPtrs are mainly made by the MoveTo() function,
used like:
\code
   PersistPtr newPtr = oldPtr->MoveTo(new_element_name);
\endcode

You can copy and assign to PersistPtrs in the normal way. They are also
testable as a bool and 'NULL' pointers are returned from some I/O functions
generally to indicate non-existence of what was requested.

To allow possible use of other I/O formats and to simplify the interface,
the I/O functions are rather generalised and some of the details of the XML
structure are ignored. So, MoveTo(new_element_name) will find both child and
sibling elements called new_element_name; ReadValue() will return the content
of either a child element or an attribute.

A 'top-level' XMLPersist object (wrapped as always by  a PersistPtr) can be 
made by the static XMLPersist function Create(). It contains a pointer to the
TinyXml Document, and the lifetimes of the two objects are synchronised as is
required the document being deleted in the XMLPersist destructor. Most
'subsidiary' XMLPersist objects do not do this and are made in normal XMLPersist
functions. The constructor of XMLPersist is private, so that external instances,
which may subvert the reference-counting mechanism, cannot be made.

**/
#endif // GUARD_Persistence_h

