//
// Molecule.cpp
//
// Author: Struan Robertson
// Date:   5/Jan/2003
//-------------------------------------------------------------------------------------------
#include "Molecule.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{

  //
  //Constructor
  //
  Molecule::Molecule(const MesmerEnv& Env, MesmerFlags& Flags):
    m_Env(Env),
    m_Flags(Flags),
    m_errorflag(0),
    m_ppPersist(NULL),
    m_Name(),
    m_Description()
  {}

  Molecule::~Molecule(){}
  /* Will need Clone() functions
  Molecule::Molecule(const Molecule& molecule) {
  // Copy constructor - define later SHR 23/Feb/2003
  }

  Molecule& Molecule::operator=(const Molecule&) {
  // Assignment operator - define later SHR 23/Feb/2003

  return *this ;
  }
  */

  //
  //Initialization
  //
  bool Molecule::InitializeMolecule(PersistPtr pp)
  {
    m_ppPersist = pp;
    const char* id= m_ppPersist->XmlReadValue("id");
    if (id) m_Name = id;
    if (m_Name.empty()) {
      cerr << "Molecular name is absent.";
      m_Name = "unknown";
      setFlag(true);
    }

    const char* desc = m_ppPersist->XmlReadValue("description");
    if (desc)
      m_Description = desc;
    // no check value for description

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    //const char* txt;

    if (getErrorFlag()){
      cerr << "Error(s) while initializing molecule: " << m_Name;
      return false;
    }
    return true;
  }

  PersistPtr  Molecule::getPersistentPointer()     {
    return m_ppPersist;
  };
  void Molecule::setPersistentPointer(PersistPtr value){
    m_ppPersist = value;
  }
  std::string Molecule::getName() const            {
    return m_Name ;
  } ;
  std::string Molecule::getDescription() const            {
    return m_Description ;
  } ;
  MesmerFlags& Molecule::getFlags()         {
    return m_Flags;
  }

  const MesmerEnv& Molecule::getEnv() const        {
    return m_Env;
  }

  int    Molecule::getErrorFlag()                       {
    return m_errorflag;
  } ;

  void   Molecule::setName(string name)            {
    m_Name = name;
  } ;
  void   Molecule::setFlag(bool value)             {
    if (value) ++m_errorflag;
  } ;



}//namespace
