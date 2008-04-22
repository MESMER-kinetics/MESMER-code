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
  Molecule::Molecule(const MesmerEnv& Env):
    m_Env(Env),
    m_flag(0),
    m_ppPersist(NULL),
    m_Name(),
    m_Description(),
    m_Mass(0.0),
    m_Sigma(sigmaDefault),
    m_Epsilon(epsilonDefault),
    m_Mass_chk(-1),
    m_Sigma_chk(-1),
    m_Epsilon_chk(-1)
  {}

  Molecule::~Molecule()
  {
    if (m_Mass_chk == 0){
      cinfo << "m_Mass is provided but not used in " << getName();
    }
    if (m_Sigma_chk == 0){
      cinfo << "m_Sigma is provided but not used in " << getName();
    }
    if (m_Epsilon_chk == 0){
      cinfo << "m_Epsilon is provided but not used in " << getName();
    }
  }

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

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      cerr << "Cannot find argument me:MW in " << getName();
      setFlag(true); // later put a function to calculate the molecular weight if the user forgot to provide it.
    }
    else { istringstream idata(txt); double mass(0.); idata >> mass; setMass(mass);}


    if (getFlag()){
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
  const MesmerEnv& Molecule::getEnv() const        {
    return m_Env;
  }

  int    Molecule::getFlag()                       {
    return m_flag;
  } ;

  void   Molecule::setMass(double value)           {
    m_Mass = value;
    m_Mass_chk = 0;
  } ;

  double Molecule::getMass()                       {
    if (m_Mass_chk >= 0){
      ++m_Mass_chk;
      return m_Mass ;
    }
    else{
      cerr << "m_Mass was not defined but requested in " << getName();
      exit(1);
    }
  } ;

  void   Molecule::setSigma(double value)          {
    m_Sigma = value;
    m_Sigma_chk = 0;
  } ;

  double Molecule::getSigma()                      {
    if (m_Sigma_chk >= 0){
      ++m_Sigma_chk;
      return m_Sigma ;
    }
    else{
      cerr << "m_Sigma was not defined but requested in " << getName()
               << ". Default value " << sigmaDefault << " is used.\n";
      //exit(1);
      return m_Sigma ;
    }
  } ;

  void   Molecule::setEpsilon(double value)        {
    m_Epsilon = value;
    m_Epsilon_chk = 0;
  } ;

  double Molecule::getEpsilon()                    {
    if (m_Epsilon_chk >= 0){
      ++m_Epsilon_chk;
      return m_Epsilon ;
    }
    else{
      cerr << "m_Epsilon was not defined but requested in " << getName()
               << ". Default value " << epsilonDefault << " is used.\n";
      //exit(1);
      return m_Epsilon ;
    }
  } ;

  void   Molecule::setName(string name)            {
    m_Name = name;
  } ;
  void   Molecule::setFlag(bool value)             {
    if (value) ++m_flag;
  } ;



}//namespace
