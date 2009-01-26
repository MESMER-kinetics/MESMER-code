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
    m_ppPersist(NULL),
    m_Name(),
    m_Description(),
    m_Mass(),
    m_Mass_chk(),
    g_bath(NULL),
    g_dos(NULL),
    g_ts(NULL),
    g_pop(NULL),
    g_coll(NULL)
  {}

  Molecule::~Molecule(){
    // delete the pointers in the reverse order.
    if (g_coll  != NULL) delete g_coll;
    if (g_pop   != NULL) delete g_pop;
    if (g_ts    != NULL) delete g_ts;
    if (g_dos  != NULL) delete g_dos;
    if (g_bath  != NULL) delete g_bath;
    if (m_Mass_chk == 0){
      cinfo << "m_Mass is provided but not used in " << getName() << endl;
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
      return false;
    }

    const char* desc = m_ppPersist->XmlReadValue("description");
    if (desc)
      m_Description = desc;
    // no check value for description

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt = ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      cerr << "Cannot find argument me:MW in " << getName() << ".";
      return false; // later put a function to calculate the molecular weight if the user forgot to provide it.
    }
    else { istringstream idata(txt); double mass(0.); idata >> mass; setMass(mass);}

    return true;
  }

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

  PersistPtr  Molecule::get_PersistentPointer()     {
    return m_ppPersist;
  };
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

  void   Molecule::setName(string name)            {
    m_Name = name;
  } ;

  bool Molecule::activateRole(string molType){
    if (molType == "bathGas" || molType == "modelled"){
      if (!g_bath){
        gBathProperties* pGbath = new gBathProperties();
        g_bath = pGbath;
        if (!g_bath->InitializeProperties(m_ppPersist, this))
          return false;
      }
    }
    if (molType == "modelled" || molType == "deficientReactant" || molType == "sink"
      || molType == "transitionState" || molType == "excessReactant"){
      if (!g_dos){
        gDensityOfStates* pGdos = new gDensityOfStates();
        g_dos = pGdos;
        if (!g_dos->InitializeProperties(m_ppPersist, this))
          return false;
      }
    }
    if (molType == "transitionState"){
      if (!g_ts  ){
        gTransitionState* pGts = new gTransitionState();
        g_ts = pGts;
        if (!g_ts->InitializeProperties(m_ppPersist, this))
          return false;
      }
    }
    if (molType == "modelled" || molType == "deficientReactant" || molType == "sink"){
      if (!g_pop ){
        gPopulation* pGpop = new gPopulation();
        g_pop = pGpop;
        if (!g_pop->InitializeProperties(m_ppPersist, this))
          return false;
      }
    }
    if (molType == "modelled"){
      if (!g_coll){
        gWellProperties* pGcoll = new gWellProperties();
        g_coll = pGcoll;
        if (!g_coll->InitializeProperties(m_ppPersist, this))
          return false;
      }
    }
    return true;
  };

  bool Molecule::roleIsActivated(string molType){
    if (molType == "modelled"){
      if (!g_bath || !g_dos || !g_ts || !g_pop || !g_coll)
        return false;
    }
    else if (molType == "deficientReactant" || "sink"){
      if (!g_dos || !g_pop)
        return false;
    }
    else if (molType == "transitionState"){
      if (!g_dos || !g_ts)
        return false;
    }
    else if (molType == "excessReactant"){
      if (!g_dos)
        return false;
    }
    else if (molType == "bathGas"){
      if (!g_bath)
        return false;
    }
    return true;
  };


  gBathProperties&        Molecule::getBath() {
    if (g_bath) return *g_bath;
    cerr << "Not a proper call to the molecular type defined in the reaction list.";
    exit(1);
  }
  gDensityOfStates&       Molecule::getDOS()  {
    if (g_dos) return *g_dos;
    cerr << "Not a proper call to the molecular type defined in the reaction list.";
    exit(1);
  }
  gTransitionState&       Molecule::getTS()   {
    if (g_ts) return *g_ts  ;
    cerr << "Not a proper call to the molecular type defined in the reaction list.";
    exit(1);
  }
  gPopulation&            Molecule::getPop()  {
    if (g_pop) return *g_pop ;
    cerr << "Not a proper call to the molecular type defined in the reaction list.";
    exit(1);
  }
  gWellProperties&        Molecule::getColl() {
    if (g_coll) return *g_coll;
    cerr << "Not a proper call to the molecular type defined in the reaction list.";
    exit(1);
  }

}//namespace
