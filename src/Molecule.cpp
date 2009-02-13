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
    g_bath(NULL),
    g_dos(NULL),
    g_ts(NULL),
    g_pop(NULL),
    g_coll(NULL),
    g_struc(NULL)
  {}

  Molecule::~Molecule(){
    // delete the pointers in the reverse order.
    if (g_struc != NULL) delete g_struc;
    if (g_coll  != NULL) delete g_coll;
    if (g_pop   != NULL) delete g_pop;
    if (g_ts    != NULL) delete g_ts;
    if (g_dos  != NULL) delete g_dos;
    if (g_bath  != NULL) delete g_bath;
    }
  
  //
  //Initialization
  //
  bool Molecule::InitializeMolecule(PersistPtr pp)
  {
    m_ppPersist = pp;
    const char* id= m_ppPersist->XmlReadValue("id", false);
    if (id) m_Name = id;
    if (m_Name.empty()) {
      cerr << "Molecular name is absent.";
      return false;
    }

    const char* desc = m_ppPersist->XmlReadValue("description", optional);
    if (desc)
      m_Description = desc;
    // no check value for description
    else
      if(!m_ppPersist->XmlReadValue("", false))
        //Has neither description attribute nor any child elements
        return false;

    return true;
  }

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

  //Make dummy calls to initialize the MolecularComponents now, during parse, before calculation.
  bool Molecule::activateRole(string molType){
    if (molType == "bathGas" || molType == "modelled")
      getBath();

    if (molType == "modelled" || molType == "deficientReactant" //  || molType == "sink" CM g_dos to be added later if needed for reverse ILT
      || molType == "transitionState" || molType == "excessReactant")
      getDOS();
    
    if (molType == "transitionState")
      getTS();

    if (molType == "modelled" || molType == "deficientReactant" || molType == "sink")
      getPop();

    if (molType == "modelled")
      getColl();

    if(molType == "modelled" || molType == "bathGas")
      getStruc();

    return true;
  };

  gBathProperties&        Molecule::getBath() {
    if (!g_bath)
      g_bath = new gBathProperties(this);
    return *g_bath;
  }
  

  gDensityOfStates&       Molecule::getDOS()  {
    if (!g_dos)
      g_dos = new gDensityOfStates(this);
     return *g_dos;
  }
  gTransitionState&       Molecule::getTS()   {
    if (!g_ts) 
      g_ts = new gTransitionState(this);
   return *g_ts  ;  }

  gPopulation&            Molecule::getPop()  {
    if (!g_pop)
      g_pop = new gPopulation(this);
    return *g_pop ;
  }
  gWellProperties&        Molecule::getColl() {
    if (!g_coll)
      g_coll = new gWellProperties(this);
     return *g_coll;
  }
  gStructure&        Molecule::getStruc() {
    if (!g_struc)
      g_struc = new gStructure(this);
     return *g_struc;
  }

}//namespace
