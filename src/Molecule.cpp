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
  Molecule::Molecule(const MesmerEnv& Env, MesmerFlags& Flags, const string& molType):
    m_Env(Env),
    m_Flags(Flags),
    m_atomNumber(0),
    m_ppPersist(NULL),
    m_Name(),
    m_Description(),
    m_molTypes(),
    g_bath(NULL),
    g_dos(NULL),
    g_ts(NULL),
    g_pop(NULL),
    g_coll(NULL),
    g_struc(NULL)
  {
    m_molTypes[molType] = true;
  }

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

    // check the number of atoms
    PersistPtr ppAtomArray = m_ppPersist->XmlMoveTo("atomArray");
    if (ppAtomArray){
      PersistPtr ppAtom = ppAtomArray->XmlMoveTo("atom");
      while (ppAtom){
        m_atomNumber++;
        ppAtom = ppAtom->XmlMoveTo("atom");
      }
    }

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

  // check whether the total vibrational frequencies (include the imaginary freq) number = 3N-6
  bool Molecule::checkFrequencies(){
    std::vector<double> vibFreq;
    if (g_dos){
      g_dos->get_VibFreq(vibFreq);
    }
    size_t numFreq = vibFreq.size();
    if (m_atomNumber > 1 && numFreq == 0){
      cinfo << "No vibrational frequencies were assigned for " << getName() << ", assuming it to be a sink term.\n";
      return true;
    }

    // if it is a transition state, assuming a imaginary frequency exists no matter user provide it or not.
    if (g_ts) numFreq++;

    // get symmety number
    std::vector<double> mmtsInt;
    int rotorType(2), isLinearRotor(0);
    if (g_dos){
      rotorType = g_dos->get_rotConsts(mmtsInt);
      if (rotorType == 0)
        isLinearRotor = 1;
    }
    else{
      // Maybe a bath gas molecule.
      return true;
    }

    switch (m_atomNumber){
      case 0:
        break; // provide no check to cases where user did not provide a frequencies vector.
      case 1:
        if (numFreq != 0)
          return false; // single atom should not have vibrational frequencies
        break;
      case 2:
        if (numFreq != (3 * m_atomNumber - 5))
          return false; // must be a linear rotor
        break;
      default:
        if (numFreq != (3 * m_atomNumber - 6 + isLinearRotor))
          return false;
    }
    return true;
  }

  bool Molecule::isCemetery()
  {
    return g_coll && g_coll->isCemetery();
  }

  const MesmerEnv& Molecule::getEnv() const        {
    return m_Env;
  }

  //Make dummy calls to initialize the MolecularComponents now, during parse, before calculation.
  bool Molecule::activateRole(string molType){
    // see if the molType is true in m_molTypes
    if (!m_molTypes[molType])
    {
      m_molTypes[molType] = true;
    }

    if (molType == "bathGas" || molType == "modelled")
      getBath();

    if (molType == "modelled" || molType == "deficientReactant" || molType == "sink" // CM g_dos to be added later if needed for reverse ILT
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
    if (!g_coll) {
      g_coll = new gWellProperties(this);

      if (!g_coll->initialization()) {
        // Throw error
      }
    }

    return *g_coll;
  }
  
  gStructure&        Molecule::getStruc() {
    if (!g_struc)
      g_struc = new gStructure(this);

     return *g_struc;
  }

}//namespace
