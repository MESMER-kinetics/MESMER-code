//
// Molecule.cpp
//
// Author: Struan Robertson
// Date:   5/Jan/2003
//-------------------------------------------------------------------------------------------
#include "Molecule.h"
#include "gWellProperties.h"
#include "gWellRadiationTransition.h"
#include "gDensityOfStates.h"
#include "gStructure.h"
#include "gBathProperties.h"
#include "gPopulation.h"
#include "gTransitionState.h"

using namespace std;
using namespace Constants;

namespace mesmer
{

  //
  //Constructor
  //
  Molecule::Molecule(const MesmerEnv& Env, MesmerFlags& Flags, const string& molType, MoleculeManager* pMoleculeManager) : m_Env(Env),
    m_Flags(Flags),
    m_atomNumber(0),
    m_ppPersist(NULL),
    m_Name(),
    m_Description(),
    m_pMoleculeManager(pMoleculeManager),
    m_molTypes(),
    g_bath(NULL),
    g_dos(NULL),
    g_ts(NULL),
    g_pop(NULL),
    g_coll(NULL),
    g_rad(NULL),
    g_struc(NULL)
  {
    m_molTypes[molType] = true;
  }

  Molecule::~Molecule() {
    // delete the pointers in the reverse order.
    if (g_struc != NULL) delete g_struc;
    if (g_coll != NULL) delete g_coll;
    if (g_rad != NULL) delete g_rad;
    if (g_pop != NULL) delete g_pop;
    if (g_ts != NULL) delete g_ts;
    if (g_dos != NULL) delete g_dos;
    if (g_bath != NULL) delete g_bath;
  }

  //
  //Initialization
  //
  bool Molecule::InitializeMolecule(PersistPtr pp)
  {
    m_ppPersist = pp;
    const char* id = m_ppPersist->XmlReadValue("id", false);
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
      if (!m_ppPersist->XmlReadValue("", false))
        //Has neither description attribute nor any child elements
        return false;

    // check the number of atoms
    PersistPtr ppAtomArray = m_ppPersist->XmlMoveTo("atomArray");
    if (ppAtomArray) {
      PersistPtr ppAtom = ppAtomArray->XmlMoveTo("atom");
      while (ppAtom) {
        m_atomNumber++;
        ppAtom = ppAtom->XmlMoveTo("atom");
      }
    }

    return true;
  }

  // Check whether the total degrees of freeedom are consistent with the number of atoms.
  bool Molecule::checkDegOfFreedom() {

    if (m_atomNumber < 1) // No atoms defined, so no test can be done.
      return true;

    unsigned int nDOF(0);
    if (g_dos) {
      nDOF = g_dos->getNoOfDegOfFreeedom();
    }
    else {
      // Maybe a bath gas molecule.
      return true;
    }
    if (m_atomNumber > 1 && nDOF == 0) {
      cinfo << "No molecular data assigned for " << getName() << ", assuming it to be a sink term.\n";
      return true;
    }

    // If species is a transition state, assume a imaginary frequency exists regardless if provided.
    if (g_ts)
      nDOF++;

    // Add translational degrees of freedom.

    nDOF += 3;

    return (nDOF == (3 * m_atomNumber));
  }

  //Make dummy calls to initialize the MolecularComponents now, during parse, before calculation.
  bool Molecule::activateRole(string molType) {
    // see if the molType is true in m_molTypes
    if (!m_molTypes[molType])
    {
      m_molTypes[molType] = true;
    }

    if (molType == "bathGas" || molType == "modelled")
      getBath();

    if (molType == "transitionState")
      getTS();

    if (molType == "modelled" || molType == "deficientReactant" //|| molType == "sink"   CM g_dos to be added later if needed for reverse ILT
      || molType == "transitionState" || molType == "excessReactant"
      || molType == "PriorCoFragment" || molType == "forThermo")
      getDOS();

    if (molType == "modelled" || molType == "deficientReactant" || molType == "sink")
      getPop();

    if (molType == "modelled")
      getColl();

    if (molType == "modelled" || molType == "bathGas")
      getStruc();

    return true;
  };

  gBathProperties& Molecule::getBath() {
    if (!g_bath)
      g_bath = new gBathProperties(this);
    return *g_bath;
  }

  gDensityOfStates& Molecule::getDOS() {
    if (!g_dos) {
      g_dos = new gDensityOfStates(this);
      if (!g_dos->initialization()) {
        cinfo << "gDensityOfStates initialization failed." << endl;
      }
    }
    return *g_dos;
  }

  gTransitionState& Molecule::getTS() {
    if (!g_ts) {
      g_ts = new gTransitionState(this);
      if (!g_ts->initialization()) {
        cinfo << "gTransitionState initialization failed." << endl;
      }
    }
    return *g_ts;
  }

  gPopulation& Molecule::getPop() {
    if (!g_pop)
      g_pop = new gPopulation(this);
    return *g_pop;
  }

  gWellProperties& Molecule::getColl() {
    if (!g_coll) {
      g_coll = new gWellProperties(this);

      if (!g_coll->initialization()) {
        string err("Collision model parameters not defined for " + m_Name + ".\n");
        throw std::runtime_error(err);
      }
    }

    return *g_coll;
  }

  gWellRadiationTransition& Molecule::getRad() {
    if (!g_rad) {
      g_rad = new gWellRadiationTransition(this);

      if (!g_rad->initialization()) {
        string err("Radiation model parameters not defined for " + m_Name + ".\n");
        throw std::runtime_error(err);
      }
    }

    return *g_rad;
  }

  gStructure& Molecule::getStruc() {
    if (!g_struc)
      g_struc = new gStructure(this);

    return *g_struc;
  }

}//namespace
