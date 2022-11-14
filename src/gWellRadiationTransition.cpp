//-------------------------------------------------------------------------------------------
//
// gWellRadiationTransition.h
//
// Author: Struan Robertson
// Date:   27/Dec/2021
//
// This header file contains the implementation of the gWellRadiationTransition class. This 
// class governs radition transitions betwen states.
//
//-------------------------------------------------------------------------------------------
#include "Molecule.h"
#include "gWellRadiationTransition.h"
#include "System.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  //
  // Constructor, destructor and initialization
  //
  gWellRadiationTransition::gWellRadiationTransition(Molecule* pMol) : MolecularComponent(), 
    m_TransitionFrequency(),
    m_EinsteinAij(),
    m_EinsteinBij(),
    m_bActivation(false),
    m_lowestBarrier(9e23)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
  }

  gWellRadiationTransition::~gWellRadiationTransition()
  {
  }

  bool gWellRadiationTransition::initialization() {

    PersistPtr pp = m_host->get_PersistentPointer();

    // Read in Einstein Aij and/or Bij coefficents.

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    const char* txt;
    if ((txt = ppPropList->XmlReadProperty("me:EinsteinAij", optional))) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_EinsteinAij.push_back(x);
    }

    if ((txt = ppPropList->XmlReadProperty("me:EinsteinBij", optional))) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_EinsteinBij.push_back(x);
    }

    if (m_EinsteinAij.size() == m_EinsteinBij.size()) {
      m_bActivation = false;
    } else if (m_EinsteinAij.size() == 0 && m_EinsteinBij.size() > 0) {
      m_bActivation = true;
    } else if (m_EinsteinAij.size() > 0 && m_EinsteinBij.size() == 0) {
      m_bActivation = false;
    } else {
      throw(runtime_error("Error reading Einstein Coefficients."));
    }

    return true;
  }

  // Set transition frequencies.
  void gWellRadiationTransition::setTransitionFrequencies(vector<double>& frequencies) {
    m_TransitionFrequency = frequencies;

    // Derive Bij Coefficients from Aij coefficients if requried.
    if (!m_bActivation && m_EinsteinBij.size() != m_EinsteinAij.size()) {
      m_EinsteinBij.clear();
      const double h = PlancksConstant_in_JouleSecond;
      const double c = SpeedOfLight_in_cm;
      const double A2B = c * c / (8.0 * M_PI * h);
      for (size_t i(0); i < m_EinsteinAij.size(); i++) {
        const double nu = m_TransitionFrequency[i] * c;
        const double Bij = A2B * m_EinsteinAij[i] / (pow(nu, 3.0));
        m_EinsteinBij.push_back(Bij);
      }
    }

    return;
  }



}//namespace
