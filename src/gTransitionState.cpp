#include <cmath>
#include <iomanip>
#include <numeric>
#include "Molecule.h"
#include "System.h"
#include "ParseForPlugin.h"
#include "gTransitionState.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //-------------------------------------------------------------------------------------------------
  // Transition state related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gTransitionState::~gTransitionState()
  {
    //if (m_ImFreq_chk == 0) cinfo << "m_ImFreq is provided but not used in " << m_host->getName() << "." << endl;
  }

  bool gTransitionState::initialization()
  {
    ErrorContext c(m_host->getName());
    PersistPtr pp = m_host->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;
    txt = ppPropList->XmlReadProperty("me:imFreqs", optional);
    if (!txt) {
      cinfo << "No imaginary vibrational frequency. Check if there is tunneling.\n";
      m_ImFreq_chk = -1;
      // Note: If there is tunneling, an imaginary frequency must also be supplied.
    }
    else {
      istringstream idata(txt); double x;
      // Make sure this number is positive.
      idata >> x;
      if (x > 0.0) {
        m_ImFreq = x;
        bool rangeSet;
        PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:imFreqs");
        ReadRdoubleRange(string(m_host->getName() + ":imFreqs"), ppProp, m_ImFreq, rangeSet);
        m_ImFreq_chk = 0;
      } else {
        string errorMsg = "Negative imaginary frequency detected for transition state " + m_host->getName() + ".";
        throw(std::runtime_error(errorMsg));
      }
    }
    return true;
  }

  double gTransitionState::get_ImFreq() {
    if (m_ImFreq_chk == -1) {
      cerr << "m_ImFreq was not defined but requested." << ". Default value " << m_ImFreq << " is given.";
      --m_ImFreq_chk;
      return m_ImFreq;
    }
    else if (m_ImFreq_chk < -1) {
      --m_ImFreq_chk;
      return m_ImFreq;
    }
    ++m_ImFreq_chk;
    return m_ImFreq;
  }

}//namespace
