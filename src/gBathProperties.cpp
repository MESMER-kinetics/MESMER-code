#include "Molecule.h"
#include "System.h"
#include "ParseForPlugin.h"
#include "gBathProperties.h"

using namespace std;
using namespace Constants;
using namespace OpenBabel;
namespace mesmer
{

  //-------------------------------------------------------------------------------------------------
  // Bath gas related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //

  gBathProperties::gBathProperties(Molecule* pMol)
    :m_Sigma(sigmaDefault),
    m_Epsilon(epsilonDefault),
    m_dafaultPrmtrs(false)
   {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // Be forgiving; we can get by without a propertyList element.

    m_Sigma   = ppPropList->XmlReadPropertyDouble("me:sigma", false);
    m_Epsilon = ppPropList->XmlReadPropertyDouble("me:epsilon", false);

    if (IsNan(m_Sigma) || IsNan(m_Epsilon)) {
      clog << "No lennard-Jones parameters found for " << pMol->getName() << "." << endl 
           << "An attempt to calculate parameters will be made. If this fails" << endl 
           << "default parameters will be used" << endl;
      m_dafaultPrmtrs = true;
      m_Sigma = ppPropList->XmlReadPropertyDouble("me:sigma");
      m_Epsilon = ppPropList->XmlReadPropertyDouble("me:epsilon");
    }
  }

}
