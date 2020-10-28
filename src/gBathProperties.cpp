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
      ppPropList = pp; //Be forgiving; we can get by without a propertyList element

    m_Sigma   = ppPropList->XmlReadPropertyDouble("me:sigma", false);
    m_Epsilon = ppPropList->XmlReadPropertyDouble("me:epsilon", false);

    if (IsNan(m_Sigma) || IsNan(m_Epsilon)) {
      m_dafaultPrmtrs = true;
      m_Sigma = ppPropList->XmlReadPropertyDouble("me:sigma");
      m_Epsilon = ppPropList->XmlReadPropertyDouble("me:epsilon");
    }
  }

}
