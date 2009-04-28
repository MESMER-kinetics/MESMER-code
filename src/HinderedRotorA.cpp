#include "Molecule.h"
#include "HinderedRotorA.h"

using namespace std;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  HinderedRotorA theHinderedRotorA("HinderedRotorA");
  //************************************************************

  //Adjusts DOS for the specified vibrations being treated as hindered rotors
  bool HinderedRotorA::countCellDOS(gDensityOfStates* pDOS, int MaximumCell, PersistPtr ppDOSC)
  {
    //if bond is needed
    const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);    
    if(bondID)
    {
      PersistPtr ppMol = pDOS->getHost()->get_PersistentPointer();
      PersistPtr ppBond = ppMol->XmlMoveTo("bondArray");
      while(ppBond=ppBond->XmlMoveTo("bond"))
      {
        const char* id = ppBond->XmlReadValue("id");
        if(id && !strcmp(id, bondID))
          break;
      }
      //if(ppBond)...
    }
    cinfo << "Hindered rotor " << bondID << endl;

    double barrier  = ppDOSC->XmlReadDouble("me:barrierZPE",optional);
    PersistPtr pp = ppDOSC->XmlMoveTo("me:barrierZPE");
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "kJ/mol";
    barrier = getConvertedEnergy(units, barrier);

    int periodicity = ppDOSC->XmlReadInteger("me:periodicity",optional);
    //...
    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    //***TODO remove contribution from vibration and add contribution from hindered rotor

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;
  }

}//namespace
