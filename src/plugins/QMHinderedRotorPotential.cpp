#include "Molecule.h"
#include "QMHinderedRotorPotential.h"

using namespace std;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  QMHinderedRotorPotential theQMHinderedRotorPotential("QMHinderedRotorPotential");
  //************************************************************

  //Adjusts DOS for the specified vibrations being treated as hindered rotors
  bool QMHinderedRotorPotential::countCellDOS(gDensityOfStates* pDOS, int MaximumCell)
  {
    //bond is needed, and the PES is recorded in the bond as an array.
    //const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);    
    const char* bondID = "";    // temporary
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
      // The PES information about the dihedral torsion of the bond should be recorded in here.
      if(ppBond){
        //ppBond->
      }
      else{
        cerr << "Cannot find the specified bond.";
      }
    }
    cinfo << "Hindered rotor " << bondID << endl;

    //double barrier  = ppDOSC->XmlReadDouble("me:barrierZPE");
    double barrier  = 0.0; // temporary
    //PersistPtr pp = ppDOSC->XmlMoveTo("me:barrierZPE");
    PersistPtr pp = NULL; // temporary
    //const char* p = pp->XmlReadValue("units", optional);
    const char* p = ""; //temporary
    string units = p ? p : "kJ/mol";
    barrier = getConvertedEnergy(units, barrier);

    //...
    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    //***TODO remove contribution from vibration and add contribution from hindered rotor

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;
  }

}//namespace
