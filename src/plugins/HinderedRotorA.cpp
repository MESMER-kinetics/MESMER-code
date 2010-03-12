#include "../Molecule.h"
#include "HinderedRotorA.h"

using namespace std;
namespace mesmer
{
//************************************************************
//Global instance, defining its id
HinderedRotorA theHinderedRotorA("HinderedRotorA");
//************************************************************

//Read data from XML and store in this instance.
bool HinderedRotorA::ReadParameters(Molecule* pMol, PersistPtr ppDOSC)
{
  //if bond is needed
  const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);    
  if(bondID)
  {
    PersistPtr ppMol = pMol->get_PersistentPointer();
    PersistPtr ppBond = ppMol->XmlMoveTo("bondArray");
    const char* id = NULL;
    while(ppBond=ppBond->XmlMoveTo("bond"))
    {
      id = ppBond->XmlReadValue("id");
      if(id && !strcmp(id, bondID))
        break;
    }
    if(id)
    {
      m_bondID = id;
      cinfo << "Hindered rotor " << m_bondID << endl;
      //...etc.
    }
  }

  m_barrier  = ppDOSC->XmlReadDouble("me:barrierZPE",optional);
  PersistPtr pp = ppDOSC->XmlMoveTo("me:barrierZPE");
  const char* p = pp->XmlReadValue("units", optional);
  string units = p ? p : "kJ/mol";
  m_barrier = getConvertedEnergy(units, m_barrier);

  m_periodicity = ppDOSC->XmlReadInteger("me:periodicity",optional);

  m_vibFreq = ppDOSC->XmlReadDouble("me:vibFreq",optional);
  pp = ppDOSC->XmlMoveTo("me:vibFreq");
  p = pp->XmlReadValue("units", optional);
  units = p ? p : "cm-1";
//...etc
  return true;
}

//Adjusts DOS for the specified vibrations being treated as hindered rotors
bool HinderedRotorA::countCellDOS(gDensityOfStates* pDOS, int MaximumCell)
{
//Possibly use http://authors.library.caltech.edu/1154/1/MCCjcp97a.pdf ?

  vector<double> cellDOS;
  if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
    return false;

  //***TODO remove contribution from vibration and add contribution from hindered rotor

  pDOS->setCellDensityOfStates(cellDOS) ;

  return true;
}
/*
double HinderedRotorA::CalcReducedMomentOfInertia()
{
}
*/
}//namespace
