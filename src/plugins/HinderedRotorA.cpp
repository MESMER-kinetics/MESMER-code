#include "../Molecule.h"
#include "HinderedRotorA.h"

using namespace std;
namespace mesmer
{
//-------------------------------------------------------------
//Global instance, defining its id
HinderedRotorA theHinderedRotorA("HinderedRotorA");
//-------------------------------------------------------------
using OpenBabel::vector3;
//Read data from XML and store in this instance.
bool HinderedRotorA::ReadParameters(Molecule* pMol, PersistPtr ppDOSC)
{
  gStructure& gs = pMol->getStruc();
  if(!gs.ReadStructure())
  {
    cerr << "A complete set of atom coordinates are required for hindered rotor calculations" <<endl;
    return false;
  }

  const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);
  if(bondID)
  {
    pair<string,string> bondats = gs.GetAtomsOfBond(bondID);
    if(bondats.first.empty())
    {
      cerr << "Unknown bond reference " << bondID << endl;
      return false;
    }
    m_bondID = bondID;
    cinfo << "Hindered rotor " << m_bondID << endl;

    vector3 coords1 = gs.GetAtomCoords(bondats.first);
    vector3 coords2 = gs.GetAtomCoords(bondats.second);

    //calc Moment of inertia about bond axis of atoms on one side of bond...
    vector<string> atomset;
    gs.GetAttachedAtoms(atomset, bondats.first, bondats.second);
    double mm1 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

    //...and the other side of the bond
    atomset.clear();
    gs.GetAttachedAtoms(atomset, bondats.second, bondats.first);
    double mm2 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

    /*
    Is the reduced moment of inertia need about the bond axis or, separately for the set of
    atoms on each side of the bond, about a parallel axis through their centre of mass?
    See:
    http://www.ccl.net/chemistry/resources/messages/2001/03/21.005-dir/index.html
    http://www.ccl.net/chemistry/resources/messages/2001/03/31.002-dir/index.html
    The bond axis is used here.
    */
    m_reducedMomentInertia = mm1 * mm2 / ( mm1 + mm2 );//units a.u.*Angstrom*Angstrom
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

  // other inputs...
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

}//namespace
