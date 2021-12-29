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
using namespace OpenBabel;

namespace mesmer
{
  //
  // Constructor, destructor and initialization
  //
  gWellRadiationTransition::gWellRadiationTransition(Molecule* pMol) : MolecularComponent(), 
    m_EinsteinBij(),
    m_lowestBarrier(9e23),
    m_numGroupedGrains(0)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
  }

  gWellRadiationTransition::~gWellRadiationTransition()
  {
  }

  bool gWellRadiationTransition::initialization() {

    PersistPtr pp = m_host->get_PersistentPointer();

    // Read in Einstein Bij coefficents.

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    const char* txt;
    if ((txt = ppPropList->XmlReadProperty("me:EinsteinBij", optional))) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_EinsteinBij.push_back(x);
    }

    // Determine the energy of the reservoir grain, first find the lowest barrier 
    // associated with the current well.
    // 
    PersistPtr ppReservoirSize = pp->XmlMoveTo("me:reservoirSize");

    m_numGroupedGrains = 0; // Reset the number of grains grouped into a reservoir grain to zero.

    if (ppReservoirSize) {

      // Check the size of the reservoir.
      double tmpvalue = pp->XmlReadDouble("me:reservoirSize");

      const char* unitsTxt = ppReservoirSize->XmlReadValue("units", false);
      string unitsInput("kJ/mol");
      if (unitsTxt) {
        unitsInput = unitsTxt;
      }
      else {
        stest << "No unit for reservoir size has been supplied, use kJ/mol." << endl;
      }

      const double value(getConvertedEnergy(unitsInput, tmpvalue));
      int grainLoc(int(value / double(m_host->getEnv().GrainSize)));
      int lowestBarrier = int(getLowestBarrier() / double(m_host->getEnv().GrainSize));

      if (grainLoc > 0) {
        if (grainLoc > lowestBarrier) {
          stest << "The reservoir size provided is too high, corrected according to the lowest barrier height." << endl;
          grainLoc = lowestBarrier;
        }
      }
      else {
        if (abs(grainLoc) > lowestBarrier) {
          stest << "The reservoir size provided is too low, corrected to zero." << endl;
          grainLoc = 0;
        }
        else {
          grainLoc += lowestBarrier;
        }
      }

    }

    return true;
  }

}//namespace
