#include "HinderedRotorInterpolation.h"

using namespace std;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  HinderedRotorInterpolation theHinderedRotorInterpolation("Hindered rotor interpolation");
  //************************************************************

  // Provide a function to define particular counts of the DOS of a molecule.
  bool HinderedRotorInterpolation::countCellDOS(gDensityOfStates* pDOS, int MaximumCell, PersistPtr ppDOSC)
  {
    // get data from XML (examples)
    // PersistPtr pLocal = ppDOSC->XmlMoveTo("local", optional); // local minimum
    // ppDOSC
    //electronic degeneracy
    vector<double> eleExc;

    const char* txt= ppDOSC->XmlReadProperty("me:electronicExcitation", optional);
    if(txt) {
      istringstream idata(txt);
      double _iele; 
      while (idata >> _iele) eleExc.push_back(_iele);
    }

    vector<double> cellDOS;
    pDOS->getCellDensityOfStates(cellDOS) ; // retreive the DOS vector.

    if (!eleExc.empty()){
      for (int j(0); j < int(eleExc.size()) ; ++j){
        int iele = static_cast<int>(eleExc[j]);
        for (int i = (MaximumCell - 1); i >= (iele - 1); --i){
          cellDOS[i] += cellDOS[i - iele + 1];
        }
      }
    }

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;
  }

}//namespace
