#include "ClassicalRotor.h"

using namespace std;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  ClassicalRotor theClassicalRotor("Classical rotors");
  //************************************************************

  // Provide a function to define particular counts of the DOS of a molecule.
  bool ClassicalRotor::countCellDOS(ModelledMolecule* pMol)
  {
    vector<double> VibFreq ; 
    pMol->get_VibFreq(VibFreq) ;

    // times the scale factor
    for (vector<double>::size_type i = 0; i < VibFreq.size(); ++i){
      VibFreq[i] *= pMol->get_scaleFactor();
    }

    const int MaximumCell = pMol->getEnv().MaxCell;
    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0) ;

    //
    // Initialize density of states array using calculated rotational
    // density of state.
    //

    //From inverse Laplace transform of rotors
    vector<double> rotConst; int rotorType = pMol->get_rotConsts(rotConst);
    double sym = pMol->get_Sym();
    double qele = pMol->getSpinMultiplicity();
    double cnt = 0.;

    for (int i = 0 ; i < MaximumCell ; ++i ) {
      cellEne[i] = double(i) + 0.5 ;
    }

    switch (rotorType){
    case 2: //3-D symmetric/asymmetric/spherical top
      cnt = qele * sqrt(4./(rotConst[0] * rotConst[1] * rotConst[2]))/sym ;
      for (int i = 0 ; i < MaximumCell ; ++i ) 
        cellDOS[i] = cnt*sqrt(cellEne[i]) ;
      break;
    case 0: //2-D linear
      cnt = qele / (rotConst[0] * sym);
      for (int i = 0 ; i < MaximumCell ; ++i ) 
        cellDOS[i] = cnt ;
      break;
    default:
      cnt = 0.;
      for (int i = 0 ; i < MaximumCell ; ++i ) 
        cellDOS[i] = cnt ;
    }

    // Implementation of the Beyer-Swinehart algorithm.
    Beyer_Swinehart(VibFreq, cellDOS);

    //electronic degeneracy
    vector<double> eleExc;
    pMol->getEleExcitation(eleExc);
    if (!eleExc.empty()){
      for (int j = 0; j < int(eleExc.size()) ; ++j){
        int iele = static_cast<int>(eleExc[j]);
        for (int i = (MaximumCell - 1); i >= (iele - 1); --i){
          cellDOS[i] += cellDOS[i - iele + 1];
        }
      }
    }

    pMol->setCellDensityOfStates(cellDOS) ;

    return true;
  }

}//namespace
