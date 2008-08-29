#include "DensityOfStates.h"
#include "Molecule.h"
//
// Test the forward microcanonical rate coefficients.
//
namespace mesmer
{

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool DensityOfStatesCalculator::countDimerCellDOS(ModelledMolecule* p_mol1,
                                                    ModelledMolecule* p_mol2, 
                                                    std::vector<double>& rctsCellDOS){
    std::vector<double> rct1CellDOS;
    std::vector<double> rct2CellDOS;
    p_mol1->getCellDensityOfStates(rct1CellDOS);
    p_mol2->getCellDensityOfStates(rct2CellDOS);
    std::vector<double> rotConsts;
    if (p_mol1->get_rotConsts(rotConsts) == -4){
      rctsCellDOS = rct2CellDOS;
    }
    else if(p_mol2->get_rotConsts(rotConsts) == -4){
      rctsCellDOS = rct1CellDOS;
    }
    else{
      FastLaplaceConvolution(rct1CellDOS, rct2CellDOS, rctsCellDOS);
    }
    return true;
  }

}//namespace
