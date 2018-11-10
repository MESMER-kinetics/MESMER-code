
//-------------------------------------------------------------------------------------------
//
// MolecularComponents.cpp
//
// Author: Struan Robertson
// Date:   29/Oct/2018
//
// This file contains the implementation of various generic methods for use along side 
// the MolecularComponents class.
//
//-------------------------------------------------------------------------------------------

#include "MolecularComponents.h"
#include "gDensityOfStates.h"

using namespace std;

namespace {

  // An anonymous namespace to contain methods to be used only within this file.

  // This method  performs a DOS convolution when one of the species is an atom 
  // and so takes advantage of its sparse density of states.
  bool FastAtomConvolution(vector<double>& atmCellDOS, vector<double>& otherCellDOS, vector<double>& rctsCellDOS) {

    rctsCellDOS.clear();
    rctsCellDOS.resize(otherCellDOS.size(), 0.0);

    for (size_t i(0); i < atmCellDOS.size(); i++) {
      if (atmCellDOS[i] > 0.0) {
        double multiplicity = atmCellDOS[i];
        for (size_t j(i), jj(0); j < rctsCellDOS.size(); j++, jj++) {
          rctsCellDOS[j] += multiplicity * otherCellDOS[jj];
        }
      }
    }

    return true;
  }

}

namespace mesmer
{

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS
  (gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, vector<double>& rctsCellDOS) {
    std::vector<double> rct1CellDOS;
    std::vector<double> rct2CellDOS;
    if (!pDOS1.getCellDensityOfStates(rct1CellDOS) || !pDOS2.getCellDensityOfStates(rct2CellDOS))
      return false;
    std::vector<double> rotConsts;
    // Check to see if one or other fragment is an atom. 
    // If so multiply by the atomic multiplicity.
    if (pDOS1.get_rotConsts(rotConsts) == ATOMIC) {
      FastAtomConvolution(rct1CellDOS, rct2CellDOS, rctsCellDOS);
    }
    else if (pDOS2.get_rotConsts(rotConsts) == ATOMIC) {
      FastAtomConvolution(rct2CellDOS, rct1CellDOS, rctsCellDOS);
    }
    else {
      FastLaplaceConvolution(rct1CellDOS, rct2CellDOS, rctsCellDOS);
    }
    return true;
  }

}