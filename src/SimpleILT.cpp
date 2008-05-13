#include "SimpleILT.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleILT theSimpleILT("Simple ILT");
  //************************************************************

  //-----
  //short note for variables:
  //dh00:
  //


  bool SimpleILT::calculateMicroRateCoeffs(Reaction* pReact)
  {
    SuperMolecule*              p_rcts = NULL;
    p_rcts = pReact->get_bi_molecularspecies();
    if (!p_rcts){
      cerr << "Not a valid bi-molecularspecies";
      return false;
    }

    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& TSFlux = pReact->get_CellFlux();
    TSFlux.resize(MaximumCell, 0.0);

    // Allocate some work space for density of states.

    vector<double> rctsCellDOS; // Density of states of equilibrim molecule.

    // Extract densities of states from molecules.

    p_rcts->getCellDensityOfStates(rctsCellDOS) ;

    // Conversion of EINF from kcal.mol^-1 to cm^-1

    const int nEinf = int(pReact->get_ThresholdEnergy()) ;
    const double preExp = pReact->get_PreExp();

    // Calculate microcanonical rate coefficients using simple ILT expression.

    for (int i = nEinf; i < MaximumCell ; ++i ) {
      TSFlux[i] = preExp * rctsCellDOS[i-nEinf];
    }

    // the flux bottom energy is equal to the well bottom of the source term
    pReact->setCellFluxBottom(p_rcts->get_relative_ZPE());

    return true;
  }

}//namespace
