#include "SimpleILT.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleILT theSimpleILT("Simple ILT");
  //************************************************************

  //
  // This method calculates the reaction flux by Laplace inversion
  // of the Arrhenius equation for a reaction proceeding in the 
  // unimolecular direction.
  //

  bool SimpleILT::calculateMicroRateCoeffs(Reaction* pReact)
  {
    vector<ModelledMolecule *> Isomers ;
    int nIsomers = pReact->get_unimolecularspecies(Isomers) ;  
    ModelledMolecule *p_rcts = Isomers[0] ;

    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& TSFlux = pReact->get_CellFlux();
    TSFlux.clear();
    TSFlux.resize(MaximumCell, 0.0);

    // Allocate some work space for and obtain density of states of the unimolecuar reactant.

    vector<double> rctsCellDOS; 
    p_rcts->getCellDensityOfStates(rctsCellDOS) ;

    // Obtain the Arrhenius parameters.

    const int nEinf = int(pReact->get_ThresholdEnergy()) ;  //<-- SHR: This should be activation energy NOT threshold energy!!
    const double preExp = pReact->get_PreExp();

    // Calculate microcanonical rate coefficients using simple ILT expression.

    for (int i = nEinf; i < MaximumCell ; ++i ) {
      TSFlux[i] = preExp * rctsCellDOS[i-nEinf];
    }

    // the flux bottom energy is equal to the well bottom of the source term
    pReact->setCellFluxBottom(p_rcts->get_zpe() - pReact->getEnv().EMin);

    return true;
  }

}//namespace
