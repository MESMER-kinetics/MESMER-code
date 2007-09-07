#include <float.h>
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


  bool SimpleILT::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc)
  {
    System* pSys = pReact->GetSys();

    vector<CollidingMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    CollidingMolecule * pReactant = unimolecularspecies[0];
    if(IsNan(pReact->get_ActivationEnergy()))
    {
      stringstream errorMsg;
      errorMsg << "To use SimpleILT for reaction " << pReact->getName()
               << " the Activation Energy needs to be set.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    // Allocate space to hold Micro-canonical rate coefficients.
    kfmc.resize(pSys->MAXCell());

    // Initialize microcanoincal rate coefficients.

    int i ;
    for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
        kfmc[i] = 0.0 ;
    }

    // Allocate some work space for density of states.

    vector<double> cellDOS(pSys->MAXCell(),0.0) ; // Density of states of equilibrim molecule.

    // Extract densities of states from molecules.

    pReactant->cellDensityOfStates(&cellDOS[0]) ;

    // Conversion of EINF from kcal.mol^-1 to cm^-1

    int nEinf = int(pReact->get_ActivationEnergy()*KcalPerMolToPerCm) ;

    // Calculate microcanonical rate coefficients using simple ILT expression.

    for (i = nEinf ; i < pSys->MAXCell() ; ++i ) {
      kfmc[i] = pReact->get_PreExp()*cellDOS[i-nEinf] / cellDOS[i] ;
    }

    return true;
  }

}//namespace
