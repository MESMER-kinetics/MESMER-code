#include "SimpleRRKM.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleRRKM theSimpleRRKM("Simple RRKM");
  //************************************************************

  bool SimpleRRKM::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc)
  {
    System* pSys = pReact->GetSys();  
    vector<CollidingMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    CollidingMolecule * pReactant = unimolecularspecies[0];

    TransitionState* pTS = pReact->get_TransitionState();
    if(!pTS)
    {
      stringstream errorMsg;
      errorMsg << "No transition state in Simple RRKM for reaction " << pReact->getName();
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }

    // Allocate space to hold Micro-canonical rate coefficients.
    kfmc.resize(pSys->MAXCell());

    // Initialize microcanoincal rate coefficients.

    int i, j ;
    for (i = 0 ; i < pSys->MAXCell() ; ++i ) {
        kfmc[i] = 0.0 ;
    }

    // Allocate some work space for density of states.

    vector<double> TScellDOS(pSys->MAXCell(),0.0) ; // Transistion state density of states.
    vector<double> cellDOS(pSys->MAXCell(),0.0) ; // Density of states of equilibrium molecule.

    // Extract densities of states from molecules.

    pReactant->cellDensityOfStates(&cellDOS[0]) ;
    pTS->cellDensityOfStates(&TScellDOS[0]) ;

    double SumOfStates  = 0.0 ;
    int thresholdEnergy = int((pTS->get_zpe() - pReactant->get_zpe()) * KcalPerMolToPerCm) ;
    for (i = thresholdEnergy, j = 0 ; i < pSys->MAXCell() ; ++i, ++j ) {

        // Integrate transition state density of states.

        SumOfStates += TScellDOS[j] ;

        // Calculate microcanonical rate coefficients using RRKM expression.

        kfmc[i] = SumOfStates / (plancksConst*cellDOS[i]) ;
    }
    return true;
  }

}//namespace
