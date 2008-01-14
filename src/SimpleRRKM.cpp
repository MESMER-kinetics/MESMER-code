#include "SimpleRRKM.h"


using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleRRKM theSimpleRRKM("Simple RRKM");
  //************************************************************

  bool SimpleRRKM::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &cellKfmc, const MesmerEnv &Env)
  {
    vector<ModelledMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    ModelledMolecule * pReactant = unimolecularspecies[0];

    TransitionState* pTS = pReact->get_TransitionState();
    if(!pTS)
    {
      stringstream errorMsg;
      errorMsg << "No transition state in Simple RRKM for reaction " << pReact->getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }

    // Allocate space to hold Micro-canonical rate coefficients.
    cellKfmc.resize(Env.MaxCell);

    // Initialize microcanoincal rate coefficients.

    int i, j ;
    for (i = 0 ; i < Env.MaxCell ; ++i ) {
        cellKfmc[i] = 0.0 ;
    }

    // Allocate some work space for density of states.

    vector<double> TScellDOS(Env.MaxCell,0.0) ; // Transistion state density of states.
    vector<double> cellDOS(Env.MaxCell,0.0) ; // Density of states of equilibrium molecule.

    // Extract densities of states from molecules.

    pReactant->getCellDensityOfStates(cellDOS) ;
    pTS->getCellDensityOfStates(TScellDOS) ;

    double SumOfStates  = 0.0 ;
    int thresholdEnergy = int((pTS->get_zpe() - pReactant->get_zpe()) * KcalPerMolToRC) ;
    for (i = thresholdEnergy, j = 0 ; i < Env.MaxCell ; ++i, ++j ) {

        // Integrate transition state density of states.

        SumOfStates += TScellDOS[j] ;

        // Calculate microcanonical rate coefficients using RRKM expression.

        cellKfmc[i] = SumOfStates * SpeedOfLight_cm / cellDOS[i];
    }
    return true;
  }

}//namespace
