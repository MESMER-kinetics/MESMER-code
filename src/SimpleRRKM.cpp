#include "SimpleRRKM.h"


using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleRRKM theSimpleRRKM("Simple RRKM");
  //************************************************************

  bool SimpleRRKM::calculateMicroRateCoeffs(Reaction* pReact)
  {
    SuperMolecule*              p_rcts = NULL;
    p_rcts = pReact->get_bi_molecularspecies();

    vector<ModelledMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);

    ModelledMolecule * p_rct = NULL;

    if (!p_rcts){
      p_rct = unimolecularspecies[0];
    }
    else{
      p_rct = p_rcts;
    }

    TransitionState* pTS = pReact->get_TransitionState();
    if(!pTS)
    {
      cerr << "Lack of transition state in reaction " << pReact->getName() << " for Simple RRKM" << endl;
      return false;
    }
    int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold Micro-canonical rate coefficients.
    pReact->m_CellKfmc.resize(MaximumCell, 0.0);

    // Allocate some work space for density of states.

    vector<double> TScellDOS; // Transistion state density of states.
    vector<double> rctCellDOS; // Density of states of equilibrium molecule.

    // Extract densities of states from molecules.
    p_rct->getCellDensityOfStates(rctCellDOS) ;
    pTS->getCellDensityOfStates(TScellDOS) ;

    int thresholdEnergy = int(pReact->get_ActivationEnergy() * kJPerMolInRC) ;
    
    if (pReact->m_pTunnelingCalculator) { // with tunneling
      pReact->m_pTunnelingCalculator->calculateTunnelingCoeffs(pReact);

      for (int i = 0; i < MaximumCell ; ++i) {
        // Integrate transition state density of states.
        double SumOfStates = 0.0;
        for (int j = 0; j <= i; ++j) SumOfStates += pReact->m_CellTunneling[i-j] * TScellDOS[j];

        // Calculate microcanonical rate coefficients using RRKM expression with tunneling correction.
        pReact->m_CellKfmc[i] = SumOfStates * SpeedOfLight_cm / rctCellDOS[i];
      }
    }
    else{ // without tunneling
      double SumOfStates = 0.0;
      for (int i = thresholdEnergy, j = 0 ; i < MaximumCell ; ++i, ++j ) {
        // Integrate transition state density of states.
        SumOfStates += TScellDOS[j];

        // Calculate microcanonical rate coefficients using RRKM expression.
        pReact->m_CellKfmc[i] = SumOfStates * SpeedOfLight_cm / rctCellDOS[i];
      }
    }

    return true;
  }

}//namespace
