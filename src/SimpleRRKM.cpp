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
    TransitionState* pTS = pReact->get_TransitionState();
    if(!pTS)
    {
      cerr << "Lack of transition state in reaction " << pReact->getName() << " for Simple RRKM" << endl;
      return false;
    }

    double ThresholdEnergy = pReact->get_relative_TSZPE() - pReact->get_relative_rctZPE();
    pReact->set_ThresholdEnergy(ThresholdEnergy);

    // Allocate some work space for density of states.
    vector<double> TScellDOS; // Transistion state density of states.
    pTS->getCellDensityOfStates(TScellDOS) ; // Extract densities of states from molecules.

    // get MaxCell from MesmerEnv structure via Reaction class
    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& TSFlux = pReact->get_CellFlux();
    TSFlux.clear();
    TSFlux.resize(MaximumCell, 0.0);

    if (pReact->thereIsTunnelling()) { // with tunneling
      int HeatOfReaction = pReact->getHeatOfReactionInt();
      const int TunnelingStart = (HeatOfReaction > 0) ? int(HeatOfReaction) : 0;

      vector<double> TunnelingProbability;
      pReact->calculateCellTunnelingCoeffs(TunnelingProbability);

      vector<double> ConvolvedSumOfStates;
      FastLaplaceConvolution(TScellDOS, TunnelingProbability, ConvolvedSumOfStates); // FFT convolution
//      Convolution(TScellDOS, TunnelingProbability, ConvolvedSumOfStates); // standard convolution

      for (int i = TunnelingStart; i < MaximumCell; ++i) {                       // Calculate TSflux using RRKM 
        TSFlux[i-TunnelingStart] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm; // with tunneling correction
      }

      // the flux bottom energy is equal to the ZPE of the higher well
      if (TunnelingStart > 0) {
        pReact->setCellFluxBottom(pReact->get_relative_pdtZPE());
      }
      else{
        pReact->setCellFluxBottom(pReact->get_relative_rctZPE());
      }

    }
    else{ // without tunneling
      double SumOfStates = 0.0;
      for (int i = 0 ; i < MaximumCell ; ++i) {
        // Integrate transition state density of states.
        SumOfStates += TScellDOS[i];

        // Calculate microcanonical rate coefficients using RRKM expression.
        TSFlux[i] = SumOfStates * SpeedOfLight_in_cm;
      }

      // the flux bottom energy is equal to the ZPE of the transition state
      pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + pReact->get_ThresholdEnergy());
    }


    return true;
  }

}//namespace
