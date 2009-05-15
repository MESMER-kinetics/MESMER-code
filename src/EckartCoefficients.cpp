#include "EckartCoefficients.h"

using namespace Constants;

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  EckartCoefficients theEckartCoefficients("Eckart");
  //************************************************************



  bool EckartCoefficients::calculateCellTunnelingCoeffs(Reaction* pReact, vector<double>& TunnelingProbability){
    std::vector<Molecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    Molecule * pReactant = unimolecularspecies[0];
    Molecule * p_Product = unimolecularspecies[1];
    Molecule * p_TransitionState = pReact->get_TransitionState();

    // PLEASE CHECK THIS SECTION FOR NUMBERS CHRIS
    const double TZPE = p_TransitionState->getDOS().get_zpe();
    const double oz1 = pReactant->getDOS().get_zpe();
    const double oz2 = p_Product->getDOS().get_zpe();
    const double ZPE0 = TZPE - oz1;
    const double ZPE1 = TZPE - oz2;
    const double diff = ZPE0 - ZPE1;
    const double odiff = pReactant->getDOS().get_zpe() - p_Product->getDOS().get_zpe();
    // PLEASE CHECK THIS SECTION FOR NUMBERS CHRIS

    //TC is the classical energy of the TS
    const double TC = p_TransitionState->getDOS().getClassicalEnergy();
    //V0 & V1 are the classical barrier heights in the forward/reverse directions
    const double V0 = TC - pReactant->getDOS().getClassicalEnergy();
    const double V1 = TC - p_Product->getDOS().getClassicalEnergy();

    //TZ is the zpe of the TS
    const double TZ = pReact->get_relative_TSZPE();
    //barrier0 & barrier1 are the zpe corrected barrier heights in the forward/reverse directions
    const int barrier0 = int(TZ) - int(pReact->get_relative_rctZPE());
    const int barrier1 = int(TZ) - int(pReact->get_relative_pdtZPE());

    //imFreq is the imaginary frequency of the TS
    const double imFreq = pReact->get_TSImFreq() * SpeedOfLight_in_cm; // convert im_freq from cm-1 -> Hz

    //get properties of vectors in which to include transmission coefficients
    const int MaximumCell = pReactant->getEnv().MaxCell;
    TunnelingProbability.clear();
    TunnelingProbability.resize(MaximumCell);

    //set transmission coefficients to 0 where no tunneling is possible;
    //where tunneling may occur, the transmission coefficients are calculated using, a, b, & c, 
    //for a 1d eckart barrier as described by W.H. Miller, JACS, 101(23), 1979

    for(int i = 0; i < MaximumCell; ++i){
      int E = i - barrier0;
      if ((E + barrier1) < 0){
        TunnelingProbability[i] = 0.0;
      }
      else
      {
        double a = (4.0 * M_PI * SpeedOfLight_in_cm /imFreq)* sqrt(E + V0) * 1.0/(1.0/sqrt(V0) + 1.0/sqrt(V1));
        double b = (4.0 * M_PI * SpeedOfLight_in_cm /imFreq)* sqrt(E + V1) * 1.0/(1.0/sqrt(V0) + 1.0/sqrt(V1));
        double c = 2.0 * M_PI * sqrt(V0 * V1 / (pow(imFreq/SpeedOfLight_in_cm ,2.0)) - 1.0/16.0);
        TunnelingProbability[i] = (sinh(a) * sinh(b)) / (pow(sinh((a+b)/2.0),2.0) + pow(cosh(c),2.0));
        // following if statement to avoid nan at small values of E
        if(IsNan(TunnelingProbability[i])) TunnelingProbability[i] = 0.0;  
      }
    }

    if (pReact->getFlags().TunnellingCoeffEnabled){
      ctest << "\nTunneling coefficients for: " << pReact->getName() 
        << "\nV0 = " << V0 << ", V1 = " << V1 
        << ", barrier0 = " << barrier0 << ", barrier1 = " << barrier1 
        << ", imFreq = " << imFreq << "\n{\n";
      for(int i = 0; i < MaximumCell; ++i){
        ctest << TunnelingProbability[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }
}//namespace

