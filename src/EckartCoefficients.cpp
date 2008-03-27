#include "EckartCoefficients.h"

using namespace Constants;

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  EckartCoefficients theEckartCoefficients("Eckart");
  //************************************************************

  bool EckartCoefficients::calculateTunnelingCoeffs(Reaction* pReact){
    std::vector<ModelledMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    ModelledMolecule * pReactant = unimolecularspecies[0];
    ModelledMolecule * p_Product = unimolecularspecies[1];
    TransitionState  * p_TransitionState = pReact->get_TransitionState();
    double TC = p_TransitionState->getClassicalEnergy();
    double TZ = p_TransitionState->get_zpe();

    double V0 = (TC - pReactant->getClassicalEnergy()) * kJPerMolInRC;
    double V1 = (TC - p_Product->getClassicalEnergy()) * kJPerMolInRC;
    double barrier0 = (TZ - pReactant->get_zpe()) * kJPerMolInRC;
    double barrier1 = (TZ - p_Product->get_zpe()) * kJPerMolInRC;
    double imFreq = pReact->m_TransitionState->get_ImFreq() * SpeedOfLight_cm; // convert im_freq from cm-1 -> Hz

    const MesmerEnv& mEnv = pReactant->getEnv();
    int MaxCell = mEnv.MaxCell;
    pReact->m_CellTunneling.resize(MaxCell);

    for(int i = 0; i < MaxCell; ++i){
      double j = double(i);
      double E = j - barrier0;
      if ((E + barrier1)<0.0){
        pReact->m_CellTunneling[i] = 0.0;
      }
      else
      {
        double a = (4.0 * M_PI * SpeedOfLight_cm /imFreq)* sqrt(E + V0) * 1.0/(1.0/sqrt(V0) + 1.0/sqrt(V1));
        double b = (4.0 * M_PI * SpeedOfLight_cm /imFreq)* sqrt(E + V1) * 1.0/(1.0/sqrt(V0) + 1.0/sqrt(V1));
        double c = 2.0 * M_PI * sqrt(V0 * V1 / (pow(imFreq/SpeedOfLight_cm ,2.0)) - 1.0/16.0);
        pReact->m_CellTunneling[i] = (sinh(a) * sinh(b)) / (pow(sinh((a+b)/2.0),2.0) + pow(cosh(c),2.0));
        if(IsNan(pReact->m_CellTunneling[i])) pReact->m_CellTunneling[i] = 0.0;  // if statement to avoid nan at small values of E
      }
    }
    return true;
  }
}//namespace

