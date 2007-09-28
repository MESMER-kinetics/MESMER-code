#include <float.h>
#include "MesmerILT.h"
#include "MesmerMath.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  MesmerILT theMesmerILT("Mesmer ILT");
  //************************************************************

  //-----
  //short note for variables: (association reaction)
  //_zpe_reactants:    zero-point energy of the reactants
  // barrierHeight:    barrier Height of the forward path
  // activationEnergy: activation energy
  //
  //       *****       -\ barrierHeight                           -\
  //A+B ***     *      -/                 -\                        \
  //             *                          \                        \ activationEnergy
  //              *                          \_zpe_reactants         /
  //               *                         /                      /
  //                *                       /                      /
  //                 **** A-B             -/                     -/
  //PES


  bool MesmerILT::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &cellKfmc)
  {

    //-----------------
    //starting variables block
    System* pSys        = pReact->GetSys();
    MesmerHP _ninf        = -0.25; // constraint: _ninf > -1.5
    double _ainf        = 2.03e-12;
    double _tinf        = 300.;  // default temperature
    double _einf        = pReact->get_ActivationEnergy();
    int MaximumCell     = pSys->MAXCell();
    double C_prime      = 3.24331e+20; // pow((2 * pi / (h * h)), 1.5)  not sure??
    //-----------------

    vector<CollidingMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    CollidingMolecule * pReactant = unimolecularspecies[0];
    CollidingMolecule * pProduct  = unimolecularspecies[1]; //not sure about this _2007_09_26__15_38_30_ Chi-Hsiu Liang

  int activationEnergy = 0;
  if(IsNan(_einf))
    {
      stringstream errorMsg;
      errorMsg << "To use MesmerILT for reaction " << pReact->getName()
               << " the Activation Energy needs to be set.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    else{
      // Conversion of EINF from kcal.mol^-1 to cm^-1
      activationEnergy = int(_einf * KcalPerMolToRC) ;
    }
    
    // Allocate space to hold Micro-canonical rate coefficients.
    cellKfmc.resize(MaximumCell); // no needs to initialize

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctCellDOS(MaximumCell,0.0) ; // Cell density of states of reactant molecule.
    vector<double> pdtCellDOS(MaximumCell,0.0) ; // Cell density of states of product molecule.
    pReactant->cellDensityOfStates(rctCellDOS) ;
    pProduct->cellDensityOfStates(pdtCellDOS) ;

    double _gamma = (MesmerHP) MesmerGamma(_ninf + 1.5);
    long double _ant = _ainf * C_prime * (edg_a * edg_b / edg_c) * powl( ( ma * mb / mc), 1.5 ) / _gamma;
    _ant /= (powl((_tinf * boltzmann_RCpK), _ninf));


    vector<double> work1(MaximumCell);
    vector<double> conv (MaximumCell);

    double pwr     = _ninf + .5;
    for (int i = 0; i < MaximumCell; ++i) {
      work1[i] = powl(m_cellEne[i], pwr);
      cellKfmc[i] = rctCellDOS[i];
    }

    //convolution
    for (int i = 0; i < MaximumCell; ++i){
      conv[i] = 0.;
      for (int j = 1; j <= i; ++j)
        conv[i] += work1[j] * cellKfmc[i + 1 - j];
    }

    for (int i = 0; i < activationEnergy; ++i)  cellKfmc[i] = 0.;

    for (int i = 0; i < (MaximumCell - activationEnergy); ++i){
      cellKfmc[i + activationEnergy] = _ant * conv[i] / pdtCellDOS[i + activationEnergy];
    }

    return true;
  }
}//namespace
