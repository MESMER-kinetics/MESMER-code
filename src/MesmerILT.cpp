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
    //additional variables
    System* pSys    = pReact->GetSys();
    double _ninf    = .0;
    int MaximumCell = pSys->MAXCell();
    //-----------------

    vector<CollidingMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    CollidingMolecule * pReactant = unimolecularspecies[0];
    if(IsNan(pReact->get_ActivationEnergy()))
    {
      stringstream errorMsg;
      errorMsg << "To use MesmerILT for reaction " << pReact->getName()
               << " the Activation Energy needs to be set.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    // Allocate space to hold Micro-canonical rate coefficients.
    cellKfmc.resize(MaximumCell);

    // Initialize microcanoincal rate coefficients.
    for (int i = 0 ; i < MaximumCell ; ++i ) cellKfmc[i] = 0.0 ;

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> cellDOS(MaximumCell,0.0) ; // Density of states of equilibrim molecule.
    pReactant->cellDensityOfStates(cellDOS) ;

    // Conversion of EINF from kcal.mol^-1 to cm^-1
    int activationEnergy = int(pReact->get_ActivationEnergy() * KcalPerMolToRC) ;

    // Calculate microcanonical rate coefficients using ILT expression.
    for (int i = activationEnergy ; i < MaximumCell ; ++i ) {
      cellKfmc[i] = pReact->get_PreExp() * cellDOS[i - activationEnergy] / cellDOS[i] ;
    }

    //-----------------
    //starting block
    double C_prime = 3.24331e+20; // powl((2 * pi / (h * h)), 1.5)  not sure??
    double pwr     = _ninf + .5;
    //-----------------


    double _gamma = (long double) MesmerGamma(_ninf + 1.5);
    long double _ant = _ainf * C_prime * (edg_a * edg_b / edg_c) * powl( ( ma * mb / mc), 1.5 ) / _gamma;
    _ant /= (powl((tinf * boltzmann_RCpK), _ninf));


    vector<double> work1(MaximumCell);
    vector<double> KCell(MaximumCell);
    vector<double> conv (MaximumCell);

    for (int i = 0; i < MaximumCell; ++i) {
      work1[i] = powl(m_cellEne[i], _ninf + .5);
      KCell[i] = cellDOS[i];
    }

    for (int i = 0; i < MaximumCell; ++i){
      conv[i] = 0.;
      for (int j = 1; j <= i; ++j)
        conv[i] += work1[j] * KCell[i + 1 - j];
    }

    for (int i = 0; i < activationEnergy; ++i)  KCell[i] = 0.;

    for (int i = 0; i < (MaximumCell - activationEnergy); ++i){
      KCell[i + activationEnergy] = _ant * conv[i] / tnr[i + activationEnergy];
    }

    int idx1 = 0; int idx2 = 0;
    for (int i = 0; i < ngrn; ++i) {
      int idx3 = idx1;
      smt = 0.;
      for (int j = 0; j < gsz; ++j) {
        ++idx1;
        smt += tnr[idx1];
      }
      if (smt > .0) {
        ++idx2;
        tn2r[idx2] = smt;
        smat = 0.;
        for (j = 0; j < gsz; ++j) {
          ++idx3;
          smat += (KCell[idx3]*tnr[idx3])/smt;
        }
        kgrn[idx2] = smat;
      }
    }


    return true;
  }

}//namespace
