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
  // avtivationEnergy: activation energy
  //
  //       *****       -\ barrierHeight                           -\
  //A+B ***     *      -/                 -\                        \
  //             *                          \                        \ avtivationEnergy
  //              *                          \_zpe_reactants         /
  //               *                         /                      /
  //                ****                    /                      /
  //                A-B                   -/                     -/
  //PES


  bool MesmerILT::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &kfmc)
  {
    System* pSys = pReact->GetSys();

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

    int barrierHeight = int(pReact->get_ActivationEnergy()*KcalPerMolToRC) ;

    // Calculate microcanonical rate coefficients using ILT expression.

    for (i = barrierHeight ; i < pSys->MAXCell() ; ++i ) {
      kfmc[i] = pReact->get_PreExp() * cellDOS[i - barrierHeight] / cellDOS[i] ;
    }

    //-----------------
    //starting block
    long double C_prime = 3.24331e+20; // powl((2 * pi / (h * h)), 1.5)
    long double pwr     = _ninf + .5;
    //-----------------

    long double _gamma = (long double) MesmerGamma(_ninf + 1.5);
    long double _ant = _ainf * C_prime * (edg_a * edg_b / edg_c) * powl( ( ma * mb / mc), 1.5 ) / _gamma;
    _ant /= (powl((tinf * boltzmann_RCpK), _ninf));
    avtivationEnergy = (int) (_zpe_reactants + barrierHeight);

    int MaximumCell = pSys->MAXCell();

    //--
    vector<double> work1(MaximumCell * 2);
    vector<double> KCell(MaximumCell * 2);
    vector<double> Corr (MaximumCell    );

    for (int i = 0; i < MaximumCell; ++i) {
      work1[i              ] = powl(m_cellEne[i], _ninf + .5);
      work1[i + MaximumCell] = .0;
      KCell[i              ] = cellDOS[i];
      KCell[i + MaximumCell] = .0;
    }

    for (i = 0; i < MaximumCell; ++i){
      Corr[i] = 0.;
      for (j=1; j<=i; j++) Corr[i] += work1[j]*KCell[i + 1 - j];
    }

    for (i = 0; i < MaximumCell; ++i){
      work1[i] = corr [i];
    }

    for (i = 0; i < avtivationEnergy; ++i){
      KCell[i] = 0.l;
    }

    for (i = 0; i < (MaximumCell - avtivationEnergy); ++i){
      KCell[i + avtivationEnergy] = ant * work1[i] / tnr[i + avtivationEnergy];
    }

    int idx1 = 0; int idx2 = 0;
    for (i = 0; i < ngrn; ++i) {
      int idx3 = idx1;
      smt = 0.0l;
      for (j = 0; j < gsz; ++j) {
        ++idx1;
        smt += tnr[idx1];
      }
      if (smt > .0) {
        ++idx2;
        tn2r[idx2] = smt;
        smat = 0.0l;
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
