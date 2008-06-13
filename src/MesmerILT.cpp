#include "MesmerILT.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  MesmerILT theMesmerILT("Mesmer ILT");
  //************************************************************

  //-------------------------------------------------
  //   Short note for variables & abbreviations in Mesmer: (REVERSIBLE association reaction)
  //
  //   zpe_react:           zero-point energy of the reactant
  //   zpe_prodt:           zero-point energy of the product
  //   barri_hgt:           Barrier height (equals to theoretical calculated threshold energy)
  //   activ_ene:           activation energy (experimental value)
  //   TS:                  transition state
  //   PES:                 potential energy surface
  //   A+B:                 molecules A and B
  //   A-B:                 complex formed by association of molecules A and B
  //            |
  //           /|\          TS
  //            |         *****       -\ barri_hgt                        -\
  //  potential |  A+B ***     *      -/               -\         activation\
  //   energy   |  (+)          *                        \          Energy   \
  //            |                *                        \                  /
  //            |                 *                       /                 /
  //            |                  *                     / zpe_react       /
  //            |               /-  **** A-B            /                -/
  //            |   zpe_prodt  /         (-)           /
  //           O|              \-                    -/
  //              ------------------------------------------------------------->
  //                             reaction coordinate
  //  PES
  //
  //   Definition of a REVERSIBLE association reaction in Mesmer:
  //
  //   1. A REVERSIBLE association reaction is going forward when the reaction is going from left to right in this
  //      potential energy surface.
  //   2. A reaction PES can change in different temperature, caused by rotational contribution to the total energy.
  //-------------------------------------------------

  bool MesmerILT::calculateMicroRateCoeffs(Reaction* pReact)
  {
    //-----------------
    //starting variables block
    const double n_infinity   = pReact->get_NInf(); // constraint: n_infinity > -1.5
    const double t_infinity   = 1. / (boltzmann_RCpK * pReact->getEnv().beta);
    const double preExp       = pReact->get_PreExp();
    const int    activ_ene    = int(pReact->get_ThresholdEnergy());
    // double tp_C = 3.24331e+20; // Defined in Constant.h, constant used in the translational partition function
    //-----------------

    SuperMolecule*              p_rcts = NULL;
    p_rcts = pReact->get_bi_molecularspecies();
    if (!p_rcts){
      cerr << "Not a valid bi-molecularspecies";
      return false;
    }

    vector<ModelledMolecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);

    ModelledMolecule*  p_pdt1 = unimolecularspecies[0];
    ModelledMolecule*  p_rct1 = p_rcts->getMember1();
    ModelledMolecule*  p_rct2 = p_rcts->getMember2();

    // Get molecular specific values
    const double edg_a = static_cast<double>(p_rct1->getSpinMultiplicity());
    const double edg_b = static_cast<double>(p_rct2->getSpinMultiplicity());
    const double edg_c = static_cast<double>(p_pdt1->getSpinMultiplicity());
    const double ma = p_rct1->getMass();
    const double mb = p_rct2->getMass();
    const double mc = p_pdt1->getMass();
    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctsCellEne; // Cell energies          of      product molecule.
    vector<double> rctsCellDOS; // Convoluted cell density of states of reactants.

    p_rcts->getCellEnergies       (rctsCellEne) ;
    p_rcts->getCellDensityOfStates(rctsCellDOS) ;

    // Allocate space to hold microcanonical rate coefficients for dissociation.
    vector<double>& TSFlux = pReact->get_CellFlux();
    TSFlux.clear();
    TSFlux.resize(MaximumCell, 0.0); 

    const double _gamma = MesmerGamma(n_infinity + 1.5);
    double _ant = preExp * tp_C * (edg_a * edg_b / edg_c) * pow( ( ma * mb / mc), 1.5 ) / _gamma;
    _ant /= (pow((t_infinity * boltzmann_RCpK), n_infinity));

    vector<double> work(MaximumCell);
    vector<double> conv(MaximumCell);

    double pwr     = n_infinity + .5;
    for (int i = 0; i < MaximumCell; ++i) {
      work[i] = pow(rctsCellDOS[i], pwr);
    }

    DOSconvolution(work, rctsCellDOS, conv);

    for (int i = activ_ene; i < MaximumCell; ++i){
      TSFlux[i] = _ant * conv[i - activ_ene];
    }

    // the flux bottom energy is equal to the well bottom of the product
    pReact->setCellFluxBottom(p_rcts->get_relative_ZPE());

    cinfo << "ILT calculation completed" << endl;

    return true;
  }

}//namespace
