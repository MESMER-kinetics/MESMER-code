#include "MesmerILT.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  MesmerILT theMesmerILT("Mesmer ILT");
  //************************************************************

  /*-------------------------------------------------//
   Short note for variables & abbreviations in Mesmer: (REVERSIBLE association reaction)

   zpe_react:           zero-point energy of the reactant
   zpe_prodt:           zero-point energy of the product
   barri_hgt:           barrier Height of the forward path
   activ_ene:           activation energy (in addition to zpe_prodt)
   TS:                  transition state
   PES:                 potential energy surface
   A+B:                 molecules A and B
   A-B:                 complex formed by association of molecules A and B
            |
           /|\          TS
            |         *****       -\ barri_hgt                        -\
  potential |  A+B ***     *      -/               -\         activation\
   energy   |  (+)          *                        \          Energy   \
            |                *                        \                  /
            |                 *                       /                 /
            |                  *                     / zpe_react       /
            |               /-  **** A-B            /                -/
            |   zpe_prodt  /         (-)           /
           O|              \-                    -/
              ------------------------------------------------------------->
                             reaction coordinate
  PES

   Definition of a REVERSIBLE association reaction in Mesmer:

   1. A REVERSIBLE association reaction will be denoted as a forward reaction (normally proceeds toward right
      as in this PES).
   2. A REVERSIBLE association reaction can either be an exothermic reaction or an endothermic reaction.
   3. A reaction PES can change in different temperature, caused by rotational contribution to the total energy.
  //-------------------------------------------------*/

  bool MesmerILT::calculateMicroRateCoeffs(Reaction* p_reaction, vector<double> &cellKfmc, const MesmerEnv &mEnv)
  {
    //-----------------
    //starting variables block
    MesmerHP _ninf        = -0.25; // constraint: _ninf > -1.5
    double   _ainf        = p_reaction->get_PreExp();
    double   _einf        = p_reaction->get_ActivationEnergy();
    double   _tinf        = mEnv.temp;
    double   C_prime      = 3.24331e+20; // pow((2 * pi / (h * h)), 1.5)  not sure?? CHL
    //-----------------

    vector<CollidingMolecule *> unimolecularspecies;
    vector<CollidingMolecule *> bi_molecularspecies;
    p_reaction->get_unimolecularspecies(unimolecularspecies);
    p_reaction->get_bi_molecularspecies(bi_molecularspecies);

    CollidingMolecule * p_pdt1 = unimolecularspecies[0];
    CollidingMolecule * p_rct1 = bi_molecularspecies[0];
    Molecule          * p_rct2 = dynamic_cast<Molecule *>(bi_molecularspecies[1]); //not sure about this _2007_10_18__12_48_34_ CHL

    // Get molecular specific values
    double edg_a = static_cast<double>(p_pdt1->get_SpinMultiplicity());
    double edg_b = static_cast<double>(p_rct1->get_SpinMultiplicity());
    double edg_c = static_cast<double>(p_rct2->get_SpinMultiplicity());
    double ma = p_pdt1->getMass();
    double mb = p_rct1->getMass();
    double mc = p_rct2->getMass();

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> pdt1CellDOS(mEnv.MaxCell,0.0) ; // Cell density of states of      product molecule.
    vector<double> pdt1CellEne(mEnv.MaxCell,0.0) ; // Cell energies          of      product molecule.
    vector<double> rct1CellDOS(mEnv.MaxCell,0.0) ; // Cell density of states of 1st reactant molecule.
    p_pdt1->cellDensityOfStates(pdt1CellDOS, mEnv) ;
    p_pdt1->cellEnergies       (pdt1CellEne, mEnv) ;
    p_rct1->cellDensityOfStates(rct1CellDOS, mEnv) ;

    // Convolution of reactant vibrational DOS onto rotational DOS


    // Allocate space to hold Micro-canonical rate coefficients.
    cellKfmc.resize(mEnv.MaxCell); // no need to initialize

    int activ_ene = 0;
    if(IsNan(_einf))
    {
      stringstream errorMsg;
      errorMsg << "To use MesmerILT for reaction " << p_reaction->getName()
               << " the Activation Energy needs to be set.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    else{
      // Conversion of EINF from kcal.mol^-1 to cm^-1
      activ_ene = int(_einf * KcalPerMolToRC) ;
    }

    double _gamma = (MesmerHP) MesmerGamma(_ninf + 1.5);
    long double _ant = _ainf * C_prime * (edg_a * edg_b / edg_c) * pow( ( ma * mb / mc), 1.5 ) / _gamma;
    _ant /= (pow((_tinf * boltzmann_RCpK), _ninf));

    vector<double> work(mEnv.MaxCell);
    vector<double> conv(mEnv.MaxCell);

    double pwr     = _ninf + .5;
    for (int i = 0; i < mEnv.MaxCell; ++i) {
      work[i] = pow(pdt1CellEne[i], pwr);
    }

    DOSconvolution(conv, work, pdt1CellDOS);

    for (int i = 0; i < activ_ene; ++i)  cellKfmc[i] = 0.;
    for (int i = 0; i < (mEnv.MaxCell - activ_ene); ++i){
      cellKfmc[i + activ_ene] = _ant * conv[i] / rct1CellDOS[i + activ_ene];
    }

    return true;
  }
}//namespace
