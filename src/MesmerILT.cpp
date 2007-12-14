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

  bool MesmerILT::calculateMicroRateCoeffs(Reaction* pReact, vector<double> &cellKfmc, const MesmerEnv &mEnv)
  {
    //-----------------
    //starting variables block
    MesmerHP _ninf        = pReact->get_NInf(); // constraint: _ninf > -1.5
    double   _ainf        = pReact->get_PreExp();
    double   _einf        = pReact->get_ActivationEnergy();
    double   _tinf        = 1. / (boltzmann_RCpK * mEnv.beta);
    double   C_prime      = 3.24331e+20; // pow((2 * pi / (h * h)), 1.5)  not sure?? CHL
    //-----------------

    SuperMolecule*              p_rcts = NULL;
    vector<ModelledMolecule *> unimolecularspecies;

    pReact->get_unimolecularspecies(unimolecularspecies);
    p_rcts = pReact->get_bi_molecularspecies();
    if (!p_rcts){
      stringstream errorMsg;
      errorMsg << "Not a valid bi-molecularspecies";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }

    ModelledMolecule*  p_pdt1 = unimolecularspecies[0];
    ModelledMolecule*  p_rct1 = p_rcts->getMember1();
    ModelledMolecule*  p_rct2 = p_rcts->getMember2();

    // Get molecular specific values
    double edg_a = static_cast<double>(p_pdt1->getSpinMultiplicity());
    double edg_b = static_cast<double>(p_rct1->getSpinMultiplicity());
    double edg_c = static_cast<double>(p_rct2->getSpinMultiplicity());
    double ma = p_pdt1->getMass();
    double mb = p_rct1->getMass();
    double mc = p_rct2->getMass();

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> pdt1CellDOS(mEnv.MaxCell,0.0) ; // Cell density of states of      product molecule.
    vector<double> pdt1CellEne(mEnv.MaxCell,0.0) ; // Cell energies          of      product molecule.
    vector<double> rctsCellDOS(mEnv.MaxCell,0.0) ; // Convoluted cell density of states of reactants.
    p_pdt1->getCellDensityOfStates(pdt1CellDOS) ;
    p_rcts->getCellDensityOfStates(rctsCellDOS) ;
    p_pdt1->getCellEnergies       (pdt1CellEne) ;

    // Allocate space to hold Micro-canonical rate coefficients.
    cellKfmc.resize(mEnv.MaxCell); // no need to initialize

    int activ_ene = 0;
    if(IsNan(_einf))
    {
      stringstream errorMsg;
      errorMsg << "To use MesmerILT for reaction " << pReact->getName()
               << " the Activation Energy needs to be set.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
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

    DOSconvolution(work, pdt1CellDOS, conv);

    for (int i = 0; i < activ_ene; ++i)  cellKfmc[i] = 0.;
    for (int i = 0; i < (mEnv.MaxCell - activ_ene); ++i){
      cellKfmc[i + activ_ene] = _ant * conv[i] / rctsCellDOS[i + activ_ene];
      //ctest << "cellKfmc[" << i + activ_ene << "] = " << cellKfmc[i + activ_ene] << endl;
    }

    {stringstream errorMsg;
    errorMsg << "MesmerILT calculation completed";
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

    return true;
  }

}//namespace
