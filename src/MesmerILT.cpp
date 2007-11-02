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
    ModelledMolecule * p_rct1 = static_cast<ModelledMolecule *>(bi_molecularspecies[0]);
    ModelledMolecule * p_rct2 = static_cast<ModelledMolecule *>(bi_molecularspecies[1]); //not sure about this _2007_10_18__12_48_34_ CHL

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
    vector<double> rctsCellDOS(mEnv.MaxCell,0.0) ; // Convoluted cell density of states of reactants.
    p_pdt1->cellDensityOfStates(pdt1CellDOS, mEnv) ;
    p_pdt1->cellEnergies       (pdt1CellEne, mEnv) ;

    // Convolution of reactant vibrational DOS onto rotational DOS
    countDimerCellDOS(p_rct1, p_rct2, rctsCellDOS, mEnv);

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

    DOSconvolution(work, pdt1CellDOS, conv);

    for (int i = 0; i < activ_ene; ++i)  cellKfmc[i] = 0.;
    for (int i = 0; i < (mEnv.MaxCell - activ_ene); ++i){
      cellKfmc[i + activ_ene] = _ant * conv[i] / rctsCellDOS[i + activ_ene];
      cout << "cellKfmc[" << i + activ_ene << " = " << cellKfmc[i + activ_ene] << endl;
    }

    return true;
  }

  // provide a function to define particular counts of the convoluted DOS of two molecules
  bool MesmerILT::countDimerCellDOS(ModelledMolecule* p_mol1, ModelledMolecule* p_mol2, vector<double> &dimerCellDOS, const MesmerEnv &mEnv){
    //----------
    // Differetiating the rotors
    // <three types of rotors: (0) non-rotor (1) 2-D linear, (2) 3-D symmetric top, (3) 3-D asymmetric top>
    // hence there are 9 combinations for convolution:
    // 0-1, 0-2, 0-3, 1-1, 1-2, 1-3, 2-2, 2-3, 3-3
    // where the first three convolutions are trivial
    //

    vector<double> rotCs1; int rotor1Type = p_mol1->get_rotConsts(rotCs1);
    vector<double> rotCs2; int rotor2Type = p_mol2->get_rotConsts(rotCs2);
    vector<double> mol1CellEne(mEnv.MaxCell,0.0);
    p_mol1->cellEnergies(mol1CellEne, mEnv) ;

    double I1 = 0., I2 = 0.,                                            //for 2-D linear rotors
           I1X = 0., I2X = 0., I1Y = 0., I2Y = 0., I1Z = 0., I2Z = 0.,  //for 3-D symmetric/asymmetric/spherical top rotors
           s1 = p_mol1->get_Sym(), s2 = p_mol2->get_Sym(),
           q1 = static_cast<double>(p_mol1->get_SpinMultiplicity()),
           q2 = static_cast<double>(p_mol2->get_SpinMultiplicity());

    /*All rotational constants are sorted --- the largest one at the head, smallest one at the tail*/


    //---------------------- Assign rotational constants ----------------------
    if      (rotor1Type == -4){}
    else if (rotor1Type ==  0){ I1 = rotCs1[0]; }// 2-D linear
    else                      { I1X = rotCs1[0]; I1Y = rotCs1[1]; I1Z = rotCs1[2];}// 3-D symmetric/asymmetric/spherical top
    if      (rotor2Type == -4){}
    else if (rotor2Type ==  0){ I2 = rotCs2[0]; }// 2-D linear
    else                      { I2X = rotCs2[0]; I2Y = rotCs1[1]; I2Z = rotCs2[2];}// 3-D symmetric/asymmetric/spherical top

    //------------------------------------------------------------------------------
    //Density of states from ILT of the product of partition functions of two rotors
    int rotorType = rotor1Type + rotor2Type;
    if      (rotorType == -8) return false;              // both not rotors
    double cnt = q1 * q2 / (s1 * s2);

    if      (rotorType <  0 && rotor1Type < rotor2Type){ // only p_mol2 a rotor
      p_mol2->cellDensityOfStates(dimerCellDOS, mEnv);
      return true;
    }
    else if (rotorType <  0 && rotor1Type > rotor2Type){ // only p_mol1 a rotor
      p_mol1->cellDensityOfStates(dimerCellDOS, mEnv);
      return true;
    }
    else if (rotorType == 0){                            // both 2-D linear
      cnt /= (I1 * I2);
      for (int i = 0; i < mEnv.MaxCell; ++i) dimerCellDOS[i] = cnt * mol1CellEne[i];
    }
    else if (rotorType == 2 && rotor1Type < rotor2Type){ // 2-D linear vs 3-D symmetric/asymmetric/spherical top
      cnt *= (4./(3.* I1 * sqrt(I2X * I2Y * I2Z)));
      for (int i = 0; i < mEnv.MaxCell; ++i) dimerCellDOS[i] = cnt * pow(mol1CellEne[i], 1.5);
    }
    else if (rotorType == 2 && rotor1Type > rotor2Type){ // 3-D symmetric/asymmetric/spherical top vs 2-D linear
      cnt *= (4./(3.* I2 * sqrt(I1X * I1Y * I1Z)));
      for (int i = 0; i < mEnv.MaxCell; ++i) dimerCellDOS[i] = cnt * pow(mol1CellEne[i], 1.5);
    }
    else if (rotorType == 4){                            // both 3-D symmetric/asymmetric/spherical top
      cnt *= (M_PI /(2. * sqrt(I1X * I1Y * I1Z * I2X * I2Y * I2Z)));
      for (int i = 0; i < mEnv.MaxCell; ++i) dimerCellDOS[i] = cnt * (mol1CellEne[i] * mol1CellEne[i]);
    }


    //-----------------------------------------------------------------------------------------------------
    // convolution of vibrational DOS onto rotational DOS -- loop through all frequencies of both molecules
    vector<double> vfMol1; p_mol1->get_VibFreq(vfMol1);
    vector<double> vfMol2; p_mol2->get_VibFreq(vfMol2);
    for (int i = 0; i < vfMol1.size(); ++i){
      int nFreq = static_cast<int>(vfMol1[i] + 0.5);
      for (int j = 0; j < (mEnv.MaxCell - nFreq); ++j){
        dimerCellDOS[nFreq + j] += dimerCellDOS[j];
      }
    }
    for (int i = 0; i < vfMol2.size(); ++i){
      int nFreq = static_cast<int>(vfMol2[i] + 0.5);
      for (int j = 0; j < (mEnv.MaxCell - nFreq); ++j){
        dimerCellDOS[nFreq + j] += dimerCellDOS[j];
      }
    }

    return true;
  }

}//namespace
