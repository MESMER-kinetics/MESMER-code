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
  //   Einf:           activation energy (experimental value)
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
    // Check to see what type of reaction we have
    if (pReact->isUnimolecular()){      // if it's unimolecular and there microrate calculation is unsuccessful
      if(!calculateUnimolecularMicroRates(pReact))  // return false
        return false;
    }
    else if(!pReact->isUnimolecular()){ // if it's not unimolecular and the microrate calculation is unsuccesful
      if(!calculateAssociationMicroRates(pReact))  // return false
        return false;
    }
    return true;
  }

  bool MesmerILT::calculateUnimolecularMicroRates(Reaction* pReact)
  {
    //starting variables block
    const double Ninf   = pReact->get_NInf(); // constraint: Ninf != 0
    const double Tinf   = pReact->get_TInf();
    const double Ainf   = pReact->get_PreExp();
    const int Einf   = (int)(pReact->get_EInf());

    ModelledMolecule*  p_rct = pReact->get_reactant();
    int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctCellEne;  // Cell energies of reactant molecule.
    vector<double> rctCellDOS; //  Cell density of states of reactant.

    getCellEnergies(MaximumCell, rctCellEne);
    p_rct->getCellDensityOfStates(rctCellDOS);

    // Allocate space to hold microcanonical rate coefficients for dissociation.
    // This line below pass the reference pointer of m_CellTSFlux to the vector by (&), so what the code does on
    // TSFlux will in fact work on m_CellTSFlux.
    vector<double>& TSFlux = pReact->get_CellFlux();
    TSFlux.clear();
    TSFlux.resize(MaximumCell, 0.0); 

    const double gammaValue = MesmerGamma(Ninf);
    const double beta0      = 1.0/(boltzmann_RCpK*Tinf);
    const double constant   = Ainf * pow(beta0,Ninf)/gammaValue;
    const double pwr        = Ninf - 1.0;
    
    /*The expression of the work vector is changed so that every member is from analytic solution of
    integral_x_to_y{E^(Ninf-1)dE}, where x and y are lower and upper energy limits of the cell respectively.*/
    vector<double> work(MaximumCell);
    for (int i = 0; i < MaximumCell; ++i){
      //work[i] = pow(rctCellEne[i], pwr);
      work[i] = (pow(i+1,Ninf)-pow(i,Ninf))/Ninf;
    }
    vector<double> conv;
    FastLaplaceConvolution(work, rctCellDOS, conv);    // FFT convolution replaces the standard convolution

    for (int i = 0; i < MaximumCell; ++i)
      TSFlux[i] = constant * conv[i];

    pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + Einf);

    cinfo << "Unimolecular ILT calculation completed" << endl;
    return true;
  }

  bool MesmerILT::calculateAssociationMicroRates(Reaction* pReact)
  {
    AssociationReaction *pAssocReaction = dynamic_cast<AssociationReaction*>(pReact) ;
    if(!pAssocReaction){
      cerr << "The Mesmer ILT method is not available for Irreversible Exchange Reactions"<< endl;
      return false ;
    }

    //starting variables block
    const double Ninf   = pAssocReaction->get_NInf(); // constraint: Ninf > -1.5
    const double Tinf   = pAssocReaction->get_TInf();
    const double Ainf   = pAssocReaction->get_PreExp();
    const int Einf   = (int)(pAssocReaction->get_EInf());
    // double tp_C = 3.24331e+20; // Defined in Constant.h, constant used in the translational partition function
    //-----------------

    vector<ModelledMolecule *> unimolecularspecies;
    pAssocReaction->get_unimolecularspecies(unimolecularspecies);

    ModelledMolecule*  p_pdt1 = unimolecularspecies[0];
    ModelledMolecule*  p_rct1 = pAssocReaction->get_pseudoIsomer();
    ModelledMolecule*  p_rct2 = pAssocReaction->get_excessReactant();

    const double ma = p_rct1->getMass();
    const double mb = p_rct2->getMass();
    const double mc = p_pdt1->getMass();
    int MaximumCell = pAssocReaction->getEnv().MaxCell;

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctsCellEne; // Cell energies          of      product molecule.
    vector<double> rctsCellDOS; // Convoluted cell density of states of reactants.

    getCellEnergies(MaximumCell, rctsCellEne) ;
    pAssocReaction->getRctsCellDensityOfStates(rctsCellDOS) ;

    // Allocate space to hold microcanonical rate coefficients for dissociation.
    // This line below pass the reference pointer of m_CellTSFlux to the vector by (&), so what the code does on
    // TSFlux will in fact work on m_CellTSFlux.
    vector<double>& TSFlux = pAssocReaction->get_CellFlux();
    TSFlux.clear();
    TSFlux.resize(MaximumCell, 0.0); 

    const double gammaValue = MesmerGamma(Ninf + 1.5);

    //double _ant = Ainf * tp_C * (edg_a * edg_b / edg_c) * pow( ( ma * mb / mc), 1.5 ) / gammaValue;
    // Replace the above line because electronic degeneracies were already accounted for in DOS calculations
    double _ant = Ainf * tp_C * pow( ( ma * mb / mc), 1.5 ) / gammaValue;
    _ant /= (pow((Tinf * boltzmann_RCpK), Ninf));

    //double pwr = Ninf + .5;
    /*The expression of the work vector is changed so that every member is from analytic solution of
    integral_x_to_y{E^(1/2+Ninf)dE}, where x and y are lower and upper energy limits of the cell respectively.*/
    vector<double> work(MaximumCell);
    for (int i = 0; i < MaximumCell; ++i){
      //work[i] = pow(rctsCellEne[i], pwr);
      work[i] = (pow(i+1,Ninf+1.5)-pow(i,Ninf+1.5))/(Ninf+1.5);
    }

    vector<double> conv;
    FastLaplaceConvolution(work, rctsCellDOS, conv);    // FFT convolution replaces the standard convolution
    //    Convolution(work, rctsCellDOS, conv);  // standard convolution

    for (int i = 0; i < MaximumCell; ++i)
      TSFlux[i] = _ant * conv[i];

    // the flux bottom energy is equal to the well bottom of the reactant
    pAssocReaction->setCellFluxBottom(pReact->get_relative_rctZPE() + Einf);

    cinfo << "Association ILT calculation completed" << endl;
    return true;
  }

}//namespace
