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
    if (pReact->isUnimolecular()){                  // if it's unimolecular 
      if(!calculateUnimolecularMicroRates(pReact))  // and the microrate calculation is unsuccessful return false
        return false;
    }
    else if(!pReact->isUnimolecular()){            // if it's not unimolecular
      if(!calculateAssociationMicroRates(pReact))  // and the microrate calculation is unsuccessful return false
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
    const int    Einf   = (int)(pReact->get_EInf());

    Molecule*  p_rct = pReact->get_reactant();
    int MaximumCell = pReact->getEnv().MaxCell;
    vector<Molecule*> v_pdts;
    const int numberOfProducts = pReact->get_products(v_pdts);

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctCellEne;  // Cell energies of reactant molecule.
    vector<double> rctCellDOS; //  Cell density of states of reactant.

    getCellEnergies(MaximumCell, rctCellEne);
    p_rct->getDOS().getCellDensityOfStates(rctCellDOS);

    // Allocate space to hold microcanonical rate coefficients for dissociation.
    // This line below pass the reference pointer of m_CellFlux to the vector by (&), so what the code does on
    // rxnFlux will in fact work on m_CellFlux.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    const double gammaValue = MesmerGamma(Ninf);
    const double beta0      = 1.0/(boltzmann_RCpK*Tinf);

    // If the activation energy specified is for the reverse direction.
    double Keq(1.0);
    if (pReact->isReverseReactionILT_Ea()){
      const double Q_R        = p_rct->getDOS().rovibronicGrnCanPrtnFn();
      const double Q_p        = pReact->pdtsRovibronicGrnCanPrtnFn();
      Keq *= Q_p / Q_R;
      // Heat of reaction contribution is included in activation energy
      if (numberOfProducts == 2){
        Keq *= translationalContribution(v_pdts[0]->getMass(), v_pdts[1]->getMass(), pReact->getEnv().beta);
      }
    }
    const double constant   = Ainf * Keq * pow(beta0,Ninf)/gammaValue;

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
      rxnFlux[i] = constant * conv[i];

    if (pReact->isReverseReactionILT_Ea()){
      // exp(-(E0 + Ea)/kB T)
      pReact->setCellFluxBottom(pReact->get_relative_pdtZPE() + Einf);
    }
    else{
      pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + Einf);
    }

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

    vector<Molecule *> unimolecularspecies;
    pAssocReaction->get_unimolecularspecies(unimolecularspecies);

    Molecule*  p_pdt1 = unimolecularspecies[0];
    Molecule*  p_rct1 = pAssocReaction->get_pseudoIsomer();
    Molecule*  p_rct2 = pAssocReaction->get_excessReactant();

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
    // This line below pass the reference pointer of m_CellFlux to the vector by (&), so what the code does on
    // rxnFlux will in fact work on m_CellFlux.
    vector<double>& rxnFlux = pAssocReaction->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

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
      rxnFlux[i] = _ant * conv[i];

    // the flux bottom energy is equal to the well bottom of the reactant
    pAssocReaction->setCellFluxBottom(pReact->get_relative_rctZPE() + Einf);

    cinfo << "Association ILT calculation completed" << endl;
    return true;
  }

}//namespace
