#include "MesmerILT.h"
#include "AssociationReaction.h"
#include "Rdouble.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  MesmerILT theMesmerILT("MesmerILT");
  MesmerILT oldMesmerILT("Mesmer ILT");
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
    //
    // Need to determine if the supplied Arrhenius parameters are for the association 
    // or dissociation direction and then invoke the appropriate algorithm.
    //
    double relative_ZPE(0.0) ;
    if (pReact->isReverseReactionILT_Ea()){

      vector<Molecule *> products ; 
      int numberOfProducts = pReact->get_products(products) ;

      if (numberOfProducts != 2) 
        return false ;

      Molecule* p_rct1 = pReact->get_reactant();
      Molecule* p_pdt1 = products[0];
      Molecule* p_pdt2 = products[1];

      const double ma = p_pdt1->getStruc().getMass();
      const double mb = p_pdt2->getStruc().getMass();
      const double mc = p_rct1->getStruc().getMass();

      // Allocate some work space for density of states and extract densities of states from molecules.
      vector<double> pdtsCellDOS; // Convoluted cell density of states of reactants.

      if(!countDimerCellDOS(p_pdt1->getDOS(), p_pdt2->getDOS(), pdtsCellDOS))
        return false;

      BimolecularConvolution(pReact, pdtsCellDOS, ma, mb, mc) ;

      relative_ZPE = pReact->get_relative_pdtZPE();

    } else {

      UnimolecularConvolution(pReact) ;

      relative_ZPE = pReact->get_relative_rctZPE() ;

    }
    pReact->setCellFluxBottom(relative_ZPE + pReact->get_EInf());

    cinfo << "Unimolecular ILT calculation completed" << endl;
    return true;
  }

  bool MesmerILT::calculateAssociationMicroRates(Reaction* pReact)
  {
    AssociationReaction *pAssocReaction = dynamic_cast<AssociationReaction*>(pReact) ;
    if(!pAssocReaction){
      cerr << "The MesmerILT method is not available for Irreversible Exchange Reactions"<< endl;
      return false ;
    }

    vector<Molecule *> unimolecularspecies;
    pAssocReaction->get_unimolecularspecies(unimolecularspecies);

    Molecule*  p_pdt1 = unimolecularspecies[0];
    Molecule*  p_rct1 = pAssocReaction->get_pseudoIsomer();
    Molecule*  p_rct2 = pAssocReaction->get_excessReactant();

    const double ma = p_rct1->getStruc().getMass();
    const double mb = p_rct2->getStruc().getMass();
    const double mc = p_pdt1->getStruc().getMass();

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctsCellDOS; // Convoluted cell density of states of reactants.

    pAssocReaction->getRctsCellDensityOfStates(rctsCellDOS) ;

    // Perform convolution.

    BimolecularConvolution(pReact, rctsCellDOS, ma, mb, mc) ;

    // the flux bottom energy is equal to the well bottom of the reactant
    pAssocReaction->setCellFluxBottom(pReact->get_relative_rctZPE() + pReact->get_EInf());

    cinfo << "Association ILT calculation completed" << endl;

    return true;
  }

  bool MesmerILT::UnimolecularConvolution(Reaction* pReact)
  {
    //
    // Obtain Arrhenius parameters. Note constraint: Ninf >= 0.0
    //
    const double Ninf = pReact->get_NInf(); 
    const double Tinf = pReact->get_TInf();
    const double Ainf = pReact->get_PreExp();

    Molecule* p_rct = pReact->get_reactant();
    int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate some work space for density of states and extract densities of states from reactant.
    vector<double> rctCellDOS; 
    if(!p_rct->getDOS().getCellDensityOfStates(rctCellDOS))
      return false;

    //
    // Initialize reaction flux vector.
    //
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    const double gammaValue = MesmerGamma(Ninf);
    const double beta0      = 1.0/(boltzmann_RCpK*Tinf);
    const double constant   = Ainf * pow(beta0,Ninf)/gammaValue;

    vector<double> conv;
    if (fabs(Ninf) > 0.0 ) {
      //
      // The expression held in the elements of the vector work has been altered from the
      // simple mean to analytic integral_x_to_y{E^(Ninf-1)dE}, where x and y are lower and
      // upper energy limits of the cell respectively.
      //
      vector<double> work(MaximumCell, 0.0);

      for (int i = 0; i < MaximumCell; ++i){
        work[i] = (pow(i+1,Ninf)-pow(i,Ninf))/Ninf;
      }
      FastLaplaceConvolution(work, rctCellDOS, conv);    // FFT convolution replaces the standard convolution
    } else {
      // Need to account for the case when Ninf is zero.
      for (int i = 0; i < MaximumCell; ++i){
        conv[i] = rctCellDOS[i] ;
      }
    }

    for (int i = 0; i < MaximumCell; ++i)
      rxnFlux[i] = constant * conv[i];

    return true;
  }

  bool MesmerILT::BimolecularConvolution(Reaction* pReact, vector<double>& ConvolvedCellDOS, double ma, double mb, double mc)
  {
    //
    // Obtain Arrhenius parameters. Note constraint: Ninf > -1.5
    //
    const double Ninf = pReact->get_NInf(); 
    const double Tinf = pReact->get_TInf();
    const double Ainf = pReact->get_PreExp();

    //
    // Initialize reaction flux vector.
    //
    int MaximumCell = pReact->getEnv().MaxCell;
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    const double gammaValue = MesmerGamma(Ninf + 1.5);

    // Note electronic degeneracies were already accounted for in DOS calculations.
    // tp_C = 3.24331e+20: defined in Constant.h, constant used in the translational
    // partition function.

    double _ant = Ainf * tp_C * pow( ( ma * mb / mc), 1.5 ) / gammaValue;
    _ant /= (pow((Tinf * boltzmann_RCpK), Ninf));

    //
    // The expression held in the elements of the vector work has been altered from the
    // simple power of the mean value to analytic integral_x_to_y{E^(Ninf-1)dE}, where
    // x and y are lower and upper energy limits of the cell respectively.
    //
    vector<double> work(MaximumCell);
    for (int i = 0; i < MaximumCell; ++i){
      work[i] = (pow(i+1,Ninf+1.5)-pow(i,Ninf+1.5))/(Ninf+1.5);
    }

    vector<double> conv;
    FastLaplaceConvolution(work, ConvolvedCellDOS, conv);    // FFT convolution replaces the standard convolution
    //    Convolution(work, rctsCellDOS, conv);  // standard convolution

    for (int i = 0; i < MaximumCell; ++i)
      rxnFlux[i] = _ant * conv[i];

    return true;
  }


}//namespace
