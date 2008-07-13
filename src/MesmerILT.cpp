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
		//
		// Check input to see if we have the correct type of reaction.
		//
        AssociationReaction *pAssocReaction = dynamic_cast<AssociationReaction*>(pReact) ;
        if (!pAssocReaction) {
            cerr << "The Mesmer ILT method is valid for association reactions only." ;
            return false ;
        }

        //-----------------
        //starting variables block
        const double Ninf   = pAssocReaction->get_NInf(); // constraint: Ninf > -1.5
        const double Tinf   = 1. / (boltzmann_RCpK * pAssocReaction->getEnv().beta);
        const double Ainf   = pAssocReaction->get_PreExp();
        const int    Einf   = int(pAssocReaction->get_ThresholdEnergy());
        // double tp_C = 3.24331e+20; // Defined in Constant.h, constant used in the translational partition function
        //-----------------

        SuperMolecule* p_rcts = pAssocReaction->get_bi_molecularspecies();
        if (!p_rcts){
            cerr << "Not a valid bi-molecularspecies";
            return false;
        }

        vector<ModelledMolecule *> unimolecularspecies;
        pAssocReaction->get_unimolecularspecies(unimolecularspecies);

        ModelledMolecule*  p_pdt1 = unimolecularspecies[0];
        ModelledMolecule*  p_rct1 = pAssocReaction->get_pseudoIsomer();
        ModelledMolecule*  p_rct2 = pAssocReaction->get_excessReactant();

        // Get molecular specific values
        //    const double edg_a = static_cast<double>(p_rct1->getSpinMultiplicity());
        //    const double edg_b = static_cast<double>(p_rct2->getSpinMultiplicity());
        //    const double edg_c = static_cast<double>(p_pdt1->getSpinMultiplicity());
        const double ma = p_rct1->getMass();
        const double mb = p_rct2->getMass();
        const double mc = p_pdt1->getMass();
        int MaximumCell = pAssocReaction->getEnv().MaxCell;

        // Allocate some work space for density of states and extract densities of states from molecules.
        vector<double> rctsCellEne; // Cell energies          of      product molecule.
        vector<double> rctsCellDOS; // Convoluted cell density of states of reactants.

        p_rcts->getCellEnergies       (rctsCellEne) ;
        p_rcts->getCellDensityOfStates(rctsCellDOS) ;

        // Allocate space to hold microcanonical rate coefficients for dissociation.
        vector<double>& TSFlux = pAssocReaction->get_CellFlux();
        TSFlux.clear();
        TSFlux.resize(MaximumCell, 0.0); 

        const double gammaValue = MesmerGamma(Ninf + 1.5);

        //double _ant = Ainf * tp_C * (edg_a * edg_b / edg_c) * pow( ( ma * mb / mc), 1.5 ) / gammaValue;
        // Replace the above line because electronic degeneracies were already accounted for in DOS calculations
        double _ant = Ainf * tp_C * pow( ( ma * mb / mc), 1.5 ) / gammaValue;
        _ant /= (pow((Tinf * boltzmann_RCpK), Ninf));

        vector<double> work(MaximumCell);;
        vector<double> conv(MaximumCell);

        double pwr = Ninf + .5;
        for (int i = 0; i < MaximumCell; ++i) {
            work[i] = pow(rctsCellEne[i], pwr);
        }

        FastLaplaceConvolution(work, rctsCellDOS, conv);    // FFT convolution replaces the standard convolution
        //    Convolution(work, rctsCellDOS, conv);  // standard convolution

        if (Einf >= 0){
            for (int i = Einf; i < MaximumCell; ++i){
                TSFlux[i] = _ant * conv[i - Einf];
            }
        }
        else{
            for (int i = 0; i < MaximumCell + Einf; ++i){
                TSFlux[i] = _ant * conv[i - Einf];
            }
        }

        // the flux bottom energy is equal to the well bottom of the product
        pAssocReaction->setCellFluxBottom(p_rcts->get_relative_ZPE());

        cinfo << "ILT calculation completed" << endl;

        return true;
    }

}//namespace
