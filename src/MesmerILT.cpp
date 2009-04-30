#include "MesmerILT.h"
#include "AssociationReaction.h"
#include "Rdouble.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance)
  MesmerILT theMesmerILT("MesmerILT");
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

    // Read ILT parameters
  bool MicroRateCalculator::ReadParameters(Reaction* pReact) {
    PersistPtr ppReac = pReact->get_PersistentPointer();

    // OpenBabel outputs <rateParameters> <A> <n> <E>
    // Attempt to read these first; 
    // if not present read the mesmer version which will add the default if necessary.
    //**TODO preExponential
    PersistPtr ppActEne, ppPreExponential;
    const char* pActEnetxt=NULL, *pPreExptxt=NULL;
    PersistPtr ppRateParams = ppReac->XmlMoveTo("rateParameters") ;

    if(ppRateParams) {
      ppActEne = ppRateParams->XmlMoveTo("E") ;
      pActEnetxt = ppRateParams->XmlReadValue("E", optional);
      ppPreExponential = ppRateParams->XmlMoveTo("A") ;
      pPreExptxt = ppRateParams->XmlReadValue("A");
    }
    else {
      ppActEne = ppReac->XmlMoveTo("me:activationEnergy") ;
      pActEnetxt = ppReac->XmlReadValue("me:activationEnergy");
      ppPreExponential = ppReac->XmlMoveTo("me:preExponential") ;
      pPreExptxt = ppReac->XmlReadValue("me:preExponential");
    }

    if (pActEnetxt)
    {
      pReact->set_revILT(ppActEne->XmlReadBoolean("reverse")); //specify the direction of the following ILT parameters
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput;
      if (unitsTxt){
        unitsInput = unitsTxt;
      }
      else{
        unitsInput = "kJ/mol";
      }
      const char* pLowertxt = ppActEne->XmlReadValue("lower", optional);
      const char* pUppertxt = ppActEne->XmlReadValue("upper", optional);
      const char* pStepStxt = ppActEne->XmlReadValue("stepsize", optional);
      double value(getConvertedEnergy(unitsInput, tmpvalue));
      if (pLowertxt && pUppertxt){
        double tmpvalueL(0.0), tmpvalueU(0.0), tmpstepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt);
        s3 >> tmpvalueL; s4 >> tmpvalueU; s5 >> tmpstepsize;
        double valueL(getConvertedEnergy(unitsInput, tmpvalueL));
        double valueU(getConvertedEnergy(unitsInput, tmpvalueL));
        double stepsize(getConvertedEnergy(unitsInput, tmpstepsize));
        //set_EInf(valueL, valueU, stepsize);
        //SETRANGE(pReact->set_EInf,valueL,valueU,stepsize);
        pReact->set_EInf(NaN);
        ActiveRdoubles.back()->set_range(valueL,valueU,stepsize);
        cinfo << " Set range of EInf" << endl; 
      }
      else{
        pReact->set_EInf(value);
      }
    }
    else{
      cerr << "Specifying ILT without activation energy provided in reaction "
           << this->getName() << ". Please correct input file.";
      return false;
    }

    if (pPreExptxt)
    {
      double value = 0.0;
      stringstream s2(pPreExptxt); s2 >> value ;
      const char* pLowertxt = ppPreExponential->XmlReadValue("lower", optional);
      const char* pUppertxt = ppPreExponential->XmlReadValue("upper", optional);
      const char* pStepStxt = ppPreExponential->XmlReadValue("stepsize", optional);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt);
        s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        //set_PreExp(valueL, valueU, stepsize);
        //SETRANGE(pReact->set_PreExp,valueL,valueU,stepsize);
        pReact->set_PreExp(NaN);
        ActiveRdoubles.back()->set_range(valueL,valueU,stepsize);
        cinfo << " Set range of PreExp" << endl; 

      }
      else{
        pReact->set_PreExp(value);
      }
    }
    else{
      cerr << "Specifying ILT without pre-exponential term provided in reaction " << this->getName() << ". Please correct input file.";
      return false;
    }

    const char* pNInftxt = ppReac->XmlReadValue("me:nInfinity", optional);
    if (pNInftxt)
    {
      PersistPtr ppNInf = ppReac->XmlMoveTo("me:nInfinity") ;
      double value = 0.0;
      stringstream s2(pNInftxt); s2 >> value ;
      const char* pLowertxt = ppNInf->XmlReadValue("lower", false);
      const char* pUppertxt = ppNInf->XmlReadValue("upper", false);
      const char* pStepStxt = ppNInf->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        //set_NInf(valueL, valueU, stepsize);
        //SETRANGE(pReact->set_NInf,valueL,valueU,stepsize);

        pReact->set_EInf(NaN);
        ActiveRdoubles.back()->set_range(valueL,valueU,stepsize);
        cinfo << " Set range of NInf" << endl; 
      }
      else{
        pReact->set_NInf(value);
      }
    }

    double TInf = ppReac->XmlReadDouble("me:TInfinity");
    if(TInf <= 0) {
      cinfo << "Tinfinity is less than or equal to 0; set to the default value of 298 K" << endl;
      TInf = 298;
    }
    else
      pReact->set_TInf(TInf);         // else set Tinf to the value in the input
    
    return true;
  }

}//namespace
