#include "SimpleILT.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  SimpleILT theSimpleILT("SimpleILT");
  SimpleILT oldSimpleILT("Simple ILT");
  //************************************************************

  //
  // This method calculates the reaction flux by Laplace inversion
  // of the Arrhenius equation for a reaction proceeding in the 
  // unimolecular direction.
  //

  bool SimpleILT::calculateMicroRateCoeffs(Reaction* pReact)
  {
    vector<Molecule *> Isomers ;
    int nIsomers = pReact->get_unimolecularspecies(Isomers) ;  
    Molecule *p_rcts = Isomers[0] ;

    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    // Allocate some work space for and obtain density of states of the unimolecuar reactant.

    vector<double> rctsCellDOS; 
    if(!p_rcts->getDOS().getCellDensityOfStates(rctsCellDOS))
      return false;

    // Obtain the Arrhenius parameters.

    const int nEinf = int(pReact->get_ThresholdEnergy()) ;  //<-- SHR: This should be activation energy NOT threshold energy!!
    const double preExp = m_PreExp ;

    // Calculate microcanonical rate coefficients using simple ILT expression.

    for (int i = nEinf; i < MaximumCell ; ++i ) {
      rxnFlux[i] = preExp * rctsCellDOS[i-nEinf];
    }

    // the flux bottom energy is equal to the well bottom of the source term
    pReact->setCellFluxBottom(p_rcts->getDOS().get_zpe() - pReact->getEnv().EMin);

    return true;
  }

  //
  // This function retrieves the activation/threshold energy for an association reaction.
  //
  double SimpleILT::get_ThresholdEnergy(Reaction* pReac) {

    //if (m_EInf < 0.0) now checked during parsing
    //  cerr << "Providing negative E_infinity in Reaction " << getName() << " is invalid.";

    double RxnHeat = pReac->getHeatOfReaction(); 
    double threshold(0.0) ;

    if (m_isRvsILTpara){
      threshold = (m_EInf > 0.0) ? m_EInf + RxnHeat : RxnHeat ;
    } else {
      threshold = m_EInf ;
    }

    if (threshold < RxnHeat){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }

    return threshold ;

  }
}//namespace
