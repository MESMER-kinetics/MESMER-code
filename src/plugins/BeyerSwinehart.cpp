#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

using namespace std;
namespace mesmer
{
  class BeyerSwinehart : public DensityOfStatesCalculator
  {
  public:

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell);

    // Provide a function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class calculates a complete DOS: it is not an extra class. 
    BeyerSwinehart(const std::string& id) : DensityOfStatesCalculator(id, false){}

    virtual ~BeyerSwinehart() {}
    virtual BeyerSwinehart* Clone() { return new BeyerSwinehart(*this); }

  } ;

  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  BeyerSwinehart theBeyerSwinehart("BeyerSwinehart");
  //************************************************************


  // Provide a function to define particular counts of the DOS of a molecule.
  bool BeyerSwinehart::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {
    vector<double> VibFreq ;
    pDOS->get_VibFreq(VibFreq) ;

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    // Implementation of the Beyer-Swinehart algorithm.
    for (size_t j(0) ; j < VibFreq.size() ; ++j ) {
      size_t freq = static_cast<size_t>(VibFreq[j]) ;
      for (size_t i(0) ; i < MaximumCell - freq ; ++i ){
        cellDOS[i + freq] += cellDOS[i] ;
      }
    }
	pDOS->setCellDensityOfStates(cellDOS) ;

    return true;
  }

  // Calculate contribution to canonical partition function.
  double BeyerSwinehart::canPrtnFnCntrb(gDensityOfStates* gdos, double beta) {

    double qtot(1.0) ; 
    vector<double> vibFreq; 
    gdos->get_VibFreq(vibFreq);
    for (size_t j(0) ; j < vibFreq.size() ; ++j ) {
      qtot /= (1.0 - exp(-beta*vibFreq[j])) ;
    }

	return qtot ;
  }

}//namespace
