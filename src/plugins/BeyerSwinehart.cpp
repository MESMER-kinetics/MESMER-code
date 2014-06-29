#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

using namespace std;
namespace mesmer
{
  class BeyerSwinehart : public DensityOfStatesCalculator
  {
  public:

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env);

    // Function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Provide a function to calculate the zero point energy of a molecule.
	virtual double ZeroPointEnergy(gDensityOfStates* gdos) ;

    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class calculates a complete DOS: it is not an extra class. 
    BeyerSwinehart(const char* id) : m_id(id) { Register(); }

    virtual ~BeyerSwinehart() {}
    virtual const char* getID()  { return m_id; }
    virtual BeyerSwinehart* Clone() { return new BeyerSwinehart(*this); }

  private:
    const char* m_id;
  } ;

  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  BeyerSwinehart theBeyerSwinehart("BeyerSwinehart");
  //************************************************************


  // Provide a function to define particular counts of the DOS of a molecule.
  bool BeyerSwinehart::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
	const size_t MaximumCell = env.MaxCell ;

    vector<double> VibFreq ;
    pDOS->get_VibFreq(VibFreq) ;

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    // Implementation of the Beyer-Swinehart algorithm.
    for (size_t j(0) ; j < VibFreq.size() ; ++j ) {
      size_t freq = static_cast<size_t>(nint(VibFreq[j])) ;
	  if (freq > MaximumCell) {
		// This is to catch those occassional cases where the first excited 
		// vibrational state is above the cutoff, which can occur at low 
		// temperatures. 
		continue ;
	  }
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

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int BeyerSwinehart::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> vibFreq; 
    gdos->get_VibFreq(vibFreq);

    return vibFreq.size() ;
  }

  // Provide a function to calculate the zero point energy of a molecule.
  double BeyerSwinehart::ZeroPointEnergy(gDensityOfStates* gdos) {

	vector<double> vibFreq; 
    gdos->get_VibFreq(vibFreq);

	double ZPE(0.0) ;
    for (size_t j(0) ; j < vibFreq.size() ; ++j ) {
      ZPE += vibFreq[j] ;
    }
    return ZPE*0.5 ;
  }

}//namespace
