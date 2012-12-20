#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

using namespace std;
namespace mesmer
{
  class ClassicalRotor : public DensityOfStatesCalculator
  {
  public:

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell);

    // Function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in TopPlugin
    // This class calculates a complete DOS: it is not an extra class. 
    ClassicalRotor(const char* id) : m_id(id){ Register(); }

    virtual ~ClassicalRotor() {}
    virtual const char* getID()  { return m_id; }

    virtual ClassicalRotor* Clone() { return new ClassicalRotor(*this); }
  private:
    const char* m_id;
  } ;

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  ClassicalRotor theClassicalRotor("ClassicalRotors");
  //************************************************************


  // Provide a function to define particular counts of the DOS of a molecule.
  bool ClassicalRotor::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {
    vector<double> VibFreq ;
    pDOS->get_VibFreq(VibFreq) ;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0) ;

    //
    // Initialize density of states array using calculated rotational
    // density of state from inverse Laplace transform of rotors.
    //
    vector<double> rotConst;
    RotationalTop rotorType = pDOS->get_rotConsts(rotConst);
    double sym = pDOS->get_Sym();
    double qele = pDOS->getSpinMultiplicity();
    double cnt = 0.;

    switch (rotorType){
      case NONLINEAR: //3-D symmetric/asymmetric/spherical top
        cnt = qele * sqrt(4./(rotConst[0] * rotConst[1] * rotConst[2]))/sym ;
        for (size_t i(0) ; i < MaximumCell ; ++i )
          cellDOS[i] = cnt*sqrt(cellEne[i]) ;
        break;
      case LINEAR: //2-D linear
        cnt = qele / (rotConst[0] * sym);
        for (size_t i(0) ; i < MaximumCell ; ++i )
          cellDOS[i] = cnt ;
        break;
      default: // Assume atom.
        cellDOS[0] = qele  ;
        break;
    }

    // Electronic excited states.
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
	vector<double> tmpCellDOS(cellDOS);
    for (size_t j(0) ; j < eleExc.size() ; ++j){
      size_t nr = int(eleExc[j]) ;
      if (nr < MaximumCell) {
        for (size_t i(0) ; i < MaximumCell - nr ; i++ ) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i] ;
        }
      }
    }

    pDOS->setCellDensityOfStates(tmpCellDOS) ;

    return true;
  }

  // Calculate contribution to canonical partition function.
  double ClassicalRotor::canPrtnFnCntrb(gDensityOfStates* gdos, double beta) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);
    double sym = gdos->get_Sym();

    double qtot(1.0) ; 
    qtot *= double(gdos->getSpinMultiplicity());

	switch(rotorType){
      case NONLINEAR://3-D symmetric/asymmetric/spherical top
        qtot *= (sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/sym) ;
        break;
      case LINEAR://2-D linear
        qtot /= (rotConst[0]*sym*beta) ;
        break;
      default:
        break; // Assume atom.
    }

    // Electronic excited states.
    vector<double> eleExc;
    gdos->getEleExcitation(eleExc);
    for (size_t j(0) ; j < eleExc.size() ; ++j){
	  qtot += qtot*exp(-beta*eleExc[j]) ;
    }

	return qtot ;
  }
  
  // Function to return the number of degrees of freedom associated with this count.
  unsigned int ClassicalRotor::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    unsigned int nDOF(0) ;
    switch(rotorType){
      case NONLINEAR:
        nDOF = 3 ;
        break;
      case LINEAR:
        nDOF = 2 ;
        break;
      default:
        // Assume atom.
        break; 
    }

    return nDOF ;
  }

}//namespace
