#ifndef GUARD_QMHinderedRotorPotential_h
#define GUARD_QMHinderedRotorPotential_h

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

namespace mesmer
{
  class QMHinderedRotorPotential : public DensityOfStatesCalculator
  {
  public:

    //Read data from XML. Some is stored hear and some in a MolecularComponent class.
    virtual bool ReadParameters(Molecule* pMol, PersistPtr ppDOSC=NULL);

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell);

    // Provide a function to calculate contribution to canonical partition function.
    // (Mostly for testing purposes.)
    virtual double canPrtnFnCntrb(const double beta) { return 1.0 ;}
    
    virtual QMHinderedRotorPotential* Clone() { return new QMHinderedRotorPotential(*this); }
    
    ///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    //This class is an extra DOS class: there needs to be a non-extra DensityOfStatesCalculator class
    QMHinderedRotorPotential(const std::string& id) : DensityOfStatesCalculator(id, true){}

    virtual ~QMHinderedRotorPotential() {}

  private:
    std::string     m_bondID;
   	vector<double>  ak; // Fourier coefficients for sine functions
		vector<double>  bk; // Fourier coefficients for cosine functions
		double a0; 					// 
    
		//store the PES in a Fourier series. According to the fineness of the 
    //data set, the complexity of the PES can be achieved.
      
    int             m_numberGridPoint;
    double          m_vibFreq;
    double          m_reducedMomentInertia;

  } ;

}//namespace

#endif // GUARD_QMHinderedRotorPotential_h
