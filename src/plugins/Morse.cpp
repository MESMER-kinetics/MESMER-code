
//-------------------------------------------------------------------------------------------
//
// Morse.cpp
//
// Author: Struan Robertson
// Date:   08/Jul/2012
//
// This file contains the implementation of the methods for calculating and testing the 
// density of states of a set of decoupled morse oscilators.
//
//-------------------------------------------------------------------------------------------

#include "../MolecularComponents.h"

using namespace std;
namespace mesmer
{
  class Morse : public DensityOfStatesCalculator
  {
  public:

    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell);

    // Function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Provide a function to calculate the zero point energy of a molecule.
	virtual double ZeroPointEnergy(gDensityOfStates* gdos) ;

	///Constructor which registers with the list of DensityOfStatesCalculators in the base class
    Morse(const char* id) : m_id(id) { Register(); }

    virtual ~Morse() {}
    virtual const char* getID()  { return m_id; }

    virtual Morse* Clone() { return new Morse(*this); }

  private :
    const char* m_id;
    PersistPtr m_ppDOSC ;
    vector<double> m_vibFreq ;  // The 0->1 transition of each Morse oscillator in cm-1.
    vector<double> m_anharmty ; // The associated anharmonicity.

  } ;

  //************************************************************
  //Global instance, defining its id
  Morse theMorse("Morse");
  //************************************************************
  // Read data from XML.
  // e.g. <me:MorseParameters vibrationalFrequency="3161.925" anharmonicity="-37.593"/>
  // All vib frequencies must have already appeared in <property dictRef="me:vibFreqs">
  // Any frequency not specifed here has anharmonicity = 0.0
  // The BeyerSwinehart method is removed when the Morse method is specified
  bool Morse::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC) {
    vector<double> vec;
    gdos->get_VibFreq(m_vibFreq); //Copy. These are scaled values
    double scale = gdos->get_scaleFactor();
    m_anharmty.assign(m_vibFreq.size(), 0.0);
        
    m_ppDOSC = ppDOSC ;
    PersistPtr pp = m_ppDOSC ;
    while(pp = pp->XmlMoveTo("me:MorseParameters")) {
      double vibFreq  = pp->XmlReadDouble("vibrationalFrequency", optional);
      double anharmty = pp->XmlReadDouble("anharmonicity",         true);
      for(unsigned i=0; i<=m_vibFreq.size(); ++i) {
        if(abs(m_vibFreq[i]/scale - vibFreq) < 0.1) { //correct back to raw frequency
          m_anharmty[i] = anharmty;
          break;
        }
      }
    }

    //Remove BeyerSwinehart
    return gdos->RemoveDOSCalculator("BeyerSwinehart");
  }

  // Provide a function to define particular counts of the DOS of a molecule.
  bool Morse::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {
    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    for (size_t nFrq(0) ; nFrq < m_vibFreq.size() ; nFrq++ )
    {
      double vibFreq  = m_vibFreq[nFrq];
      double anharmty = m_anharmty[nFrq];

      // Maximum bound energy.

      int nmax(0) ;
      if (anharmty < 0.0) {
        nmax = int(-0.5*(vibFreq + anharmty)/anharmty)  ;
      } else {
        nmax = int(double(MaximumCell)/vibFreq) ;
      }

      // Calculate energy levels.

      vector<double> energyLevels ;
      for (int n(0) ; n <= nmax ; n++ ) {
        double nu = double(n) ;
        double energy = nu*vibFreq + nu*(nu + 1)*anharmty ;
        energyLevels.push_back(energy) ;
      } 

      // Convolve with the density of states for the other degrees of freedom.
      // (Essentially the Stein-Rabinovitch algorithm).

      vector<double> tmpCellDOS(cellDOS) ;
      for (size_t k(1) ; k < energyLevels.size() ; k++ ) {
        size_t nr = nint(energyLevels[k]) ;
        if (nr < MaximumCell) {
          for (size_t i(0) ; i < MaximumCell - nr ; i++ ) {
            tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i] ;
          }
        }
      }

      cellDOS = tmpCellDOS ;
    }

    // Replace existing density of states.   

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;
  }

  // Calculate contribution to canonical partition function.
  double Morse::canPrtnFnCntrb(gDensityOfStates* gdos, double beta) {

    double qtot(1.0) ; 
    for (size_t nFrq(0) ; nFrq < m_vibFreq.size() ; nFrq++ )
    {
      double vibFreq  = m_vibFreq[nFrq];
      double anharmty = m_anharmty[nFrq];

      // Maximum bound energy.

      int nmax(0) ;
      if (anharmty < 0.0) {
        nmax = int(-0.5*(vibFreq + anharmty)/anharmty)  ;
      } else {
        nmax = int(15.0*log(10.0)/(beta*vibFreq)) ;
      }

      // Calculate canonical partition function.

      double qtmp(0.0) ;
      for (int n(0) ; n <= nmax ; n++ ) {
        double nu = double(n) ;
        double energy = nu*vibFreq + nu*(nu + 1)*anharmty ;
        qtmp += exp(-beta*energy) ;
      }

      qtot *= qtmp ;
    }

    return qtot ;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int Morse::NoDegOfFreedom(gDensityOfStates* gdos) {
    return m_vibFreq.size() ;
  }

  // Provide a function to calculate the zero point energy of a molecule.
  double Morse::ZeroPointEnergy(gDensityOfStates* gdos) {

	double ZPE(0.0) ;
    for (size_t j(0) ; j < m_vibFreq.size() ; ++j ) {
      ZPE += m_vibFreq[j] + m_anharmty[j]*0.5 ;
    }
    return ZPE*0.5 ;
  }

}//namespace
