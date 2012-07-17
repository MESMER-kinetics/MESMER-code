
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

#include "../DensityOfStates.h"
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

		///Constructor which registers with the list of DensityOfStatesCalculators in the base class
		//This class calculates a complete DOS: it is not an extra class. 
		Morse(const std::string& id) : DensityOfStatesCalculator(id, false){}

		virtual ~Morse() {}
		virtual Morse* Clone() { return new Morse(*this); }

	private :

		PersistPtr m_ppDOSC ;

	} ;

	//************************************************************
	//Global instance, defining its id (usually the only instance) but here with an alternative name
	Morse theMorse("Morse");
	//************************************************************

	//Read data from XML. 
	bool Morse::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC) {
		m_ppDOSC = ppDOSC ;
		return true ;
	}

	// Provide a function to define particular counts of the DOS of a molecule.
	bool Morse::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
	{

		PersistPtr pp = m_ppDOSC ;
		while(pp = pp->XmlMoveTo("me:PotentialPoint"))
		{
		}

		vector<double> VibFreq ;
		pDOS->get_VibFreq(VibFreq) ;

		vector<double> cellDOS;
		if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
			return false;

		// Implementation of the Beyer-Swinehart algorithm.
		for (size_t j(0) ; j < VibFreq.size() ; ++j ) {
			size_t freq = static_cast<size_t>(VibFreq[j]) ;
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
	double Morse::canPrtnFnCntrb(gDensityOfStates* gdos, double beta) {

		double qtot(1.0) ; 
		vector<double> vibFreq; 
		gdos->get_VibFreq(vibFreq);
		for (size_t j(0) ; j < vibFreq.size() ; ++j ) {
			qtot /= (1.0 - exp(-beta*vibFreq[j])) ;
		}

		return qtot ;
	}

	// Function to return the number of degrees of freedom associated with this count.
	unsigned int Morse::NoDegOfFreedom(gDensityOfStates* gdos) {

		vector<double> vibFreq; 
		gdos->get_VibFreq(vibFreq);

		return vibFreq.size() ;
	}

}//namespace
