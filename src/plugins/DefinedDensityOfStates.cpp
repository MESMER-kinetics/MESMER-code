
//-------------------------------------------------------------------------------------------
//
// DefinedDensityOfStates.cpp
//
// This file contains the implementation of the methods for importing and testing the 
// density of states calculated or determine by other means.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../gDensityOfStates.h"
#include "../Molecule.h"
#include "../Constants.h"
#include "../Spline.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class DefinedDensityOfStates : public DensityOfStatesCalculator
  {
  public:
    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Function to define particular counts of the DOS of a molecule.
		virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env) ;

    // Function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
		DefinedDensityOfStates(const char* id) : m_id(id),
			m_logSpline(true),
			m_grainSize(1.0),
      m_energy(), 
      m_DoS() 
    { Register(); }

    virtual ~DefinedDensityOfStates() {}
    virtual const char* getID()  { return m_id; }
    virtual bool includesRotations(){return true;}
    virtual DefinedDensityOfStates* Clone() { return new DefinedDensityOfStates(*this); }

  private:
    const char* m_id;
		bool m_logSpline;
		double m_grainSize;
		vector<double> m_energy;   // Energies of the states and their associated degeneracies.
    vector<double> m_DoS;

  } ;

  //-------------------------------------------------------------
  //Global instance, defining its id
  DefinedDensityOfStates theDefinedDensityOfStates("DefinedDensityOfStates");
  //-------------------------------------------------------------

  //Read data from XML and store in this instance.
  bool DefinedDensityOfStates::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
  {
    // Read in defined states data.

    PersistPtr pp = ppDOSC->XmlMoveTo("me:States") ;

    if (pp) {

      const char* p = pp->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";

			m_logSpline = !(pp->XmlReadBoolean("nologSpline"));

			// If number of states rather than density has been added it is necessary to 
			// know the energy interval for which the number of stated have been calculated. 
			double grainSize = pp->XmlReadDouble("grainSize", optional);
			if (!IsNan(grainSize) && grainSize > 0.0)
				m_grainSize = grainSize;

			while(pp = pp->XmlMoveTo("me:State"))
      {
        double energy = pp->XmlReadDouble("energy", optional);
        energy = getConvertedEnergy(units, energy);
				m_energy.push_back(energy) ;
				double DoS = pp->XmlReadDouble("degeneracy", optional);
				if (m_logSpline)
					DoS = log(max(DoS, 1.0));
				m_DoS.push_back(DoS) ;
      }

    } else {

      // Default : free rotor.

      cinfo << "No rotory-electronic states defined for " << gdos->getHost()->getName() << "." << endl;

			m_energy.push_back(0.0) ;
			m_DoS.push_back(1) ;
    }

    return true;
  }

	//
	// Calculate quantum mechanical 1D rotor densities of states of a free 
	// rotor and convolve them with the main density of states.
	//
	bool DefinedDensityOfStates::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
	{
		const size_t MaximumCell = env.MaxCell;

		vector<double> cellDOS(MaximumCell, 0.0);

		// Spline the sum of states.

		Spline spline;
		spline.Initialize(m_energy, m_DoS);

		double cellSize = env.CellSize;
		double maxEnergy = m_energy[m_energy.size() - 1];

		// Calculate the cell reaction flux from the spline.

		for (size_t i(0); i < MaximumCell; ++i) {
			double DensityOfStates = spline.Calculate(min(double(i)*cellSize, maxEnergy));
			if (m_logSpline) DensityOfStates = exp(DensityOfStates);
			cellDOS[i] = DensityOfStates/m_grainSize;
		}

		// Initialise density of states.   

		pDOS->setCellDensityOfStates(cellDOS);

		return true;
	}

  //
  // Provide a function to calculate contribution to canonical partition function.
  //
  void DefinedDensityOfStates::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne)
  {
    double Qrot(0.0), Erot(0.0), vErot(0.0) ;

    double zeroPointEnergy(m_energy[0]) ; 
    for (size_t i(0) ; i < m_energy.size() ; i++ ) {
      double ene = m_energy[i] - zeroPointEnergy ;
			double DoS = (m_logSpline) ? exp(m_DoS[i]) : m_DoS[i];
      Qrot  += DoS*exp(-beta*ene) ;
      Erot  += ene*DoS*exp(-beta*ene) ;
      vErot += ene*ene*DoS*exp(-beta*ene) ;
    }
    Erot /= Qrot ;
    vErot = vErot/Qrot - Erot*Erot ;

    PrtnFn   *= Qrot ;
    IntrlEne += Erot ;
    varEne   += vErot ;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int DefinedDensityOfStates::NoDegOfFreedom(gDensityOfStates* gdos) {

    unsigned int nDOF(1) ;

    return nDOF ;
  }

}//namespace
