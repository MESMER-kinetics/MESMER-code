//-------------------------------------------------------------------------------------------
//
// FourierGrid1D.cpp
//
// Authors: Struan Robertson, Chi-Hsiu Liang
// Date:    16/Oct/2016
//
// This file contains the implementation of the methods for calculating and testing the 
// density of states of a one dimensional fourier grid model.
//
// This method solves one dimensional Schrodinger equation for bound state eigenvalues
// and eigenfunctions corresponding to a potential V(x).
//
// This method is fully described in:
// C.C. Marston and G.G.Balint-Kurti, J.Chem. Phys.,91,3571(1989).
// and
// G.G. Balint-Kurti, R.N. Dixon and C.C. Marston, Internat. Rev.
// Phys. Chem.,11, 317 (1992).
//
// The eigenvalues of the Hamiltonian matrix which lie below the
// asymptotic (large x) value of V(x) correspond to the bound state
// energies. The corresponding eigenvectors are the representation
// of the bound state wavefunctions on a regular one dimensional grid.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"
#include "../Molecule.h"
#include "../Constants.h"
#include "../gDensityOfStates.h"
#include "../gStructure.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
	class FourierGrid1D : public DensityOfStatesCalculator
	{
	public:
		//Read data from XML. 
		virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC = NULL);

		// Provide a function to define particular counts of the DOS of a molecule.
		virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env);

		// Provide a function to calculate contribution to canonical partition function.
		virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne);

		// Function to return the number of degrees of freedom associated with this count.
		virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) { return 1; };

		// Provide a function to calculate the zero point energy of a molecule.
		virtual double ZeroPointEnergy(gDensityOfStates* gdos) { return m_ZPE; };

		// Constructor which registers with the list of DensityOfStatesCalculators in the base class
		// This class is an extra DOS class.
		FourierGrid1D(const char* id) : m_id(id),
			m_ZPE(0.0),
			m_reducedMass(0.0),
			m_x(),
			m_energyLevels(),
			m_potential(),
		  m_plotStates(false),
			m_writeStates(false) { Register(); }

		virtual const char* getID() { return m_id; }

		virtual ~FourierGrid1D() {}
		virtual FourierGrid1D* Clone() { return new FourierGrid1D(*this); }

	private:

		// Provide data for plotting states against potential.
		void outputPlotData() const;

		// Print the hindered rotor states.
		void outputStateData() const;

		const char* m_id;

		double m_ZPE;                   // Zero point energy.
		double m_reducedMass;           // reduced mass.

		vector<double> m_energyLevels;	// The energies of the hindered rotor states.
		vector<double> m_x;           	// Linear coordinate.
		vector<double> m_potential;	    // Potential associated with linear coordinate.

		bool m_plotStates;              // If true output data for plotting. 
		bool m_writeStates;             // If true energy levels written to output. 
		PersistPtr m_ppConfigData;
	};

	//-------------------------------------------------------------
	//Global instance, defining its id
	FourierGrid1D theFourierGrid1D("FourierGrid1D");
	//-------------------------------------------------------------

	using OpenBabel::vector3;
	//Read data from XML and store in this instance.
	bool FourierGrid1D::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
	{
		gStructure& gs = gdos->getHost()->getStruc();
		if (!gs.ReadStructure())
		{
			cerr << "A complete set of atom coordinates are required for hindered rotor calculations" << endl;
			return false;
		}

		//Remove the vibrational frequency that this hindered rotation replaces
		const char* vibFreq = ppDOSC->XmlReadValue("me:replaceVibFreq", optional);
		if (vibFreq)
		{
			if (!gdos->removeVibFreq(atof(vibFreq)))
			{
				cerr << "Cannot find vibrational frequency " << vibFreq << " to replace it with hindered rotor" << endl;
				return false;
			}
			cinfo << " replacing vib freq " << vibFreq;
		}
		cinfo << endl;

		// Determine reduced moment of inertia function. These can either 
		// be read in or calculated if the coordinates are available.

		PersistPtr pp = ppDOSC->XmlMoveTo("me:InternalRotorInertia");

		if (pp) {
			const char* p = pp->XmlReadValue("units", optional);
			string units = p ? p : "cm-1";

			vector<int> indicies;
			vector<double> coefficients;
			int maxIndex(0);
			while (pp = pp->XmlMoveTo("me:InertiaPoint"))
			{
				int index = pp->XmlReadInteger("index", optional);
				indicies.push_back(index);
				maxIndex = max(maxIndex, index);

				double coefficient = pp->XmlReadDouble("coefficient", optional);
				if (IsNan(coefficient))
					coefficient = 0.0;
				coefficients.push_back(coefficient);
			}

			if (units == "amuA^2") {
			}

		}
		else {

			// No inertia terms red in so let us see if they are to be calculated.
			// If so, set the appropriate flags.

			pp = ppDOSC->XmlMoveTo("me:CalculateInternalRotorInertia");
			if (pp) {
			}
			else {

				// No information supplied, so simply calculate reduced moment of inertia
				// for current configuration.

			}
		}

		// Read in potential information.

		pp = ppDOSC->XmlMoveTo("me:Potential");

		if (pp) {

			const char* p = pp->XmlReadValue("format", true);
			string format(p);

			p = pp->XmlReadValue("units", optional);
			string units = p ? p : "kJ/mol";

			// Numerical potential.

			vector<double> potential;
			vector<double> angle;

			while (pp = pp->XmlMoveTo("me:PotentialPoint"))
			{
				double anglePoint = pp->XmlReadDouble("angle", optional);
				if (IsNan(anglePoint))
					anglePoint = 0.0;
				angle.push_back(anglePoint);

				double potentialPoint = pp->XmlReadDouble("potential", optional);
				if (IsNan(potentialPoint))
					potentialPoint = 0.0;
				potentialPoint = getConvertedEnergy(units, potentialPoint);
				potential.push_back(potentialPoint);
			}

		}

		// Check if there is a Hessian and knock out the frequency
		// associated with this mode.

		if (gdos->hasHessian()) {
			// The following vector, "mode", will be used to hold the internal rotation 
			// mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010). 
			vector<double> mode(3 * gs.NumAtoms(), 0.0);
			// gs.internalRotationVector(get_BondID(), mode);
			if (!gdos->projectMode(mode)) {
				cerr << "Failed to project out internal rotation." << endl;
				return false;
			}
		}

		// Check if is data for plotting are required.

		pp = ppDOSC->XmlMoveTo("me:PlotStates");
		if (pp) {
			m_plotStates = true;
		}

		// Check if is energy level values are to be written.

		pp = ppDOSC->XmlMoveTo("me:WriteStates");
		if (pp) {
			m_writeStates = true;
		}

		// Check if configuration data are required.

		pp = ppDOSC->XmlMoveTo("me:ConfigurationalData");
		if (pp) {
			gs.set_Verbose(true);
			m_ppConfigData = pp;
		}

		return true;
	}

	//
	// Calculate quantum mechanical 1D rotor densities of states of an 
	// internal rotor and convolve them with the main density of states.
	//
	bool FourierGrid1D::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
	{
		const size_t MaximumCell = env.MaxCell;

		vector<double> cellDOS;
		if (!pDOS->getCellDensityOfStates(cellDOS, false)) // Retrieve the DOS vector without re-calculating.
			return false;

		vector<double> tmpCellDOS(cellDOS);

		// Calculate the (angle dependent) internal moment of inertia. This
		// has to be done at this point in order to capture any interaction
		// with other internal rotors.

		//if (calIntrlIrt) {
		//	gStructure& gs = pDOS->getHost()->getStruc();
		//	size_t nAngle(36);
		//	vector<double> angle(nAngle, 0.0), redInvMOI;
		//	gs.reducedMomentInertiaAngular(get_BondID(), get_Phase(), angle, redInvMOI, m_ppConfigData);  // Units a.u.*Angstrom*Angstrom.
		//	FourierCosCoeffs(angle, redInvMOI, m_kineticCosCoeff, get_Expansion());
		//	FourierSinCoeffs(angle, redInvMOI, m_kineticSinCoeff, get_Expansion());
		//}

		// Find maximum quantum No. for rotor. To ensure convergence basis functions 
		// that span a range twice that request are used with a minimum of 100000
		// wavenumbers which appears in the definition of root below. 

		int lmax = int(100);
		size_t nstates = lmax;

		dMatrix hamiltonian(nstates);

		// Add diagonal potential terms first.

		for (int n(0); n <= nstates; n++) {
			hamiltonian[n][n] = m_potential[n];
		}

		// Add off-diagonal kinetic terms.

			for (int l(0); l < lmax ; l++) {
				double Tl = double(l*l)/(2.0*m_reducedMass);
				for (size_t i(0); i < nstates; i++) {
					for (size_t j(i+1); j < nstates; j++) {
						hamiltonian[i][j] += Tl*cos(l*(i-j)) ;
						hamiltonian[j][i] += hamiltonian[i][j];
					}
				}
			}

		// Now diagonalize hamiltonian matrix to determine energy levels.
		vector<double> eigenvalues(nstates, 0.0);

		hamiltonian.diagonalize(&eigenvalues[0]);

		// Save energy levels for partition function calculations.

		m_energyLevels = eigenvalues;

		// Shift eigenvalues by the zero point energy and convolve with the 
		// density of states for the other degrees of freedom.

		m_ZPE = m_energyLevels[0];
		for (size_t k(1); k < nstates; k++) {
			size_t nr = nint((m_energyLevels[k] - m_ZPE) / env.CellSize);
			if (nr < MaximumCell) {
				for (size_t i(0); i < MaximumCell - nr; i++) {
					tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i];
				}
			}
		}

		// Replace existing density of states.   

		pDOS->setCellDensityOfStates(tmpCellDOS);

		// If required, created graphical date.
		if (m_plotStates)
			outputPlotData();

		// If required, created graphical date.
		if (m_writeStates)
			outputStateData();

		return true;

	}

	//
	// Provide a function to calculate contribution to canonical partition function.
	// (Mostly for testing purposes.)
	//
	void FourierGrid1D::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne)
	{
		double Qintrot(0.0), Eintrot(0.0), varEintrot(0.0);

		double zeroPointEnergy(m_energyLevels[0]);
		for (size_t i(0); i < m_energyLevels.size(); i++) {
			double ene = m_energyLevels[i] - zeroPointEnergy;
			Qintrot += exp(-beta*ene);
			Eintrot += ene*exp(-beta*ene);
			varEintrot += ene*ene*exp(-beta*ene);
		}
		Eintrot /= Qintrot;
		varEintrot = varEintrot / Qintrot - Eintrot*Eintrot;

		PrtnFn *= Qintrot;
		IntrlEne += Eintrot;
		varEne += varEintrot;
	}

	// Provide data for plotting states against potential.
	void FourierGrid1D::outputPlotData() const {

		ctest << endl << "Vibrational state data for plotting." << endl << endl;
		//int npoints(500);
		//double dAngle = M_PI / double(npoints);
		//for (int i(-npoints); i < npoints; ++i) {
		//	double angle = double(i)*dAngle;
		//	double potential = CalculatePotential(angle);
		//	ctest << formatFloat(angle, 6, 15) << ", " << formatFloat(potential, 6, 15) << endl;
		//}
		//ctest << endl;

		//for (size_t i(0); i < m_energyLevels.size(); i++) {
		//	ctest << formatFloat(-M_PI, 6, 15) << ", " << formatFloat(m_energyLevels[i], 6, 15) << endl;
		//	ctest << formatFloat(M_PI, 6, 15) << ", " << formatFloat(m_energyLevels[i], 6, 15) << endl;
		//	ctest << endl;
		//}

	}

	// Print the vibrational states.
	void FourierGrid1D::outputStateData() const {

		ctest << "\nFourier grid vibrational states (cm-1).\n" << endl;
		for (size_t i(0); i < m_energyLevels.size(); i++) {
			ctest << formatFloat(m_energyLevels[i], 6, 15) << endl;
		}
		ctest << endl;
	}

};
