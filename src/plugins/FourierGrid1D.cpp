//-------------------------------------------------------------------------------------------
//
// FourierGrid1D.cpp
//
// Authors: Struan Robertson (based, in part, on earlier work by Chi-Hsiu Liang).
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
// The program uses an even number of grid points. The Hamiltonian matrix elements are
// calculated using analytic formula described in the second of the above references.
//
// The analytical Hamiltonian expression given in the reference contains a small error.
// The formula should read:
//
//            H(i,j) = {(h**2)/(4*m*(L**2)} *
//                      [(N-1)(N-2)/6 + (N/2)] + V(Xi),   if i=j
//
//            H(i,j) = {[(-1)**(i-j)] / m } *
//                      { h/[2*L*sin(M_PI*(i-j)/N)]}**2 ,     if i#j
//
// The eigenvalues of the Hamiltonian matrix which lie below the asymptotic (large x) value
// of V(x) correspond to the bound state energies. The corresponding eigenvectors are the
// representation of the bound state wavefunctions on a regular one dimensional grid.
//
// Note this method is the conjugate to that used for the hindered rotors, where the basis
// is taken to be eigenfunctions of the kinetic energy operator.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"
#include "../Molecule.h"
#include "../Constants.h"
#include "../gDensityOfStates.h"
#include "../gStructure.h"
#include "Potential1D.h"

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
			m_ppotential(NULL),
			m_nGridPoints(100),
      m_ZPE(0.0),
      m_reducedMass(0.0),
      m_energyLevels(),
      m_plotStates(false),
      m_writeStates(false) {
      Register();
    }

    virtual const char* getID() { return m_id; }

		virtual ~FourierGrid1D() { delete m_ppotential; }
    virtual FourierGrid1D* Clone() { return new FourierGrid1D(*this); }
    virtual const char* Description() {
      return
        "Solves one dimensional Schrodinger equation for bound state\n"
        "  eigenvalues and eigenfunctions corresponding to a potential V(x).\n  "
        "\n";
    };

  private:

    // Provide data for plotting states against potential.
    void outputPlotData(vector<double> &abscissa, vector<double> &potential) const;

    // Print the hindered rotor states.
    void outputStateData() const;

    const char* m_id;

    Potential1D* m_ppotential;       // Class holding the potential representation.
		size_t m_nGridPoints;            // Number of gridpoints to be used to construct Hamiltonian.

    double m_ZPE;                    // Zero point energy.
    double m_reducedMass;            // reduced mass.
    string m_format;                 // Format that the potential is represented in.

    vector<double> m_energyLevels;	 // The energies of the hindered rotor states.

    bool m_plotStates;               // If true output data for plotting. 
    bool m_writeStates;              // If true energy levels written to output. 
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

		// Read the number of grid points it be used and make sure it is even. 

		m_nGridPoints  = ppDOSC->XmlReadInteger("me:NumGridPnts", optional);
		m_nGridPoints += (m_nGridPoints % 2) ? 1 : 0;

    // Determine reduced moment mass. These can either be read in or 
    // calculated if the coordinates are available. Always take the 
    // user defined reduced mass if present, if not then calculate. 

		PersistPtr ppRM = ppDOSC->XmlMoveTo("me:reducedMass");
    if (!ppRM)
      ppRM = ppDOSC->XmlMoveTo("reducedMass");
    const char* bondref = ppDOSC->XmlReadValue("bondRef", optional);
    if (!bondref)
      bondref = ppDOSC->XmlReadValue("me:bondRef", optional);

    // The following vector, "mode", will be used to hold the internal 
    // motion vector (analogus to the the internal rotor case). 
    vector<double> mode(3 * gs.NumAtoms(), 0.0);
    if (ppRM) {
      const char* p = ppRM->XmlReadValue("units", optional);
      string units = p ? p : "amu";
      m_reducedMass = ppDOSC->XmlReadDouble("me:reducedMass", optional);
    }
    else if (bondref && *bondref != '\0') {
      // Read bond IDs
      stringstream iss(bondref);
      vector<string> bondIDs((istream_iterator<std::string>(iss)), istream_iterator<std::string>());

      size_t nBonds = bondIDs.size();

      // Determine type of interaction

      switch (nBonds) {
      case 1:
        m_reducedMass = gs.bondStretchReducedMass(bondIDs, mode);
        break;
      case 2:
        m_reducedMass = gs.angleBendReducedMass(bondIDs, mode);
        break;
      case 3:
        m_reducedMass = gs.inversionReducedMass(bondIDs, mode);
        break;
      default:
        cerr << "Bond definition for FGH term in " << gdos->getHost()->getName() << " incorrectly defined." << endl;
        return false;
      }

    }
    else {
      cerr << "No reduced mass found or means to calculate it." << endl;
      return false;
    }

    // Read in potential information.

    PersistPtr pp = ppDOSC->XmlMoveTo("me:vibrationalPotential");

    if (pp) {

      const char* p = pp->XmlReadValue("format", true);
      m_format = string(p);

      if (m_format == "analytical") {
        m_ppotential = new AnalyticalPotential(m_nGridPoints);
      }
      else if (m_format == "numerical") {
        m_ppotential = new NumericalPotential(m_nGridPoints);
      }
      else {
				cerr << "__FUNCTION__: Unknown potential type." << endl;
				return false;
			}
      m_ppotential->InitializePotential(pp);
    }

    // Check if there is a Hessian and knock out the frequency
    // associated with this mode.

    if (gdos->hasHessian()) {
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
  // Calculate quantum mechanical 1D densities of states of a
  // vibration and convolve them with the main density of states.
  //
  bool FourierGrid1D::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
    const size_t MaximumCell = env.MaxCell;

    vector<double> cellDOS;
    if (!pDOS->getCellDensityOfStates(cellDOS, false)) // Retrieve the DOS vector without re-calculating.
      return false;

    vector<double> tmpCellDOS(cellDOS);

    // Find maximum quantum No. for rotor. To ensure convergence basis functions 
    // that span a range twice that request are used with a minimum of 100000
    // wavenumbers which appears in the definition of root below. 

    vector<double> abscissa;
    vector<double> potential;
    m_ppotential->calculatePotential(abscissa, potential);

    dMatrix hamiltonian(m_nGridPoints);

    // Add kinetic energy terms.
		const double cfctr = M_PI / double(m_nGridPoints);
		const double chlen = m_ppotential->get_characteristicLength();
		const double tfctr = PlancksConstant_in_JouleSecond * AvogadroC *1.e+23 / (4.0*SpeedOfLight_in_cm*m_reducedMass*chlen*chlen);
		const double dgnEl = tfctr*((double(m_nGridPoints)*double(m_nGridPoints) + 2.0) / 6.0);
		for (int i(0) ; i < int(m_nGridPoints) ; i++) {
			hamiltonian[i][i] = dgnEl;
      for (int j(i+1) ; j < int(m_nGridPoints) ; j++) {
				int dij     = i - j;
				double tmp  = sin(double(dij)*cfctr) ;
				double dsgn = (dij) % 2 ? -1.0 : 1.0;
				hamiltonian[i][j]  = tfctr*dsgn /(tmp * tmp) ;
				hamiltonian[j][i]  = hamiltonian[i][j];
      }
    }

		// Add diagonal potential terms.

		for (size_t n(0); n < m_nGridPoints; n++) {
			hamiltonian[n][n] += potential[n];
		}

		// Now diagonalize hamiltonian matrix to determine energy levels.
    vector<double> eigenvalues(m_nGridPoints, 0.0);

    hamiltonian.diagonalize(&eigenvalues[0]);

    // Save energy levels for partition function calculations.

    m_energyLevels = eigenvalues;

    // Shift eigenvalues by the zero point energy and convolve with the 
    // density of states for the other degrees of freedom.

    m_ZPE = m_energyLevels[0];
    for (size_t k(1); k < m_energyLevels.size() ; k++) {
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
      outputPlotData(abscissa, potential);

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
  void FourierGrid1D::outputPlotData(vector<double> &abscissa, vector<double> &potential) const {

    ctest << endl << "Vibrational state data for plotting." << endl << endl;

    size_t npoints(abscissa.size());
    for (size_t i(0); i < npoints; ++i) {
    	ctest << formatFloat(abscissa[i], 6, 15) << ", " << formatFloat(potential[i], 6, 15) << endl;
    }
    ctest << endl;

    for (size_t i(0); i < m_energyLevels.size(); i++) {
    	ctest << formatFloat(abscissa[0], 6, 15) << ", " << formatFloat(m_energyLevels[i], 6, 15) << endl;
    	ctest << formatFloat(abscissa[npoints-1], 6, 15) << ", " << formatFloat(m_energyLevels[i], 6, 15) << endl;
    	ctest << endl;
    }

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
