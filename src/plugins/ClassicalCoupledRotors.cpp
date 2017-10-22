
//-------------------------------------------------------------------------------------------
//
// ClassicalCoupledRotors.cpp
//
// Author: Struan Robertson
// Date:   07/Aug/2017
//
// This file contains the implementation of the methods for calculating and testing the 
// density of states of a set of internal rotors coupled with each other and overall rotation.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../Molecule.h"
#include "../gDensityOfStates.h"
#include "../gStructure.h"
#include "../Sobol.h"
#include "HinderedRotorUtils.h"

using namespace std;
namespace mesmer
{
  class ClassicalCoupledRotors : public HinderedRotorUtils
  {
  public:

    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC = NULL);

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env);

    // Function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne);

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos);

    // Constructor which registers with the list of DensityOfStatesCalculators in TopPlugin.
    // This class calculates a complete DOS: it is not an extra class. 
    ClassicalCoupledRotors(const char* id) : m_id(id), 
			HinderedRotorUtils(id),
			m_Sym(1.0),
			m_OpticalSym(1.0), 
			m_SpinMultiplicity(1), 
			m_MCPnts(100), 
			m_bondIDs(),
			m_periodicities(), 
			m_knmtcFctr(), 
			m_potential(),
			m_potentialCosCoeffs(),
			m_potentialSinCoeffs() { Register(); }

    virtual ~ClassicalCoupledRotors() {}
    virtual const char* getID() { return m_id; }

    // Included only in a subset of DensityOfStatesCalculators.
    // Otherwise the baseclass function returns false.
    virtual bool includesRotations() { return true; }

    virtual ClassicalCoupledRotors* Clone() { return new ClassicalCoupledRotors(*this); }

  private:

		// Calculates the potential energy associated with torsions.
		double TorsionalPotential(const vector<double>& angles);

    const char* m_id;

    double m_Sym;              // Rotational Symmetry Number.
    double m_OpticalSym;       // Transition states may have 2 optical isomers, which
                               // is accounted for by an additionsl symmetry number.
    size_t m_SpinMultiplicity; // Spin multiplicity.
		size_t m_MCPnts;           // Number of Monte-Carlo points to use.

		vector<string> m_bondIDs;       // The IDs of the binds to be used in the coupled model.
		vector<size_t> m_periodicities; // The periodicities of the hindered rotors.

		vector<double> m_knmtcFctr; // Configuration kinematic factor.
		vector<double> m_potential; // Torsion configuration potential.

		vector<vector<double> > m_potentialCosCoeffs; // Torsion potential cosine Coefficients.
		vector<vector<double> > m_potentialSinCoeffs; // Torsion potential sine Coefficients.

	};

  //************************************************************
  //Global instance, defining its id
  ClassicalCoupledRotors theClassicalCoupledRotors("ClassicalCoupledRotors");
  //************************************************************

  //Read data from XML and store in this instance.
  bool ClassicalCoupledRotors::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
  {
    PersistPtr pp = gdos->getHost()->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    // Spin multiplicity (Note spin can be defined as an attribute on a molecule).
    m_SpinMultiplicity = ppPropList->XmlReadPropertyInteger("me:spinMultiplicity", optional);
    if (m_SpinMultiplicity == 0)
      m_SpinMultiplicity = pp->XmlReadInteger("spinMultiplicity");

    // Rotational Symmetry Number
    m_Sym = ppPropList->XmlReadPropertyDouble("me:symmetryNumber");

    // Transition states may have 2 optical isomers, which
    // is accounted for by an additionsl symmetry number.
    if (gdos->getHost()->isMolType("transitionState")) {
      m_OpticalSym = ppPropList->XmlReadPropertyDouble("me:TSOpticalSymmetryNumber", optional);
      if (!IsNan(m_OpticalSym)) {
        if (m_OpticalSym > 2.0) {
          string name(gdos->getHost()->getName());
          cinfo << "Warning: The transition state " << name << " has an optical symmetry number greater than 2." << endl;
        }

        // Adjust Rotational Symmetry Number by Optical Symmetry Number if necessary.
        if (m_OpticalSym > 1.0)
          m_Sym /= m_OpticalSym;
      }
    }

		// Check there is structural data present.
		gStructure& gs = gdos->getHost()->getStruc();
		if (!gs.ReadStructure())
		{
			cerr << "A complete set of atom coordinates are required for a coupled hindered rotor calculations" << endl;
			return false;
		}

		// Read in rotor information.
		PersistPtr ppDOS = pp->XmlMoveTo("me:DOSCMethod");
		int MCpnts = ppDOS->XmlReadInteger("me:MCPoints", optional);
		PersistPtr ppRotor(ppDOS->XmlMoveTo("me:RotorArray"));
		while (ppRotor = ppRotor->XmlMoveTo("me:Rotor")) {
			const char* bondID = ppRotor->XmlReadValue("bondRef");
			if (!bondID || *bondID == '\0') {
				cerr << "No <bondRef> specified for the coupled hindered rotating bond" << endl;
				return false;
			}
			else {
				m_bondIDs.push_back(string(bondID));
				gs.addRotBondID(string(bondID));
			}
			int periodicity = ppRotor->XmlReadInteger("periodicity", optional);
			m_periodicities.push_back((periodicity > 0) ? periodicity:1); // set periodicity to unity by default.

			// Read in potential parameters. 
			vector<double> potentialCosCoeff;
			vector<double> potentialSinCoeff;
			ReadPotentialParameters(ppRotor, string(bondID), potentialCosCoeff, potentialSinCoeff);
      m_potentialCosCoeffs.push_back(potentialCosCoeff);
      m_potentialSinCoeffs.push_back(potentialSinCoeff);

			// Check if there is a Hessian and knock out the frequency
			// associated with this internal rotation.

			if (gdos->hasHessian()) {
				// The following vector, "mode", will be used to hold the internal rotation 
				// mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010). 
				vector<double> mode(3 * gs.NumAtoms(), 0.0);
				gs.internalRotationVector(string(bondID), mode);
				if (!gdos->projectMode(mode)) {
					cerr << "Failed to project out internal rotation." << endl;
					return false;
				}
			}
		}

		// Check if configuration data are required.

		pp = ppDOS->XmlMoveTo("me:ConfigurationalData");
		if (pp) {
			gs.set_Verbose(true);
			// m_ppConfigData = pp;
		}

		return true;
  }

  // Provide a function to define particular counts of the DOS of a molecule.
  bool ClassicalCoupledRotors::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
    const size_t MaximumCell = env.MaxCell;
    const double cellSize = env.CellSize;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellSize, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0);

		gStructure& gs = pDOS->getHost()->getStruc();

		// Instantiate a random vector generator.
		Sobol sobol;

		// Configuration loop.
		long long seed(1);
		double twoPi = 2.0*M_PI;
		m_knmtcFctr.resize(m_MCPnts, 0.0);
		m_potential.resize(m_MCPnts, 0.0);
		for (size_t i(0); i < m_MCPnts; ++i) {

			// Select angular coordinates.
			vector<double> angles(m_bondIDs.size(), 0.0);
			sobol.sobol(angles.size(), &seed, angles);
			for (size_t j(0); j < angles.size(); j++)
				angles[j] *= twoPi;

			// Calculate the determinant of the Wilson G Matrix.
			double det = gs.getGRITDeterminant(angles);  // Units a.u.*Angstrom*Angstrom.
			m_knmtcFctr[i] = sqrt(fabs(det));

			// Calculate potential energy.
			m_potential[i] = TorsionalPotential(angles) ;
		}

		// Heavyside function integration.
		for (size_t i(0); i < m_MCPnts; ++i) {
			double kfctr = m_knmtcFctr[i];
			double ptnl  = m_potential[i];
			size_t ll = size_t(ptnl / cellSize) ;
			for (size_t j(ll); j < cellDOS.size(); ++j) {
				cellDOS[j] += kfctr;
			}
		}

		// Conversion, spin and symmetry number factor.
		double cnt = size_t(m_SpinMultiplicity)/pow(conMntInt2RotCnt, 0.5*double(m_bondIDs.size() + 3));
		size_t SymNo = size_t(m_Sym);
		for (size_t j(0); j < m_bondIDs.size(); ++j)
			SymNo *= m_periodicities[j];
		cnt /= double(SymNo);
		cnt /= double(m_MCPnts);
		cnt /= ((m_bondIDs.size() > 1) ? MesmerGamma(0.5*double(m_bondIDs.size() + 1)) : 1.0);
		for (size_t j(0); j < cellDOS.size(); ++j) {
			cellDOS[j] *= cnt;
		}

		// Convolve with remaining energy contributions.
		vector<double> tmpCellDOS(cellDOS.size(), 0.0);
		double pwr = 0.5*(m_bondIDs.size() - 1);
		for (size_t j(0); j < tmpCellDOS.size(); ++j) {
			tmpCellDOS[j] += pow(cellEne[j], pwr);
		}
		vector<double> hndrRtrDOS(cellDOS.size(), 0.0);
		FastLaplaceConvolution(cellDOS, tmpCellDOS, hndrRtrDOS);

    // Electronic excited states.
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
    tmpCellDOS = hndrRtrDOS;
    for (size_t j(0); j < eleExc.size(); ++j) {
      size_t nr = nint(eleExc[j] / cellSize);
      if (nr < MaximumCell) {
        for (size_t i(0); i < MaximumCell - nr; i++) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i];
        }
      }
    }

    pDOS->setCellDensityOfStates(tmpCellDOS);

    return true;
  }

  // Calculate contribution to canonical partition function.
  void  ClassicalCoupledRotors::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) {

		// Conversion, spin and symmetry number factor.
		double cnt = size_t(m_SpinMultiplicity) / pow(conMntInt2RotCnt, 0.5*double(m_bondIDs.size() + 3));
		size_t SymNo = size_t(m_Sym);
		for (size_t j(0); j < m_bondIDs.size(); ++j)
			SymNo *= m_periodicities[j];
		cnt /= double(SymNo);
		cnt /= double(m_MCPnts);

		// Heavyside function integration.
		double qtot(0.0), ene(0.0), var(0.0);
		for (size_t i(0); i < m_MCPnts; ++i) {
			const double ptnl = m_potential[i];
			const double tmp  = m_knmtcFctr[i]*exp(-beta*ptnl);
			qtot += tmp;
			ene  += ptnl*tmp;
			var  += ptnl*ptnl*tmp;
		}
		var  /= qtot;
		ene  /= qtot;
		var  -= ene*ene ;
		var  += 0.5*double(m_bondIDs.size() + 3) / (beta*beta);
		ene  += 0.5*double(m_bondIDs.size() + 3) / beta ;
		qtot *= cnt / pow(beta, 0.5*double (m_bondIDs.size() + 3));

    // Electronic excited states.
    vector<double> eleExc;
    gdos->getEleExcitation(eleExc);
    double qelec(1.0), Eelec(0.0), varEelec(0.0);
    for (size_t j(0); j < eleExc.size(); ++j) {
      qelec += exp(-beta*eleExc[j]);
      Eelec += eleExc[j] * exp(-beta*eleExc[j]);
      varEelec += eleExc[j] * eleExc[j] * exp(-beta*eleExc[j]);
    }
    Eelec /= qelec;
    varEelec = varEelec / qelec - Eelec*Eelec;

    PrtnFn   *= qtot*qelec;
    IntrlEne += ene + Eelec;
    varEne   += var + varEelec;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int ClassicalCoupledRotors::NoDegOfFreedom(gDensityOfStates* gdos) {
    unsigned int nDOF(m_bondIDs.size() + 3); // The +3 is for the external rotors.
    return nDOF;
  }

	// Calculates the potential energy associated with torsions.
	double ClassicalCoupledRotors::TorsionalPotential(const vector<double>& angles) {

		// SHR, 15/Oct/2017: At present the potential is modelled as the sum of separate 
		// torsion terms, one for each angle. However, a better potential would take 
		// account of non-bond terms.

		double ptnl(0.0);

		for (size_t i(0); i < angles.size(); ++i) {
			double theta = angles[i] ;
			const vector<double>& cosCoeff = m_potentialCosCoeffs[i];
			for (size_t k(0); k < cosCoeff.size(); ++k) {
				ptnl += cosCoeff[k] * cos(double(k)*theta);
			}
		}

		return ptnl;
	}

}//namespace
