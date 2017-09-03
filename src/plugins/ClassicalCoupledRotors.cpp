
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

using namespace std;
namespace mesmer
{
  class ClassicalCoupledRotors : public DensityOfStatesCalculator
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

    // Constructor which registers with the list of DensityOfStatesCalculators in TopPlugin
    // This class calculates a complete DOS: it is not an extra class. 
    ClassicalCoupledRotors(const char* id) : m_id(id), m_Sym(1.0), m_OpticalSym(1.0), m_SpinMultiplicity(1) { Register(); }

    virtual ~ClassicalCoupledRotors() {}
    virtual const char* getID() { return m_id; }

    // Included only in a subset of DensityOfStatesCalculators.
    // Otherwise the baseclass function returns false.
    virtual bool includesRotations() { return true; }

    virtual ClassicalCoupledRotors* Clone() { return new ClassicalCoupledRotors(*this); }

  private:

    const char* m_id;

    double m_Sym;              // Rotational Symmetry Number.
    double m_OpticalSym;       // Transition states may have 2 optical isomers, which
                               // is accounted for by an additionsl symmetry number. 
    size_t m_SpinMultiplicity; // Spin multiplicity.

		vector<string> m_bondIDs;  // The IDs of the binds to be used in the coupled model.

  };

  //************************************************************
  //Global instance, defining its id
  ClassicalCoupledRotors theClassicalCoupledRotors("ClassicalCoupledRotorss");
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

		// Read bond IDs
		const char* bondIDs = ppDOSC->XmlReadValue("me:bondRef", optional);
		if (!bondIDs || *bondIDs == '\0')
		{
			cerr << "No <bondRef> specified for the coupled hindered rotating bond" << endl;
			return false;
		}
		else {
			string idata(bondIDs), delim(" ");
			size_t prev = 0, pos = 0;
			do
			{
				pos = idata.find(delim, prev);
				if (pos == string::npos) pos = idata.length();
				string token = idata.substr(prev, pos - prev);
				if (!token.empty()) m_bondIDs.push_back(token);
				prev = pos + delim.length();
			} while (pos < idata.length() && prev < idata.length());
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

		// Main loop

		const size_t npnts(2000);
		long long seed(1);
		for (size_t i(0); i < npnts; ++i) {
			vector<double> rndmd(m_bondIDs.size(), 0.0);
			sobol.sobol(rndmd.size(), &seed, rndmd);

			// gs.reducedMomentInertiaAngular(get_BondID(), get_Phase(), angle, redInvMOI, m_ppConfigData);  // Units a.u.*Angstrom*Angstrom.


		}

    // Electronic excited states.
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
    vector<double> tmpCellDOS(cellDOS);
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

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    double qtot(1.0), ene(0.0), var(0.0);
    qtot *= double(m_SpinMultiplicity);

    switch (rotorType) {
    case NONLINEAR://3-D symmetric/asymmetric/spherical top
      qtot *= (sqrt(M_PI / (rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta, -1.5)) / m_Sym);
      ene = 1.5 / beta;
      var = 1.5 / (beta*beta);
      break;
    case LINEAR://2-D linear
      qtot /= (rotConst[0] * m_Sym*beta);
      ene = 1.0 / beta;
      var = 1.0 / (beta*beta);
      break;
    default:
      break; // Assume atom.
    }

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

    PrtnFn *= qtot*qelec;
    IntrlEne += ene + Eelec;
    varEne += var + varEelec;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int ClassicalCoupledRotors::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    unsigned int nDOF(3);
    return nDOF;
  }

}//namespace
