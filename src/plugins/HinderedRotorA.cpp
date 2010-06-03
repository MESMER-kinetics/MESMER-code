#include "../Molecule.h"
#include "HinderedRotorA.h"

using namespace std;
namespace mesmer
{
  //-------------------------------------------------------------
  //Global instance, defining its id
  HinderedRotorA theHinderedRotorA("HinderedRotorA");
  //-------------------------------------------------------------

  // Convertion factor need to obtain rotational constant. Moment of Inertia must be in amu Ang^2.
  static const double conMntInt2RotCnt = 16.85917 ;

  using OpenBabel::vector3;
  //Read data from XML and store in this instance.
  bool HinderedRotorA::ReadParameters(Molecule* pMol, PersistPtr ppDOSC)
  {
    gStructure& gs = pMol->getStruc();
    if(!gs.ReadStructure())
    {
      cerr << "A complete set of atom coordinates are required for hindered rotor calculations" <<endl;
      return false;
    }

    const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);
    if(bondID)
    {
      pair<string,string> bondats = gs.GetAtomsOfBond(bondID);
      if(bondats.first.empty())
      {
        cerr << "Unknown bond reference " << bondID << endl;
        return false;
      }
      m_bondID = bondID;
      cinfo << "Hindered rotor " << m_bondID << endl;

      vector3 coords1 = gs.GetAtomCoords(bondats.first);
      vector3 coords2 = gs.GetAtomCoords(bondats.second);

      //calc Moment of inertia about bond axis of atoms on one side of bond...
      vector<string> atomset;
      gs.GetAttachedAtoms(atomset, bondats.first, bondats.second);
      double mm1 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

      //...and the other side of the bond
      atomset.clear();
      gs.GetAttachedAtoms(atomset, bondats.second, bondats.first);
      double mm2 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

      /*
      Is the reduced moment of inertia need about the bond axis or, separately for the set of
      atoms on each side of the bond, about a parallel axis through their centre of mass?
      See:
      http://www.ccl.net/chemistry/resources/messages/2001/03/21.005-dir/index.html
      http://www.ccl.net/chemistry/resources/messages/2001/03/31.002-dir/index.html
      The bond axis is used here.
      */
      m_reducedMomentInertia = mm1 * mm2 / ( mm1 + mm2 );//units a.u.*Angstrom*Angstrom
    }

    m_barrier  = ppDOSC->XmlReadDouble("me:barrierZPE",optional);
    PersistPtr pp = ppDOSC->XmlMoveTo("me:barrierZPE");
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "kJ/mol";
    m_barrier = getConvertedEnergy(units, m_barrier);

    m_periodicity = ppDOSC->XmlReadInteger("me:periodicity",optional);

    m_vibFreq = ppDOSC->XmlReadDouble("me:vibFreq",optional);
    pp = ppDOSC->XmlMoveTo("me:vibFreq");
    p = pp->XmlReadValue("units", optional);
    units = p ? p : "cm-1";

    // other inputs...
    return true;
  }

  //
  // Calculate quantum mechanical 1D rotor densities of states of a free 
  // rotor and convolve them with the main density of states.
  //
  bool HinderedRotorA::countCellDOS(gDensityOfStates* pDOS, int MaximumCell)
  {

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    vector<double> tmpCellDOS(cellDOS) ;

    // Find maximum quantum No. for rotor.

    double bint = conMntInt2RotCnt/m_reducedMomentInertia ;
    double root = sqrt(double(MaximumCell)/bint) ;
    int kmax    = int(root + 1.0) ;
	int nstates = 2*kmax +1 ;

	dMatrix hamiltonian(nstates) ;

	for (int k(1), i(1); k <= kmax ; k++) {
		double energy = bint*double(k*k);
		hamiltonian[i][i] = energy ;
		i++ ;
		hamiltonian[i][i] = energy ;
		i++ ;
	}

	// Now call eigenvalue solvers.

	vector<double> eigenvalues(nstates,0.0) ;

	hamiltonian.diagonalize(&eigenvalues[0]);

    // Now convolve with the density of states for the other degrees of freedom.

    for (int k(1) ; k < nstates ; k++ ) {
      int nr = int(eigenvalues[k]) ;
      if (nr < MaximumCell) {
        for (int i(0) ; i < MaximumCell - nr ; i++ ) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i] ;
        }
      }
    }

    // Apply symmetry number.

    for (int i(0) ; i < MaximumCell ; i++ ) {
      tmpCellDOS[i] /= double(m_periodicity) ;
    }

    // Replace existing density of states.   

    pDOS->setCellDensityOfStates(tmpCellDOS) ;

    return true;

  }

  //
  // Provide a function to calculate contribution to canonical partition function.
  // (Mostly for testing purposes.)
  //
  double HinderedRotorA::canPrtnFnCntrb(const double beta)
  {
    return sqrt(M_PI*m_reducedMomentInertia/conMntInt2RotCnt/beta)/double(m_periodicity) ;
  }

}//namespace
