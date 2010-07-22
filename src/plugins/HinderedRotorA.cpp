
#include "HinderedRotorA.h"

#include "../Molecule.h"
#include "../Constants.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  //-------------------------------------------------------------
  //Global instance, defining its id
  HinderedRotorA theHinderedRotorA("HinderedRotorA");
  //-------------------------------------------------------------

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

	// Read in potential information.

	m_periodicity = max(m_periodicity, ppDOSC->XmlReadInteger("me:periodicity",optional));

	PersistPtr pp = ppDOSC->XmlMoveTo("me:HinderedRotorPotential") ;

	if (pp) {

	  const char* p = pp->XmlReadValue("format", true);
	  string format(p) ;

	  p = pp->XmlReadValue("units", optional);
	  string units = p ? p : "kJ/mol";

	  if (format == "analytical") {

		// Analytical potential.

		vector<int> indicies ;
		vector<double> coefficients ;
		int maxIndex(0) ;
		while(pp = pp->XmlMoveTo("me:PotentialPoint"))
		{
		  int index = pp->XmlReadInteger("index", optional);
		  indicies.push_back(index) ;
		  maxIndex = max(maxIndex,index) ;

		  double coefficient = pp->XmlReadDouble("coefficient", optional);
		  coefficient = getConvertedEnergy(units, coefficient);
		  coefficients.push_back(coefficient) ;
		}

		// As coefficients can be supplied in any order, they are sorted here.
		m_potentialCosCoeff.resize(++maxIndex) ;
		for (size_t i(0) ; i < coefficients.size() ; i++ ) {
		  m_potentialCosCoeff[indicies[i]] = coefficients[i] ;
		}

	  } else if (format == "numerical") {

		// Numerical potential.

		vector<double> potential ;
		vector<double> angle ;
		m_expansion = pp->XmlReadInteger("expansionSize",optional);
		while(pp = pp->XmlMoveTo("me:PotentialPoint"))
		{
		  double anglePoint = pp->XmlReadDouble("angle", optional);
		  angle.push_back(anglePoint) ;

		  double potentialPoint = pp->XmlReadDouble("potential", optional);
		  potentialPoint = getConvertedEnergy(units, potentialPoint);
		  potential.push_back(potentialPoint) ;
		}

		CosineFourierCoeffs(angle, potential) ;

	  } else {

		// Unknown format.

		cinfo << "Unknown hindering potential format for " << bondID << ", assuming free rotor." <<endl;

		m_potentialCosCoeff.push_back(0.0) ;
		
	  }

	} else {

	  // Default : free rotor.

	  cinfo << "No potential defined for " << bondID << ", assuming free rotor." <<endl;

	  m_potentialCosCoeff.push_back(0.0) ;

	}

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

	// Add diagonal kinetic and potential terms first.

	hamiltonian[0][0] = m_potentialCosCoeff[0] ;
	for (int k(1), i(1); k <= kmax ; k++) {
	  double energy = bint*double(k*k) + m_potentialCosCoeff[0] ;
	  hamiltonian[i][i] = energy ;
	  i++ ;                         // Need to account for the two directions of rotation.
	  hamiltonian[i][i] = energy ;
	  i++ ;
	}

	// Add off-diagonal potential terms.

	for (int n(1); n < int(m_potentialCosCoeff.size()) && n <= kmax ; n++) {
	  int idx = 2*n - 1 ;
	  double matrixElement = m_potentialCosCoeff[n]/2.0 ; 
	  hamiltonian[idx][0] = hamiltonian[0][idx] = matrixElement ;
	  idx++ ;
	  hamiltonian[idx][0] = hamiltonian[0][idx] = matrixElement ;
	  idx++ ;
	  for (int k(1) ; idx < nstates ; idx++, k++) {
		hamiltonian[idx][k] = hamiltonian[k][idx] = matrixElement ;
	  }
	}

	// Now diagonalize hamiltonian matrix to determine energy levels.

	vector<double> eigenvalues(nstates,0.0) ;

	hamiltonian.diagonalize(&eigenvalues[0]);

	// Shift eigenvalues by the zero point energy and convolve with the 
	// density of states for the other degrees of freedom.

	double zeroPointEnergy(eigenvalues[0]) ;
	for (int k(1) ; k < nstates ; k++ ) {
	  int nr = int(eigenvalues[k] - zeroPointEnergy) ;
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
	double Qintrot = sqrt(M_PI*m_reducedMomentInertia/conMntInt2RotCnt/beta)/double(m_periodicity) ;

    Qintrot = exp(-beta*m_potentialCosCoeff[0]) ;
	for (size_t n(1); n < m_potentialCosCoeff.size() ; n++) {
      Qintrot *= ModifiedBessalFuncion(beta*m_potentialCosCoeff[n]) ;
	}

	return Qintrot ;
  }

  //
  // Calculation of the modifed Bessel function, Io(x), for real x.
  //
  double HinderedRotorA::ModifiedBessalFuncion(const double x) const 
  {
	static double p1(1.0),          p2(3.5156229),     p3(3.0899424) ;
	static double p4(1.2067492),    p5(0.2659732),     p6(0.360768e-1) ;
	static double p7(0.45813e-2) ;

	static double q1(0.39894228),   q2(0.1328592e-1),  q3(0.225319e-2) ;
	static double q4(-0.157565e-2), q5(0.916281e-2),   q6(-0.2057706e-1) ;
	static double q7(0.2635537e-1), q8(-0.1647633e-1), q9(0.392377e-2) ;

	//  Bessel index

	if (abs(x) < 3.75) {

	  double y1 = (x/3.75) ;
	  y1 *= y1 ;

	  return (p1+y1* (p2+y1* (p3+y1* (p4+y1* (p5+y1* (p6+ y1*p7)))))) ;

	} else {

	  double ax = abs(x) ;
	  double y1 = 3.75/ax ;
	  double b1 = (exp(ax)/sqrt(ax)) ;
	  double fd = (q6+y1* (q7+y1* (q8+y1*q9))) ;
	  double b2 = (q1+y1* (q2+y1* (q3+y1* (q4+y1* (q5+y1*fd))))) ;

	  return b1*b2 ;

	}
  }

  //
  // Calculate cosine coefficients from potential data points.
  //
  void HinderedRotorA::CosineFourierCoeffs(vector<double> &angle, vector<double> &potential)
  {
	size_t ndata = potential.size() ;

	// Locate the potential minimum and shift to that minimum.

	double vmin(potential[0]), amin(angle[0]) ;
	for (size_t i(1); i < ndata; ++i) {
	  if (potential[i] < vmin){
		vmin = potential[i] ;
		amin = angle[i] ;
	  }
	}

	for (size_t i(0); i < ndata; ++i) {
	  potential[i] -= vmin ;
	  angle[i]     -= amin ;
	  angle[i]     *= M_PI/180. ;
	}

	// Now determine the cosine coefficients.

	for(size_t k(0); k < m_expansion; ++k) {
	  double sum(0.0) ;
	  for(size_t i(0); i < ndata; ++i) {
		double nTheta = double(k) * angle[i];
		sum += potential[i] * cos(nTheta);
	  }
	  m_potentialCosCoeff.push_back(2.0*sum/double(ndata)) ;
	}
	m_potentialCosCoeff[0] /= 2.0 ;

	// Test potential

	cinfo << endl ;
	cinfo << "          Angle      Potential         Series" << endl ;
	cinfo << endl ;
	for (size_t i(0); i < ndata; ++i) {
	  double sum(0.0) ;
	  for(size_t k(0); k < m_expansion; ++k) {
		double nTheta = double(k) * angle[i];
		sum += m_potentialCosCoeff[k] * cos(nTheta);
	  }
	  cinfo << formatFloat(angle[i], 6, 15) << formatFloat(potential[i], 6, 15) << formatFloat(sum, 6, 15) << endl ;
	}
	cinfo << endl ;

	return ;

  }

}//namespace
