
//-------------------------------------------------------------------------------------------
//
// HinderedRotorQM1D.cpp
//
// This file contains the implementation of the methods for calculating and testing the 
// density of states of a one dimensional quantum mechanical hindered rotor.
//
// The Hamiltonian is represented in a basis of one-dimensional free rotor functions and 
// then diagonalized.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"
#include "../Molecule.h"
#include "../Constants.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class HinderedRotorQM1D : public DensityOfStatesCalculator
  {
  public:
    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, int MaximumCell);

    // Provide a function to calculate contribution to canonical partition function.
    // (Mostly for testing purposes.)
    virtual double canPrtnFnCntrb(const double beta) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorQM1D(const std::string& id) : DensityOfStatesCalculator(id, true),
      m_bondID(),
      m_reducedMomentInertia(0.0),
      m_periodicity(1),
      m_potentialCosCoeff(),
      m_expansion(4),
	  m_energyLevels() {}

    virtual ~HinderedRotorQM1D() {}
    virtual HinderedRotorQM1D* Clone() { return new HinderedRotorQM1D(*this); }

  private:

    // Calculate cosine coefficients from potential data points.
    void CosineFourierCoeffs(vector<double> &angle, vector<double> &potential) ;

    std::string m_bondID;

    double m_reducedMomentInertia;
    int    m_periodicity;

	vector<double> m_potentialCosCoeff ; // The cosine coefficients of the hindered rotor potential.

    size_t m_expansion ;                 // Number of coefficients in the cosine expansion.

	vector<double> m_energyLevels ;	     // The energies of the hindered rotor states.

  } ;

  //-------------------------------------------------------------
  //Global instance, defining its id
  HinderedRotorQM1D theHinderedRotorQM1D("HinderedRotorQM1D");
  //-------------------------------------------------------------

  using OpenBabel::vector3;
  //Read data from XML and store in this instance.
  bool HinderedRotorQM1D::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
  {
    gStructure& gs = gdos->getHost()->getStruc();
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
      cinfo << "Hindered rotor " << m_bondID;  

      //Remove the vibrational frequency that this hindered rotation replaces
      const char* vibFreq = ppDOSC->XmlReadValue("me:replaceVibFreq",optional);
      if(vibFreq)
      {
        if(!gdos->removeVibFreq(atof(vibFreq)))
        {
          cerr << "Cannot find vibrational frequency " << vibFreq << " to replace it with hindered rotor" <<endl;
          return false;
        }
        cinfo << " replacing vib freq " << vibFreq;      
      }
      cinfo << endl;

      vector3 coords1 = gs.GetAtomCoords(bondats.first);
      vector3 coords2 = gs.GetAtomCoords(bondats.second);

      //Calc moment of inertia about bond axis of atoms on one side of bond...
      vector<string> atomset;
      atomset.push_back(bondats.second); //will not look beyond this atom on the other side of the bond
      gs.GetAttachedAtoms(atomset, bondats.first);
      atomset.erase(atomset.begin()); //the other side of the bond is not in this set
      double mm1 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

      //...and the other side of the bond
      atomset.clear();
      atomset.push_back(bondats.first);
      gs.GetAttachedAtoms(atomset, bondats.second);
      atomset.erase(atomset.begin());
      double mm2 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);

      /*
      Is the reduced moment of inertia needed about the bond axis or, separately for the set of
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
          if(IsNan(coefficient))
            coefficient = 0.0;
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
          if(IsNan(anglePoint))
            anglePoint = 0.0;
          angle.push_back(anglePoint) ;

          double potentialPoint = pp->XmlReadDouble("potential", optional);
          if(IsNan(potentialPoint))
            potentialPoint = 0.0;
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
  bool HinderedRotorQM1D::countCellDOS(gDensityOfStates* pDOS, int MaximumCell)
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

    // Save energy levels for partition function calculations.

    m_energyLevels = eigenvalues ;

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
  double HinderedRotorQM1D::canPrtnFnCntrb(const double beta)
  {
    double Qintrot(0.0) ;

    double zeroPointEnergy(m_energyLevels[0]) ; 
    for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
      Qintrot += exp(-beta*(m_energyLevels[i] - zeroPointEnergy)) ;
    }

    return Qintrot/double(m_periodicity) ;
  }

  //
  // Calculate cosine coefficients from potential data points.
  //
  void HinderedRotorQM1D::CosineFourierCoeffs(vector<double> &angle, vector<double> &potential)
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
