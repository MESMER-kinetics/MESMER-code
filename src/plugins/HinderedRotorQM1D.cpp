
//-------------------------------------------------------------------------------------------
//
// HinderedRotorQM1D.cpp
//
// This file contains the implementation of the methods for calculating and testing the 
// density of states of a one dimensional quantum mechanical hindered rotor.
//
// The Hamiltonian is represented in a basis of one-dimensional free rotor functions and 
// then diagonalized. (The implementation given here is closely related to that described
// by J.D. Lewis in J. Mol. Struct. p. 427, Vol. 12 (1972).)
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
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell);

    // Provide a function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) {return 1 ; } ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorQM1D(const char* id) : m_id(id),
      m_bondID(),
      m_reducedMomentInertia(0.0),
      m_periodicity(1),
      m_potentialCosCoeff(),
      m_potentialSinCoeff(),
      m_expansion(4),
      m_energyLevels(),
      m_plotStates(false),
      m_useSinTerms(false)  
    { Register(); }

    virtual ~HinderedRotorQM1D() {}
    virtual const char* getID()  { return m_id; }
    virtual HinderedRotorQM1D* Clone() { return new HinderedRotorQM1D(*this); }

  private:
    const char* m_id;

    // Calculate the Fourier coefficients from potential data points.
    void FourierCoeffs(vector<double> &angle, vector<double> &potential) ;

    // Provide data for plotting states against potential.
    void outputPlotData() ;

    std::string m_bondID;

    double m_reducedMomentInertia;
    int    m_periodicity;

    vector<double> m_potentialCosCoeff ; // The cosine coefficients of the hindered rotor potential.
    vector<double> m_potentialSinCoeff ; // The sine coefficients of the hindered rotor potential.

    size_t m_expansion ;                 // Number of coefficients in the cosine expansion.

    vector<double> m_energyLevels ;	     // The energies of the hindered rotor states.

    bool m_plotStates ;                  // If true output data for plotting. 

    bool m_useSinTerms ;                 // If true sine terms are used in the representation of the potential.
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

	// The following vector, "mode", will be used to hold the ther internal rotation 
	// mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010). 
    vector<double> mode(3*gs.NumAtoms(), 0.0);

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
      if (vibFreq)
      {
        if(!gdos->removeVibFreq(atof(vibFreq)))
        {
          cerr << "Cannot find vibrational frequency " << vibFreq << " to replace it with hindered rotor" <<endl;
          return false;
        }
        cinfo << " replacing vib freq " << vibFreq;      
      }
      cinfo << '\n';

      vector3 coords1 = gs.GetAtomCoords(bondats.first);
      vector3 coords2 = gs.GetAtomCoords(bondats.second);

      //Calc moment of inertia about bond axis of atoms on one side of bond...
      vector<string> atomset;
      atomset.push_back(bondats.second); //will not look beyond this atom on the other side of the bond
      gs.GetAttachedAtoms(atomset, bondats.first);
      atomset.erase(atomset.begin()); //the other side of the bond is not in this set
      double mm1 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);
      gs.CalcInternalRotVec(atomset, coords1, coords2, mode) ;

      //...and the other side of the bond
      atomset.clear();
      atomset.push_back(bondats.first);
      gs.GetAttachedAtoms(atomset, bondats.second);
      atomset.erase(atomset.begin());
      double mm2 = gs.CalcMomentAboutAxis(atomset, coords1, coords2);
      gs.CalcInternalRotVec(atomset, coords2, coords1, mode) ;

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
        m_potentialSinCoeff.resize(++maxIndex) ;
        for (size_t i(0) ; i < coefficients.size() ; i++ ) {
          m_potentialCosCoeff[indicies[i]] = coefficients[i] ;
          m_potentialSinCoeff[i]           = 0.0 ;
        }

      } else if (format == "numerical") {

        // Numerical potential.

        vector<double> potential ;
        vector<double> angle ;
        m_expansion = pp->XmlReadInteger("expansionSize",optional);

        // Check if sine terms are to be used.

        const char *pUseSineTerms(pp->XmlReadValue("UseSineTerms",optional)) ;
        if (pUseSineTerms && string(pUseSineTerms) == "yes") {
          m_useSinTerms = true ;
        }

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

        FourierCoeffs(angle, potential) ;

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

    // Check if there is a Hessian and knock out the frequency
    // associated with this internal rotation.

    if (gdos->hasHessian()) {
      bool status = gdos->projectMode(mode) ;
    }

    // Check if is data for plotting are required.

    pp = ppDOSC->XmlMoveTo("me:PlotStates") ;
    if (pp) {
      m_plotStates = true ;
    }

    return true;
  }

  //
  // Calculate quantum mechanical 1D rotor densities of states of a free 
  // rotor and convolve them with the main density of states.
  //
  bool HinderedRotorQM1D::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    vector<double> tmpCellDOS(cellDOS) ;

    // Find maximum quantum No. for rotor.

    double bint    = conMntInt2RotCnt/m_reducedMomentInertia ;
    double root    = sqrt(double(MaximumCell)/bint) ;
    int kmax       = int(root + 1.0) ;
    size_t nstates = 2*kmax +1 ;

    // Check if sine terms are required and if so use the augmented matrix approach. See NR Sec. 11.4.

    size_t msize = (m_useSinTerms) ? 2*nstates : nstates ; 

    dMatrix hamiltonian(msize) ;

    vector<int> stateIndicies(nstates,0) ;

    // Add diagonal kinetic and potential terms first.

    hamiltonian[0][0] = m_potentialCosCoeff[0] ;
    for (int k(1), i(1); k <= kmax ; k++) {
      double energy = bint*double(k*k) + m_potentialCosCoeff[0] ;
      hamiltonian[i][i] = energy ;
      stateIndicies[i]  = -k ;
      i++ ;                         // Need to account for the two directions of rotation.
      hamiltonian[i][i] = energy ;
      stateIndicies[i]  = k ;
      i++ ;
    }

    // Add off-diagonal cosine potential terms.

    for (int n(1); n < int(m_potentialCosCoeff.size()) && n <= kmax ; n++) {
      double matrixElement = m_potentialCosCoeff[n]/2.0 ; 
      for (size_t i(0) ; i < nstates; i++) {
        for (size_t j(0) ; j < nstates; j++) {
          hamiltonian[i][j] += matrixElement*(((abs(stateIndicies[j] - stateIndicies[i]) - n) == 0) ? 1.0 : 0.0)  ;
        }
      }
    }

    if (m_useSinTerms) {

      // Following the augmented matrix approach, first copy the cosine part to the lower right hand block.

      for (size_t i(0), ii(nstates); i < nstates; i++, ii++) {
        for (size_t j(0), jj(nstates); j < nstates; j++, jj++) {
          hamiltonian[ii][jj] = hamiltonian[i][j] ;
        }
      }

      // Now, construct the off-diagonal sine potential terms, placing result in the lower left off-diagoanl block.

      for (int n(1); n < int(m_potentialSinCoeff.size()) && n <= kmax ; n++) {
        double matrixElement = m_potentialSinCoeff[n]/2.0 ; 
        for (size_t i(0) ; i < nstates; i++) {
          for (size_t j(0) ; j < nstates; j++) {
            hamiltonian[nstates + i][j] += matrixElement*( 
              (((stateIndicies[j] - stateIndicies[i] - n) == 0) ? 1.0 : 0.0)
              -	(((stateIndicies[j] - stateIndicies[i] + n) == 0) ? 1.0 : 0.0) ) ;
          }
        }
      }

      // Now, copy the negated off-diagonal sine potential terms to the upper right off-diagonal block.

      for (size_t i(0), ii(nstates); i < nstates; i++, ii++) {
        for (size_t j(0), jj(nstates); j < nstates; j++, jj++) {
          hamiltonian[i][jj] = -hamiltonian[ii][j] ;
        }
      }

    }

    // Now diagonalize hamiltonian matrix to determine energy levels.

    vector<double> eigenvalues(msize,0.0) ;

    hamiltonian.diagonalize(&eigenvalues[0]);

    // Save energy levels for partition function calculations.

    if (m_useSinTerms) {
      for (size_t j(0) ; j < nstates; j++) {
        m_energyLevels.push_back(eigenvalues[2*j]) ;
      }
    } else {
      m_energyLevels = eigenvalues ;
    }


    // Shift eigenvalues by the zero point energy and convolve with the 
    // density of states for the other degrees of freedom.

    double zeroPointEnergy(m_energyLevels[0]) ;
    for (size_t k(1) ; k < nstates ; k++ ) {
      size_t nr = int(m_energyLevels[k] - zeroPointEnergy) ;
      if (nr < MaximumCell) {
        for (size_t i(0) ; i < MaximumCell - nr ; i++ ) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i] ;
        }
      }
    }

    // Apply symmetry number.

    for (size_t i(0) ; i < MaximumCell ; i++ ) {
      tmpCellDOS[i] /= double(m_periodicity) ;
    }

    // Replace existing density of states.   

    pDOS->setCellDensityOfStates(tmpCellDOS) ;

    // If required, created graphical date.
    if (m_plotStates) 
      outputPlotData() ;

    return true;

  }

  //
  // Provide a function to calculate contribution to canonical partition function.
  // (Mostly for testing purposes.)
  //
  double HinderedRotorQM1D::canPrtnFnCntrb(gDensityOfStates* gdos, double beta)
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
  void HinderedRotorQM1D::FourierCoeffs(vector<double> &angle, vector<double> &potential)
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

    // Determine the cosine coefficients.

    for(size_t k(0); k < m_expansion; ++k) {
      double sum(0.0) ;
      for(size_t i(0); i < ndata; ++i) {
        double nTheta = double(k) * angle[i];
        sum += potential[i] * cos(nTheta);
      }
      m_potentialCosCoeff.push_back(2.0*sum/double(ndata)) ;
    }
    m_potentialCosCoeff[0] /= 2.0 ;

    // Determine the sine coefficients.

    if (m_useSinTerms) {
      for(size_t k(0); k < m_expansion; ++k) {
        double sum(0.0) ;
        for(size_t i(0); i < ndata; ++i) {
          double nTheta = double(k) * angle[i];
          sum += potential[i] * sin(nTheta);
        }
        m_potentialSinCoeff.push_back(2.0*sum/double(ndata)) ;
      }
      m_potentialSinCoeff[0] = 0.0 ;
    } else {
      for(size_t k(0); k < m_expansion; ++k) {
        m_potentialSinCoeff.push_back(0.0) ;
      }
    }

    // Test potential
    cinfo << "          Angle         Potential          Series\n";
    for (size_t i(0); i < ndata; ++i) {
      double sum(0.0) ;
      for(size_t k(0); k < m_expansion; ++k) {
        double nTheta = double(k) * angle[i];
        sum += m_potentialCosCoeff[k] * cos(nTheta);
        sum += m_potentialSinCoeff[k] * sin(nTheta);
      }
      cinfo << formatFloat(angle[i], 6, 15) << ", " <<  formatFloat(potential[i], 6, 15) << ", " <<  formatFloat(sum, 6, 15) <<'\n' ;
    }
    cinfo << endl ;

    return ;

  }

  // Provide data for plotting states against potential.
  void HinderedRotorQM1D::outputPlotData() {

    cinfo << endl << "Hindered rotor data for plotting." << endl << endl ;
    int npoints(500) ;
    double dAngle = M_PI/double(npoints) ;
    for (int i(-npoints); i < npoints; ++i) {
      double sum(0.0) ;
      double angle = double(i)*dAngle ;
      for(size_t k(0); k < m_expansion; ++k) {
        double nTheta = double(k) * angle;
        sum += m_potentialCosCoeff[k] * cos(nTheta);
        sum += m_potentialSinCoeff[k] * sin(nTheta);
      }
      cinfo << formatFloat(angle, 6, 15) << ", "<< formatFloat(sum, 6, 15) << endl ;
    }
    cinfo << endl ;

    for (size_t i(0); i < m_energyLevels.size() ; i++) {
      cinfo << formatFloat(-M_PI, 6, 15) << ", "<< formatFloat(m_energyLevels[i], 6, 15) << endl ;
      cinfo << formatFloat( M_PI, 6, 15) << ", "<< formatFloat(m_energyLevels[i], 6, 15) << endl ;    
      cinfo << endl ;
    }

  }

}//namespace
