
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
#include "HinderedRotorUtils.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class HinderedRotorQM1D : public DensityOfStatesCalculator, HinderedRotorUtils
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

    // Provide a function to calculate the zero point energy of a molecule.
    virtual double ZeroPointEnergy(gDensityOfStates* gdos) {return m_ZPE ; } ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorQM1D(const char* id) : DensityOfStatesCalculator(),
      HinderedRotorUtils(),
      m_id(id),
      m_bondID(),
      m_periodicity(1),
      m_kineticCosCoeff(),
      m_energyLevels(),
      m_ZPE(0.0),
      m_plotStates(false),
      m_writeStates(false)
    { Register(); }

    virtual ~HinderedRotorQM1D() {}
    virtual const char* getID()  { return m_id; }
    virtual HinderedRotorQM1D* Clone() { return new HinderedRotorQM1D(*this); }

  private:
    const char* m_id;

    // Shift potential to origin.
    void ShiftPotential() ;

    // Provide data for plotting states against potential.
    void outputPlotData() const ;

    // Print the hindered rotor states.
    void outputStateData() const ;

    std::string m_bondID;

    int    m_periodicity;

    vector<double> m_kineticCosCoeff ;   // The cosine coefficients of the internal rotor inertia term.

    vector<double> m_energyLevels ;	     // The energies of the hindered rotor states.
    double m_ZPE ;                       // Zero point energy. 

    bool m_plotStates ;                  // If true output data for plotting. 
    bool m_writeStates ;                 // If true energy levels written to output. 

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

    // The following vector, "mode", will be used to hold the internal rotation 
    // mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010). 
    vector<double> mode(3*gs.NumAtoms(), 0.0);

    const char* bondID = ppDOSC->XmlReadValue("bondRef",optional);
    if(!bondID)
      bondID = ppDOSC->XmlReadValue("me:bondRef",optional);
    if (!bondID || *bondID=='\0')
    {
      cerr << "No <bondRef> specified for the hindered rotating bond" <<endl;
      return false;
    }

    pair<string,string> bondats = gs.GetAtomsOfBond(bondID);
    if(bondats.first.empty())
    {
      cerr << "Unknown bond reference " << bondID << endl;
      return false;
    }
    m_bondID = bondID;
    cinfo << "Hindered rotor " << m_bondID ;  

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
    cinfo << endl;

    // Determine reduced moment of inertia function. These can either 
    // be read in or calculated if the coordinates are available.

    m_kineticCosCoeff.clear() ;
    PersistPtr pp = ppDOSC->XmlMoveTo("me:InternalRotorInertia") ;

    if (pp) {
      const char* p = pp->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";

      vector<int> indicies ;
      vector<double> coefficients ;
      int maxIndex(0) ;
      while(pp = pp->XmlMoveTo("me:InertiaPoint"))
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
      m_kineticCosCoeff.resize(++maxIndex) ;
      for (size_t i(0) ; i < coefficients.size() ; i++ ) {
        m_kineticCosCoeff[indicies[i]] = coefficients[i] ;
      }
    } else {

      // Calculate reduced moment of inertia.

      double reducedMomentInertia(reducedMomentInertia(gs, bondats, mode)) ;  //units a.u.*Angstrom*Angstrom
      m_kineticCosCoeff.push_back(conMntInt2RotCnt/reducedMomentInertia) ;
    }

    // Read in potential information.

    m_periodicity = max(m_periodicity, ppDOSC->XmlReadInteger("me:periodicity",optional));

    pp = ppDOSC->XmlMoveTo("me:HinderedRotorPotential") ;

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

        // Shift potential so that lowest minimum is at zero.

        ShiftPotential() ;

      } else if (format == "numerical") {

        // Numerical potential.

        vector<double> potential ;
        vector<double> angle ;
        m_expansion = pp->XmlReadInteger("expansionSize",optional);

        // Check if sine terms are to be used.

        const char *pUseSineTerms(pp->XmlReadValue("useSineTerms",optional)) ;
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
      if(!gdos->projectMode(mode)) {
        cerr << "Failed to project out internal rotation." <<endl;
        return false;
      }
    }

    // Check if is data for plotting are required.

    pp = ppDOSC->XmlMoveTo("me:PlotStates") ;
    if (pp) {
      m_plotStates = true ;
    }

    // Check if is energy level values are to be written.

    pp = ppDOSC->XmlMoveTo("me:WriteStates") ;
    if (pp) {
      m_writeStates = true ;
    }

    return true;
  }

  //
  // Calculate quantum mechanical 1D rotor densities of states of an 
  // internal rotor and convolve them with the main density of states.
  //
  bool HinderedRotorQM1D::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    vector<double> tmpCellDOS(cellDOS) ;

    // Find maximum quantum No. for rotor.

    double bint    = m_kineticCosCoeff[0] ;
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

    // Add off-diagonal cosine kinetic terms.

    if (m_kineticCosCoeff.size() > 1) {
      for (int m(1); m < int(m_kineticCosCoeff.size()) ; m++) {
        double matrixElement = m_kineticCosCoeff[m]/2.0 ; 
        for (size_t i(0) ; i < nstates; i++) {
          int k = stateIndicies[i] ;
          for (size_t j(0) ; j < nstates; j++) {
            int jj = stateIndicies[j] ;
            hamiltonian[i][j] += matrixElement*( (((k - jj + m) == 0) ? double(k*(k + m)) : 0.0) +
                                                 (((k - jj - m) == 0) ? double(k*(k - m)) : 0.0) ) ;
          }
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

      // Add off-diagonal sine kinetic terms.

      if (m_kineticCosCoeff.size() > 1) {
        // Hmmm... Do these occur? Are they large? 
      }
    }

    // Now diagonalize hamiltonian matrix to determine energy levels.

    vector<double> eigenvalues(msize,0.0) ;

    hamiltonian.diagonalize(&eigenvalues[0]);

    // Save energy levels for partition function calculations.

    if (m_useSinTerms) {
      m_energyLevels.clear() ;
      for (size_t j(0) ; j < nstates ; j++) {
        m_energyLevels.push_back(eigenvalues[2*j]) ;
      }
    } else {
      m_energyLevels = eigenvalues ;
    }


    // Shift eigenvalues by the zero point energy and convolve with the 
    // density of states for the other degrees of freedom.

    m_ZPE = m_energyLevels[0] ;
    for (size_t k(1) ; k < nstates ; k++ ) {
      size_t nr = int(m_energyLevels[k] - m_ZPE) ;
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

    // If required, created graphical date.
    if (m_writeStates) 
      outputStateData() ;

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

  // Shift potential to origin.
  void HinderedRotorQM1D::ShiftPotential() {

    // A coarse search for minima is done over intervals of 1 degree a minima 
    // being located when the gradient changes sign from -ve to +ve. Then a 
    // Newton-Raphson type iteration is applied to get a better estimate of
    // the location of the minimum. Finally, the potential is shifted.

    double minPotential(1.e10) ; 
    double anga(0.0), grda(0.0) ; 
    size_t nDegrees(360) ;
    for (size_t i(0); i <= nDegrees; ++i) {
      double angb = double(i) * M_PI/180. ;
      double grdb = CalculateGradient(angb) ;
      if (grdb < 0.0) {
        anga = angb ;
        grda = grdb ;
      } else if (grda < 0.0) {
        double angc(0.0), grdc(1.0), tol( 1.e-05) ;
        int n(0) ;
        while (fabs(grdc) > tol && n < 10) {
          angc = (grdb*anga - grda*angb)/(grdb - grda) ;
          grdc = CalculateGradient(angc) ;
          if (grdc > 0.0) {
            anga = angc ;
            grda = grdc ;
          } else {
            angb = angc ;
            grdb = grdc ;
          }
          n++ ;
        }
        double potential = CalculatePotential(angc) ;
        minPotential = min(minPotential, potential) ;
      }
    }

    m_potentialCosCoeff[0] -= minPotential ;

    return ;
  }

  // Provide data for plotting states against potential.
  void HinderedRotorQM1D::outputPlotData() const {

    ctest << endl << "Hindered rotor data for plotting." << endl << endl ;
    int npoints(500) ;
    double dAngle = M_PI/double(npoints) ;
    for (int i(-npoints); i < npoints; ++i) {     
      double angle = double(i)*dAngle ;
      double potential = CalculatePotential(angle) ;
      ctest << formatFloat(angle, 6, 15) << ", "<< formatFloat(potential, 6, 15) << endl ;
    }
    ctest << endl ;

    for (size_t i(0); i < m_energyLevels.size() ; i++) {
      ctest << formatFloat(-M_PI, 6, 15) << ", "<< formatFloat(m_energyLevels[i], 6, 15) << endl ;
      ctest << formatFloat( M_PI, 6, 15) << ", "<< formatFloat(m_energyLevels[i], 6, 15) << endl ;    
      ctest << endl ;
    }

  }

  // Print the hindered rotor states.
  void HinderedRotorQM1D::outputStateData() const {

    ctest << endl << "Hindered rotor states (cm-1)." << endl << endl ;
    for (size_t i(0); i < m_energyLevels.size() ; i++) {
      ctest << formatFloat(m_energyLevels[i], 6, 15) << endl ;
    }
    ctest << endl ;
  }

}//namespace
