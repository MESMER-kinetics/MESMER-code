
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
#include "../gDensityOfStates.h"
#include "../gStructure.h"
#include "HinderedRotorUtils.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class HinderedRotorQM1D : public HinderedRotorUtils
  {
  public:
    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol,  const MesmerEnv& env);

    // Provide a function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) {return 1 ; } ;

    // Provide a function to calculate the zero point energy of a molecule.
    virtual double ZeroPointEnergy(gDensityOfStates* gdos) {return m_ZPE ; } ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorQM1D(const char* id) : 
    HinderedRotorUtils(id),
      m_periodicity(1),
      m_kineticCosCoeff(),
      m_kineticSinCoeff(),
      m_energyLevels(),
      m_ZPE(0.0),
      m_plotStates(false),
      m_writeStates(false)
    { }

    virtual ~HinderedRotorQM1D() {}
    virtual HinderedRotorQM1D* Clone() { return new HinderedRotorQM1D(*this); }

  private:

    // Shift potential to origin.
    void ShiftPotential() ;

    // Provide data for plotting states against potential.
    void outputPlotData() const ;

    // Print the hindered rotor states.
    void outputStateData() const ;

    // Determine the Fourier components of the reciprocal moment of inertia.
    void RMoIFourierCoeffs() ;

    int    m_periodicity;

    vector<double> m_kineticCosCoeff ;   // The cosine coefficients of the internal rotor inertia term.
    vector<double> m_kineticSinCoeff ;   // The sine coefficients of the internal rotor inertia term.

    vector<double> m_energyLevels ;	     // The energies of the hindered rotor states.
    double m_ZPE ;                       // Zero point energy. 

    bool m_plotStates ;                  // If true output data for plotting. 
    bool m_writeStates ;                 // If true energy levels written to output. 
    PersistPtr m_ppConfigData;
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
    if(!bondID)
      bondID = ppDOSC->XmlReadValue("me:bondRef",optional);
    if (!bondID || *bondID=='\0')
    {
      cerr << "No <bondRef> specified for the hindered rotating bond" <<endl;
      return false;
    }

    // Save rotatable bond ID for calculation of GRIT.
    gs.addRotBondID(string(bondID)) ;

    pair<string,string> bondats = gs.GetAtomsOfBond(bondID);
    if(bondats.first.empty())
    {
      cerr << "Unknown bond reference " << bondID << endl;
      return false;
    }
    set_BondID(bondID) ;
    cinfo << "Hindered rotor " << get_BondID() ;  

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
    m_kineticSinCoeff.clear() ;
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
        coefficients.push_back(coefficient) ;
      }

      // As coefficients can be supplied in any order, they are sorted here.
      m_kineticCosCoeff.resize(++maxIndex) ;
      for (size_t i(0) ; i < coefficients.size() ; i++ ) {
        m_kineticCosCoeff[indicies[i]] = coefficients[i] ;
      }

      if (units == "amuA^2") {
        RMoIFourierCoeffs() ;
      }

    } else {

	  // No inertia terms red in so let us see if they are to be calculated.
	  // If so, set the appropriate flags.

      pp = ppDOSC->XmlMoveTo("me:CalculateInternalRotorInertia") ;
      if (pp) {
        set_CalIntrlIrt(true) ;
        double phase(0.0) ;
        phase = pp->XmlReadDouble("phaseDifference", optional);
        set_Phase(phase) ;
      } else {

        // No information supplied, so simply calculate reduced moment of inertia
		// for current configuration.

        double reducedMoI(gs.reducedMomentInertia(bondats)) ;  // Units a.u.*Angstrom*Angstrom.
        m_kineticCosCoeff.push_back(conMntInt2RotCnt/reducedMoI) ;
      }
    }

    // Read in potential information.

    vector<double>&  potentialCosCoeff = get_PotentialCosCoeff() ;
    vector<double>&  potentialSinCoeff = get_PotentialSinCoeff() ;

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
        potentialCosCoeff.resize(++maxIndex) ;
        potentialSinCoeff.resize(++maxIndex) ;
        for (size_t i(0) ; i < coefficients.size() ; i++ ) {
          potentialCosCoeff[indicies[i]] = coefficients[i] ;
          potentialSinCoeff[i]           = 0.0 ;
        }

        // Shift potential so that lowest minimum is at zero.

        ShiftPotential() ;

      } else if (format == "numerical") {

        // Numerical potential.

        vector<double> potential ;
        vector<double> angle ;
        set_Expansion(pp->XmlReadInteger("expansionSize",optional));

        // Check if sine terms are to be used.
        set_UseSinTerms(pp->XmlReadBoolean("useSineTerms") || pp->XmlReadBoolean("UseSineTerms"));

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

        PotentialFourierCoeffs(angle, potential) ;

      } else {

        // Unknown format.

        cinfo << "Unknown hindering potential format for " << bondID << ", assuming free rotor." <<endl;

        potentialCosCoeff.push_back(0.0) ;

      }

    } else {

      // Default : free rotor.

      cinfo << "No potential defined for " << bondID << ", assuming free rotor." <<endl;

      potentialCosCoeff.push_back(0.0) ;

    }

    // Check if there is a Hessian and knock out the frequency
    // associated with this internal rotation.

    if (gdos->hasHessian()) {
      // The following vector, "mode", will be used to hold the internal rotation 
      // mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010). 
      vector<double> mode(3*gs.NumAtoms(), 0.0);
      gs.internalRotationVector(get_BondID(), mode) ;
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

    // Check if configuration data are required.

    pp = ppDOSC->XmlMoveTo("me:ConfigurationalData") ;
    if (pp) {
      gs.set_Verbose(true) ;
      m_ppConfigData = pp;
    }

    return true;
  }

  //
  // Calculate quantum mechanical 1D rotor densities of states of an 
  // internal rotor and convolve them with the main density of states.
  //
  bool HinderedRotorQM1D::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
    const size_t MaximumCell = env.MaxCell ;

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // Retrieve the DOS vector without re-calculating.
      return false;

    vector<double> tmpCellDOS(cellDOS) ;

    const vector<double>&  potentialCosCoeff = get_PotentialCosCoeff() ;
    const vector<double>&  potentialSinCoeff = get_PotentialSinCoeff() ;
    const bool useSinTerms = get_UseSinTerms() ;
    const bool calIntrlIrt = get_CalIntrlIrt() ;

    // Calculate the (angle dependent) internal moment of inertia. This
    // has to be done at this point in order to capture any interaction
    // with other internal rotors.

    if (calIntrlIrt) {
      gStructure& gs = pDOS->getHost()->getStruc();
      size_t nAngle(36) ;
      vector<double> angle(nAngle,0.0), redInvMOI ;
      gs.reducedMomentInertiaAngular(get_BondID(), get_Phase(), angle, redInvMOI, m_ppConfigData) ;  // Units a.u.*Angstrom*Angstrom.
      FourierCosCoeffs(angle, redInvMOI, m_kineticCosCoeff, get_Expansion()) ;
      FourierSinCoeffs(angle, redInvMOI, m_kineticSinCoeff, get_Expansion()) ;
    }

    // Find maximum quantum No. for rotor. To ensure convergence basis functions 
    // that span a range twice that request are used with a minimum of 100000
    // wavenumbers which appears in the definition of root below. 

    double bint    = m_kineticCosCoeff[0] ;
    double root    = sqrt(double(max(2*MaximumCell,size_t(100000)))/bint) ;
    int kmax       = int(root + 1.0) ;
    size_t nstates = 2*kmax +1 ;

    // Check if sine terms are required and if so use the augmented matrix approach. See NR Sec. 11.4.

    size_t msize = (useSinTerms) ? 2*nstates : nstates ; 

    dMatrix hamiltonian(msize) ;

    vector<int> stateIndicies(nstates,0) ;

    // Add diagonal kinetic and potential terms first.

    hamiltonian[0][0] = potentialCosCoeff[0] ;
    for (int k(1), i(1); k <= kmax ; k++) {
      double energy = bint*double(k*k) + potentialCosCoeff[0] ;
      hamiltonian[i][i] = energy ;
      stateIndicies[i]  = -k ;
      i++ ;                         // Need to account for the two directions of rotation.
      hamiltonian[i][i] = energy ;
      stateIndicies[i]  = k ;
      i++ ;
    }

    // Add off-diagonal cosine potential terms.

    for (int n(1); n < int(potentialCosCoeff.size()) && n <= kmax ; n++) {
      double matrixElement = potentialCosCoeff[n]/2.0 ; 
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
            hamiltonian[i][j] += matrixElement*( 
			  (((k - jj + m) == 0) ? double(k*(k + m)) : 0.0) +
              (((k - jj - m) == 0) ? double(k*(k - m)) : 0.0) ) ;
          }
        }
      }
    }

    if (useSinTerms) {

      // Following the augmented matrix approach, first copy the cosine part to the lower right hand block.

      for (size_t i(0), ii(nstates); i < nstates; i++, ii++) {
        for (size_t j(0), jj(nstates); j < nstates; j++, jj++) {
          hamiltonian[ii][jj] = hamiltonian[i][j] ;
        }
      }

      // Now, construct the off-diagonal sine potential terms,
	  // placing result in the lower left off-diagoanl block.

      for (int n(1); n < int(potentialSinCoeff.size()) && n <= kmax ; n++) {
        double matrixElement = potentialSinCoeff[n]/2.0 ; 
        for (size_t i(0) ; i < nstates; i++) {
          for (size_t j(0) ; j < nstates; j++) {
            hamiltonian[nstates + i][j] += matrixElement*( 
              (((stateIndicies[j] - stateIndicies[i] - n) == 0) ? 1.0 : 0.0)
            - (((stateIndicies[j] - stateIndicies[i] + n) == 0) ? 1.0 : 0.0) ) ;
          }
        }
      }

      // Add off-diagonal sine kinetic terms.

      if (m_kineticSinCoeff.size() > 1) {
		for (int m(1); m < int(m_kineticSinCoeff.size()) && m <= kmax ; m++) {
		  double matrixElement = -m_kineticSinCoeff[m]/2.0 ; 
		  for (size_t i(0) ; i < nstates; i++) {
			int k = stateIndicies[i] ;
			for (size_t j(0) ; j < nstates; j++) {
			  int jj = stateIndicies[j] ;
			  hamiltonian[nstates + i][j] += matrixElement*( 
				(((k - jj + m) == 0) ? double(k*(k + m)) : 0.0) -
				(((k - jj - m) == 0) ? double(k*(k - m)) : 0.0) ) ;
			}
		  }
		}
      }

	  // Now, copy the negated off-diagonal sine potential and kinetic terms
	  // to the upper right off-diagonal block.

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

    if (useSinTerms) {
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
      size_t nr = nint((m_energyLevels[k] - m_ZPE)/env.CellSize) ;
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
  void HinderedRotorQM1D::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne)
  {
    double Qintrot(0.0), Eintrot(0.0), varEintrot(0.0) ;

    double zeroPointEnergy(m_energyLevels[0]) ; 
    for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
      double ene = m_energyLevels[i] - zeroPointEnergy ;
      Qintrot    += exp(-beta*ene) ;
      Eintrot    += ene*exp(-beta*ene) ;
      varEintrot += ene*ene*exp(-beta*ene) ;
    }
    Eintrot   /= Qintrot ;
    varEintrot = varEintrot/Qintrot - Eintrot*Eintrot ;
    Qintrot   /= double(m_periodicity) ;

    PrtnFn   *= Qintrot ;
    IntrlEne += Eintrot ;
    varEne   += varEintrot ;
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

    vector<double>&  potentialCosCoeff = get_PotentialCosCoeff() ;
    potentialCosCoeff[0] -= minPotential ;

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

    ctest << "\nHindered rotor states (cm-1).\n" << endl ;
    for (size_t i(0); i < m_energyLevels.size() ; i++) {
      ctest << formatFloat(m_energyLevels[i], 6, 15) << endl ;
    }
    ctest << endl ;
  }

  // Determine the Fourier components of the reciprocal moment of inertia.
  void HinderedRotorQM1D::RMoIFourierCoeffs() {

    size_t ndata = m_kineticCosCoeff.size() ;

    // Calculate the Moment of inertia at a set of points.

    const size_t nAngle(360) ;
    const double dAngle(2.0*M_PI/double(nAngle)) ;

    vector<double> MoI(nAngle, 0.0) ;
    vector<double> angle(nAngle, 0.0) ;
    for (size_t i(0); i < nAngle; ++i) {
      angle[i] = double(i)*dAngle ;
      double sum(0.0) ;
      for (size_t j(0); j < ndata; ++j) {
        double nTheta = double(j) * angle[i];
        sum += m_kineticCosCoeff[j] * cos(nTheta);
      }
      MoI[i] = conMntInt2RotCnt/sum ;
    }

    // Determine the cosine coefficients.

    for (size_t j(0); j < ndata; ++j) {
      double sum(0.0) ;
      for (size_t i(0); i < nAngle; ++i) {
        double nTheta = double(j) * angle[i];
        sum += MoI[i] * cos(nTheta);
      }
      m_kineticCosCoeff[j] = 2.0*sum/double(nAngle) ;
    }
    m_kineticCosCoeff[0] /= 2.0 ;
  }

}//namespace
