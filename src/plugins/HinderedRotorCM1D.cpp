
//-------------------------------------------------------------------------------------------
//
// HinderedRotorCM1D.cpp
//
// This file contains the implementation of the methods for calculating and testing the 
// density of states of a one dimensional classical mechanical hindered rotor.
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
  class HinderedRotorCM1D : public DensityOfStatesCalculator, protected HinderedRotorUtils
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
    HinderedRotorCM1D(const char* id) : m_id(id),
      m_bondID(),
      m_reducedMomentInertia(0.0),
      m_periodicity(1),
      m_energyLevels()
    { Register(); }

    virtual ~HinderedRotorCM1D() {}
    virtual const char* getID()  { return m_id; }
    virtual HinderedRotorCM1D* Clone() { return new HinderedRotorCM1D(*this); }

  private:
    const char* m_id;

    std::string m_bondID;

    double m_reducedMomentInertia;
    int    m_periodicity;

    vector<double> m_energyLevels ;	     // The energies of the hindered rotor states.

  } ;

  //-------------------------------------------------------------
  //Global instance, defining its id
  HinderedRotorCM1D theHinderedRotorCM1D("HinderedRotorCM1D");
  //-------------------------------------------------------------

  using OpenBabel::vector3;
  //Read data from XML and store in this instance.
  bool HinderedRotorCM1D::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
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
    cinfo << "Hindered rotor " << m_bondID << endl;  

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

    // Calculate reduced moment of inertia.

    m_reducedMomentInertia = reducedMomentInertia(gs, bondats, mode);  //units a.u.*Angstrom*Angstrom

    // Read in potential information.

	vector<double>&  potentialCosCoeff = get_PotentialCosCoeff() ;

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
        potentialCosCoeff.resize(++maxIndex) ;
        for (size_t i(0) ; i < coefficients.size() ; i++ ) {
          potentialCosCoeff[indicies[i]] = coefficients[i] ;
        }

      } else if (format == "numerical") {

        // Numerical potential.

        vector<double> potential ;
        vector<double> angle ;
        set_Expansion(pp->XmlReadInteger("expansionSize",optional));

		// Check if sine terms are to be used.

		const char *pUseSineTerms(pp->XmlReadValue("useSineTerms",optional)) ;
		if (pUseSineTerms && string(pUseSineTerms) == "yes") {
		  set_UseSinTerms(true) ;
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

        potentialCosCoeff.push_back(0.0) ;

      }

    } else {

      // Default : free rotor.

      cinfo << "No potential defined for " << bondID << ", assuming free rotor." <<endl;

      potentialCosCoeff.push_back(0.0) ;

    }

    return true;
  }

  //
  // Calculate classical mechanical 1D rotor densities of states of a free 
  // rotor and convolve them with the main density of states.
  //
  bool HinderedRotorCM1D::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {
    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, 0, false)) // retrieve the DOS vector without recalculating
      return false;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);

    // Calculate the one dimensional rotor constant and adjust for symmetry.
    double bint = sqrt(m_reducedMomentInertia/conMntInt2RotCnt)/double(m_periodicity) ;

    //
    // Calculate the free rotor density of states. The density of states goes as
    // (Energy)^(-1/2), howwever, to get a good estimate of the density of states
    // it is better to use the average of this function across the cell, hence the
    // slightly strange formula below.
    //
    vector<double> freeRtrDOS(MaximumCell,0.0) ;
    for (size_t i(0) ; i < MaximumCell ; i++ ) {
      const double ene = cellEne[i] - 0.5 ;
      freeRtrDOS[i] = 2.0*bint*(sqrt(ene+1.0)-sqrt(ene)) ;
    }

    //
    // Calculate the configuraton integral contribution to the density of states.
    // The configuration integral features a delta function of the differnce between
    // the required energy an the hindering potential and so it is necessary to 
    // determine the roots of this difference and the gradient of the potential at
    // these points.
    //
    // The configuation intergral evaluation is done in three steps:
    //   1) The potential is eveluated at number of angles.
    //   2) The configuration integral is determined for fine intervals of energy
    //      by bracketing and interpolating the energy and related functions.
    //   3) The cell (note cell not grain) values are determined by integration.
    //
    // Step 3) seems to be require as mid cell energies tend under estimate the 
    // the configuration integral contribution.
    //

	const vector<double>&  potentialCosCoeff = get_PotentialCosCoeff() ;

    // 1) Set-up array of potential points.
    const size_t npnts(2000) ;
    const double intvl(2.0*M_PI/double(npnts)) ;
    vector<double>  ptnl(npnts,0.0) ;
    vector<double> dptnl(npnts,0.0) ;
    for (size_t i(0); i < npnts; ++i) {
      double angle(double(i)*intvl) ;
      for(size_t k(0); k < get_Expansion() ; ++k) {
        double nTheta = double(k) * angle;
        ptnl[i]  +=  potentialCosCoeff[k] * cos(nTheta);
        dptnl[i] += -potentialCosCoeff[k] * double(k) * sin(nTheta);
      }
    }

    // 2) Locate roots via bracketing and interpolate the potential gradient.
    const int    nintvl   = 10 ;
    const double cellsize = 1.0 ;
    const double dene     = cellsize/double(nintvl) ;
    size_t       emax     = 2*nintvl*(int(potentialCosCoeff[0]) + 1) + 1 ;
    vector<double> cfgHdr(emax,0.0) ;
    for (size_t i(0); i < emax ; ++i) {
      const double ene = dene*double(i) ;
      for (size_t j(0); j < npnts; ++j) {
        const double v1(ptnl[j]), v2(ptnl[(j+1)%npnts]) ;
        if ((ene > v1 && ene < v2)||(ene > v2 && ene < v1)) {
          const double dp = fabs(dptnl[j] - (dptnl[j] - dptnl[j+1])*(v1 - ene)/(v1 - v2)) ;
          cfgHdr[i] += (dp > 0.0) ? 1.0/dp : 0.0 ;
        }
      }
      cfgHdr[i] /= 2.0*M_PI ;
    }

    // 3) Integrate using the trapezium rule to get an averages across a cells.
    vector<double> tmpCellDOS(MaximumCell,0.0) ;
    emax = 2*(int(potentialCosCoeff[0]) + 1) ;
    for (size_t i(0), idx(0); i < emax ; ++i) {
      double sum(0.0) ;
      sum += 0.5*cfgHdr[idx] ;
      for (size_t j(1); j < size_t(nintvl); ++j, ++idx) {
        sum += cfgHdr[idx] ;
      }
      sum += 0.5*cfgHdr[++idx] ;
      tmpCellDOS[i] = sum/double(nintvl) ;
    }

    //
    // Convolve free rotor and configuration integral terms to give the
    // hindered rotor density of states.
    //
    vector<double> hndrRtrDOS(MaximumCell,0.0) ;
    FastLaplaceConvolution(freeRtrDOS, tmpCellDOS, hndrRtrDOS) ;

    // Convolve one dimensional rotor states with overall density of states.
    FastLaplaceConvolution(cellDOS, hndrRtrDOS, tmpCellDOS) ;

    // Replace existing density of states.   
    pDOS->setCellDensityOfStates(tmpCellDOS) ;

    return true;

  }

  //
  // Provide a function to calculate contribution to canonical partition function.
  // (Mostly for testing purposes.)
  //
  double HinderedRotorCM1D::canPrtnFnCntrb(gDensityOfStates* gdos, double beta)
  {
    //
    // Calculate the free rotor term first.
    //
    double Qintrot = sqrt(M_PI*m_reducedMomentInertia/conMntInt2RotCnt/beta)/double(m_periodicity) ;

    //
    // Calculate the hindering potential correction via numerical integration.
    //
	const vector<double>&  potentialCosCoeff = get_PotentialCosCoeff() ;

    const size_t npnts(1000) ;
    const double intvl(2*M_PI/double(npnts)) ;
    double Qhdr(0.0) ;
    for (size_t i(0); i < npnts; ++i) {
      double ptnl(0.0) ;
      double angle(double(i)*intvl) ;
      for(size_t k(0); k < get_Expansion() ; ++k) {
        double nTheta = double(k) * angle;
        ptnl += potentialCosCoeff[k] * cos(nTheta);
      }
      Qhdr += exp(-beta*ptnl);
    }

    return Qintrot*Qhdr/double(npnts) ;
  }

}//namespace
