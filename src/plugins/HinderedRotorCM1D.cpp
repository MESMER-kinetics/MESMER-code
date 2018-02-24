
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
#include "../gDensityOfStates.h"
#include "../gStructure.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class HinderedRotorCM1D : public HinderedRotorUtils
  {
  public:
    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Provide a function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env);

    // Provide a function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) {return 1 ; } ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    // This class is an extra DOS class: a non-extra DensityOfStatesCalculator class also
    // needs to be specified.
    HinderedRotorCM1D(const char* id) : 
    HinderedRotorUtils(id),
      m_reducedMomentInertia(0.0),
      m_periodicity(1),
      m_energyLevels()
    { }

    virtual ~HinderedRotorCM1D() {}
    virtual HinderedRotorCM1D* Clone() { return new HinderedRotorCM1D(*this); }

  private:

		vector<double> m_potentialCosCoeff; // The cosine coefficients of the hindered rotor potential.
		vector<double> m_potentialSinCoeff; // The sine coefficients of the hindered rotor potential.

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

    m_reducedMomentInertia = gs.reducedMomentInertia(bondats);  //units a.u.*Angstrom*Angstrom

    // Read in potential information.

    m_periodicity = max(m_periodicity, ppDOSC->XmlReadInteger("me:periodicity",optional));

		ReadPotentialParameters(ppDOSC, string(bondID), m_potentialCosCoeff, m_potentialSinCoeff);

		// Check if there is a Hessian and knock out the frequency
		// associated with this internal rotation.

		if (gdos->hasHessian()) {
			// The following vector, "mode", will be used to hold the internal rotation 
			// mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010). 
			vector<double> mode(3 * gs.NumAtoms(), 0.0);
			gs.internalRotationVector(get_BondID(), mode);
			if (!gdos->projectMode(mode)) {
				cerr << "Failed to project out internal rotation." << endl;
				return false;
			}
		}

    return true;
  }

  //
  // Calculate classical mechanical 1D rotor densities of states of a free 
  // rotor and convolve them with the main density of states.
  //
  bool HinderedRotorCM1D::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
    const size_t MaximumCell = env.MaxCell ;
    const double cellSize = env.CellSize ;

    vector<double> cellDOS;
    if(!pDOS->getCellDensityOfStates(cellDOS, false)) // retrieve the DOS vector without recalculating
      return false;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellSize, cellEne);

    // Calculate the one dimensional rotor constant and adjust for symmetry.
    double bint = sqrt(m_reducedMomentInertia/conMntInt2RotCnt)/double(m_periodicity) ;

    //
    // Calculate the free rotor density of states. The density of states goes as
    // (Energy)^(-1/2), however, to get a good estimate of the density of states
    // it is better to use the average of this function across the cell, hence the
    // slightly strange formula below.
    //
    vector<double> freeRtrDOS(MaximumCell,0.0) ;
    for (size_t i(0) ; i < MaximumCell ; i++ ) {
      //const double ene = cellEne[i] - 0.5*cellSize ;
      //freeRtrDOS[i] = 2.0*bint*(sqrt(ene + cellSize)-sqrt(ene)) ;
			freeRtrDOS[i] = bint / sqrt(cellEne[i] - 0.37*cellSize);
		}

    //
    // Calculate the configuraton integral contribution to the density of states.
    // The configuration integral features a delta function of the differnce between
    // the required energy an the hindering potential and so it is necessary to 
    // determine the roots of this difference and the gradient of the potential at
    // these points.
    //
    // The configuation intergral evaluation is done in two steps:
    //   1) The potential is evaluated at number of angles.
    //   2) The configuration integral is determined for finite intervals of energy
    //      by bracketing and interpolating the energy and related functions.
		//

    // 1) Set-up array of potential points.
    const size_t npnts(200000) ;
    const double intvl(2.0*M_PI/double(npnts)) ;
    vector<double>  ptnl(npnts,0.0) ;
		double emax(-1.e+10);
		double emin(1.e+10);
    for (size_t i(0); i < npnts; ++i) {
      double angle(double(i)*intvl) ;
      for(size_t k(0); k < get_Expansion() ; ++k) {
        double nTheta = double(k) * angle;
        ptnl[i]  +=  m_potentialCosCoeff[k] * cos(nTheta);
      }
			emax = max(emax, ptnl[i]);
			emin = min(emin, ptnl[i]);
		}

    // 2) Locate roots via bracketing and determine the potential gradient.
    size_t nemax = 2*size_t(emax-emin) ;
		vector<double> tmpCellDOS(MaximumCell, 0.0);
    for (size_t i(0); i < nemax ; ++i) {
			const double ene = cellEne[i] + emin;
			double sum(0.0);
      for (size_t j(0); j < npnts; ++j) {
        const double v1(ptnl[j]), v2(ptnl[(j+1)%npnts]) ;
        if ((ene > v1 && ene < v2)||(ene > v2 && ene < v1)) {
					double angle = intvl*double(j) + intvl*(ene - v1)/(v2 -v1);
					double dp(0.0);
					for (size_t k(1); k < get_Expansion(); ++k) {
						double nTheta = double(k) * angle;
						dp += -m_potentialCosCoeff[k] * double(k) * sin(nTheta);
					}
          sum += fabs(1.0/dp) ;
        }
      }
			tmpCellDOS[i] = sum / (2.0*M_PI);
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
  void HinderedRotorCM1D::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) {

    //
    // Calculate the free rotor term first.
    //
    double Qintrot = 0.5*sqrt(m_reducedMomentInertia/conMntInt2RotCnt/beta/M_PI)/double(m_periodicity) ;

    //
    // Calculate the hindering potential correction via numerical integration.
    //
    const size_t npnts(10000) ;
    const double intvl(2*M_PI/double(npnts)) ;
    double Qhdr(0.0), Ehdr(0.0), varEhdr(0.0) ;
    for (size_t i(0); i < npnts; ++i) {
      double ptnl(0.0) ;
      double angle(double(i)*intvl) ;
      for(size_t k(0); k < get_Expansion() ; ++k) {
        double nTheta = double(k) * angle;
        ptnl += m_potentialCosCoeff[k] * cos(nTheta);
      }
      Qhdr    += exp(-beta*ptnl);
      Ehdr    += ptnl*exp(-beta*ptnl);
      varEhdr += ptnl*ptnl*exp(-beta*ptnl);
    }
    Ehdr    /= Qhdr ;
    varEhdr  = varEhdr/Qhdr - Ehdr*Ehdr ;

    PrtnFn   *= Qintrot*Qhdr*intvl;
    IntrlEne += 1.0/(2.0*beta) + Ehdr ;
    varEne   += 1.0/ (2.0*beta*beta) + varEhdr ;
  }

}//namespace
