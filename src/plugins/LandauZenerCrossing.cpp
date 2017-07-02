//-------------------------------------------------------------------------------------------
//
// LandauZenerCrossing.cpp
//
// Author: Dave Glowacki, based on fortran code written by Jeremy Harvey
// Date:   28-8-2009
//
// Produces Landau-Zener spin forbidden crossing coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../System.h"
#include "../gDensityOfStates.h"

using namespace Constants;
using namespace std;

namespace mesmer
{
  class LandauZenerCrossing : public MicroRateCalculator
  {
  public:
    virtual bool ParseData(PersistPtr pp);

    /********************************************************************************
    Constructor which registers this class with the list of plugins, initializes the ID
    and also does some initialization specific to this class.
    ********************************************************************************/
    LandauZenerCrossing(const char* id) : m_id(id),
      m_SOCelement(0.0),
      m_GradDiffMagnitude(0.0),
      m_ReducedMass(0.0)
    {
      Register();
    }

    virtual ~LandauZenerCrossing() {}

    virtual const char* getID() { return m_id; }

    virtual LandauZenerCrossing* Clone() { return new LandauZenerCrossing(*this); }

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

  private:

		bool calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability);

		bool ReadDoubleAndUnits(double& element, PersistPtr pp, const string identifier, const string units);

    const char* m_id;
    double m_SOCelement, m_GradDiffMagnitude, m_ReducedMass;
  };

  //************************************************************
  //Global instance, defining its id
  LandauZenerCrossing theLandauZenerCrossing("LandauZenerCrossing");
  //************************************************************

  bool LandauZenerCrossing::ParseData(PersistPtr pp)
  {
		// Check for transition state.
		if (!m_parent->get_TransitionState())
		{
			cerr << "Reaction " << m_parent->getName()
				<< " uses Landau-Zener crossing method, which should have transition state." << endl;
			return false;
		}

		bool status(true);

		status = status && ReadDoubleAndUnits(m_SOCelement, pp, "me:RMS_SOC_element", "cm-1");
		status = status && ReadDoubleAndUnits(m_GradDiffMagnitude, pp, "me:GradientDifferenceMagnitude", "a.u./Bohr");
		status = status && ReadDoubleAndUnits(m_ReducedMass, pp, "me:GradientReducedMass", "a.m.u.");

		return status;
  }

  bool LandauZenerCrossing::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability) {

    // get threshold energy:
    double ZPE_corr_barrier_height;
    ZPE_corr_barrier_height = pReact->get_ThresholdEnergy();

    const double SOCelementAU = m_SOCelement / Hartree_in_RC ;
    const double ReducedMassAU = m_ReducedMass * 1.822888e+3;

    // Get properties of vectors in which to include Crossing coefficients.
    const size_t MaximumCell = pReact->getEnv().MaxCell;
    CrossingProbability.clear();
    CrossingProbability.resize(MaximumCell, 0.0);

    // Set transmission coefficients to 0 below the ZPE corrected barrier height;
    // above the barrier, the Landau Zener transmission coefficients are calculated 
    // as described by Harvey & Aschi, Faraday Discuss, 2003 (124) 129-143.

		const double Cnst = 2.0 * M_PI * SOCelementAU*SOCelementAU / m_GradDiffMagnitude ;
    for (size_t i = 0; i < MaximumCell; ++i) {

      double E(double(i) + 0.5);  // Set E to the avg energy of the cell in cm-1.
      if (ZPE_corr_barrier_height >= 0.0) { 
        E -= ZPE_corr_barrier_height ;
      }

      if (E >= 0.0) {
        double E_AU = E / Hartree_in_RC;
        double trans_probability = exp(-Cnst * sqrt(ReducedMassAU / (2.0 * E_AU)));
        CrossingProbability[i] = (1.0 - trans_probability*trans_probability) ;
      }
    }

    if (pReact->getFlags().CrossingCoeffEnabled) {
      ctest << "\nCrossing coefficients for: " << pReact->getName() << endl;
      for (size_t i = 0; i < MaximumCell; ++i) {
        ctest << CrossingProbability[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }

  bool LandauZenerCrossing::ReadDoubleAndUnits(double& element, PersistPtr pp, const string identifier, const string units)
  {
    PersistPtr ppData = pp->XmlMoveTo(identifier);
    element = pp->XmlReadDouble(identifier); // Or default.
    if (!IsNan(element) && ppData)
    {
      const char* unitsTxt = ppData->XmlReadValue("units");
      if (string(unitsTxt) != units) {
        cinfo << "MESMER could not read units for " << identifier << "; assuming " << units << "." << endl;
        cinfo << identifier << " = " << element << " " << units << endl;
      }
    }

		return true;
  }

  bool LandauZenerCrossing::calculateMicroCnlFlux(Reaction* pReact)
  {
    Molecule* pTS = pReact->get_TransitionState();
    if (!pTS)
    {
      cerr << "Lack of transition state in reaction " << pReact->getName() << " for Landau-Zener Crossing" << endl;
      return false;
    }

    // Allocate some work space for transistion state density of states.
    vector<double> TScellDOS;
    if (!pTS->getDOS().getCellDensityOfStates(TScellDOS))
      return false;

    // Extract densities of states from molecules.
    const size_t MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    vector<double> CrossingProbability;
		calculateCellCrossingCoeffs(pReact, CrossingProbability);

    vector<double> ConvolvedSumOfStates;
    FastLaplaceConvolution(TScellDOS, CrossingProbability, ConvolvedSumOfStates); // FFT convolution

    const size_t BarrierHeight = size_t(max(0.0,pReact->get_ThresholdEnergy()));

      for (size_t i = BarrierHeight; i < MaximumCell; ++i)     // Calculate k(E)s using RRKM expression.
        rxnFlux[i - BarrierHeight] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm;

      pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + BarrierHeight);

    return true;
  }
}//namespace
