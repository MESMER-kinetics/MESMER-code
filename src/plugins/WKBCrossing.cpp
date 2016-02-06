//-------------------------------------------------------------------------------------------
//
// WKBCrossing.cpp
//
// Author: Dave Glowacki, based on fortran code written by Jeremy Harvey
// Date:   28-8-2009
//
// Produces WKB spin forbidden crossing coefficients
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
	class WKBCrossing : public MicroRateCalculator
	{
	public:
		virtual bool ParseData(PersistPtr pp);

		/********************************************************************************
		Constructor which registers this class with the list of plugins, initializes the ID
		and also does some initialization specific to this class.
		********************************************************************************/
		WKBCrossing(const char* id) : m_id(id),
			m_SOCelement(0.0),
			m_GradDiffMagnitude(0.0),
			m_ReducedMass(0.0),
			m_AverageSlope(0.0)
		{
			Register();
		}

		virtual ~WKBCrossing() {}

		virtual const char* getID() { return m_id; }

		virtual WKBCrossing* Clone() { return new WKBCrossing(*this); }

		virtual bool calculateMicroCnlFlux(Reaction* pReact);

	private:

		bool calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability);

		bool ReadDoubleAndUnits(double& element, PersistPtr pp, const string identifier, const string units);

		const char* m_id;
		double m_SOCelement, m_GradDiffMagnitude, m_ReducedMass, m_AverageSlope;
	};

	//************************************************************
	//Global instance, defining its id 
	WKBCrossing theWKBCrossing("WKBCrossing");
	//************************************************************

	bool WKBCrossing::ParseData(PersistPtr pp)
	{
		// Check for transition state.
		if (!m_parent->get_TransitionState())
		{
			cerr << "Reaction " << m_parent->getName()
				<< " uses the WKB Crossing method, which should have transition state." << endl;
			return false;
		}

		bool status(true);

		status = status && ReadDoubleAndUnits(m_SOCelement, pp, "me:RMS_SOC_element", "cm-1");
		status = status && ReadDoubleAndUnits(m_GradDiffMagnitude, pp, "me:GradientDifferenceMagnitude", "a.u./Bohr");
		status = status && ReadDoubleAndUnits(m_ReducedMass, pp, "me:GradientReducedMass", "a.m.u.");
		status = status && ReadDoubleAndUnits(m_AverageSlope, pp, "me:AverageSlope", "a.u./Bohr");

		return status;
	}

	bool WKBCrossing::calculateCellCrossingCoeffs(Reaction* pReact, vector<double>& CrossingProbability)
	{
		double ZPE_corr_barrier_height;

		// get threshold energy:
		ZPE_corr_barrier_height = pReact->get_ThresholdEnergy();

		const double SOCelementAU = m_SOCelement/Hartree_in_RC ;
		const double ReducedMassAU = m_ReducedMass * 1.822888e+3;

		//get properties of vectors in which to include Crossing coefficients
		const size_t MaximumCell = pReact->getEnv().MaxCell;
		CrossingProbability.clear();
		CrossingProbability.resize(MaximumCell, 0.0);

		//set transmission coefficients to 0 below the ZPE corrected barrier height;
		//above the barrier, the Landau Zener transmission coefficients are calculated 
		//as described by Harvey & Aschi, Faraday Discuss, 2003 (124) 129-143

		double Ai(0);

		const double DT_1 = 4.0*M_PI*M_PI*SOCelementAU*SOCelementAU*pow((2.0e+0 * ReducedMassAU / (m_GradDiffMagnitude * m_AverageSlope)), (2.0e+0 / 3.0e+0));
		const double DT_2 = pow((2.0e+0 * ReducedMassAU * pow(m_GradDiffMagnitude, 2.0e+0) / pow(m_AverageSlope, 4.0e+0)), (1.0e+0 / 3.0e+0));

		for (size_t i = 0; i < MaximumCell; ++i) {
			double E = double(i) - ZPE_corr_barrier_height;
			double E_AU = E / Hartree_in_RC;
			double xvalue = -E_AU * DT_2;
			airy2(xvalue, Ai);            //airy2 accurate over entire tunnelling regime
			// double Aip(0), Bi(0), Bip(0); //airy returns Ai, Bi, and their derivatives, Aip & Bip
			//airy(xvalue, Ai, Aip, Bi, Bip); //airy accurate in shallow tunnelling regime; less so for deep tunnelling
			CrossingProbability[i] = DT_1 * pow(Ai, 2.0e+0);
			// following if statement to avoid nan
			if (IsNan(CrossingProbability[i])) CrossingProbability[i] = 0.0;
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

	bool WKBCrossing::ReadDoubleAndUnits(double& element, PersistPtr pp, const string identifier, const string units)
	{
		PersistPtr ppData = pp->XmlMoveTo(identifier);
		element = pp->XmlReadDouble(identifier); // Or default.
		if (!IsNan(element) && ppData)
		{
			const char* unitsTxt = ppData->XmlReadValue("units", false);
			if (string(unitsTxt) != units) {
				cinfo << "MESMER could not read units for " << identifier << "; assuming " << units << "." << endl;
				cinfo << identifier << " = " << element << " " << units << endl;
			}
		}

		return true;
	}

	bool WKBCrossing::calculateMicroCnlFlux(Reaction* pReact)
	{
		Molecule* pTS = pReact->get_TransitionState();
		if (!pTS)
		{
			cerr << "Lack of transition state in reaction " << pReact->getName() << " for WKB Crossing" << endl;
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

		// The flux bottom energy is equal to the ZPE of the higher well.

		const size_t CrossingStart = size_t(max(0.0, pReact->getHeatOfReaction()));

		for (size_t i = CrossingStart; i < MaximumCell; ++i) {                       // Calculate flux using RRKM 
			rxnFlux[i - CrossingStart] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm; // with crossing correction
		}

		if (CrossingStart > 0) { pReact->setCellFluxBottom(pReact->get_relative_pdtZPE()); }
		else { pReact->setCellFluxBottom(pReact->get_relative_rctZPE()); }

		return true;
	}

}//namespace

