//-------------------------------------------------------------------------------------------
//
// SimpleRRKM.cpp
//
// Authors: Struan Robertson, David Glowacki and Alyn Liang
// Date:    Jul/2007
//
// Calculates microcanonical rate coefficients for a well defined transition state, 
// applying tunnelling corrections if requested.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  class SimpleRRKM : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    SimpleRRKM(const char* id) : m_id(id) { Register(); }
  
    virtual ~SimpleRRKM() {}
    virtual const char* getID()  { return m_id; }

    virtual SimpleRRKM* Clone() { return new SimpleRRKM(*this); }

    virtual bool ParseData(PersistPtr);
  
    virtual bool calculateMicroRateCoeffs(Reaction* pReact);

  private:
    const char* m_id;
  };

  //************************************************************
	//Global instance, defining its id (usually the only instance)
	SimpleRRKM theSimpleRRKM("SimpleRRKM");
	//************************************************************

	bool SimpleRRKM::calculateMicroRateCoeffs(Reaction* pReact)
	{
		Molecule* pTS = pReact->get_TransitionState();
		if(!pTS)
		{
			cerr << "Lack of transition state in reaction " << pReact->getName() << " for SimpleRRKM" << endl;
			return false;
		}
		// Allocate some work space for density of states.
		vector<double> TScellDOS; // Transistion state density of states.
		if(!pTS->getDOS().getCellDensityOfStates(TScellDOS))
			return false; // Extract densities of states from molecules.

		// get MaxCell from MesmerEnv structure via Reaction class
		const int MaximumCell = pReact->getEnv().MaxCell;

		// Allocate space to hold transition state flux and initialize elements to zero.
		vector<double>& rxnFlux = pReact->get_CellFlux();
		rxnFlux.clear();
		rxnFlux.resize(MaximumCell, 0.0);

		if (pReact->thereIsTunnelling()) { // with tunneling
			int HeatOfReaction = pReact->getHeatOfReactionInt();
			const int TunnelingStart = (HeatOfReaction > 0) ? int(HeatOfReaction) : 0;

			vector<double> TunnelingProbability;
			pReact->calculateCellTunnelingCoeffs(TunnelingProbability);

			vector<double> ConvolvedSumOfStates;
			FastLaplaceConvolution(TScellDOS, TunnelingProbability, ConvolvedSumOfStates); // FFT convolution
			//      Convolution(TScellDOS, TunnelingProbability, ConvolvedSumOfStates); // standard convolution

			for (int i = TunnelingStart; i < MaximumCell; ++i) {                       // Calculate flux using RRKM 
				rxnFlux[i-TunnelingStart] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm; // with tunneling correction
			}

			// the flux bottom energy is equal to the ZPE of the higher well
			if (TunnelingStart > 0) {
				pReact->setCellFluxBottom(pReact->get_relative_pdtZPE());
			}
			else{
				pReact->setCellFluxBottom(pReact->get_relative_rctZPE());
			}

		}
		else if(pReact->thereIsCrossing()){ //with spin forbidden crossing

			vector<double> CrossingProbability;
			pReact->calculateCellCrossingCoeffs(CrossingProbability);

			vector<double> ConvolvedSumOfStates;
			FastLaplaceConvolution(TScellDOS, CrossingProbability, ConvolvedSumOfStates); // FFT convolution

			const bool tunnelling = pReact->thereIsCrossingWithTunnelling();

			// the flux bottom energy is equal to the ZPE of the higher well if there is tunnelling & crossing
			if(tunnelling){

				int HeatOfReaction = pReact->getHeatOfReactionInt();
				const int CrossingStart = (HeatOfReaction > 0) ? int(HeatOfReaction) : 0;

				for (int i = CrossingStart; i < MaximumCell; ++i) {                       // Calculate flux using RRKM 
					rxnFlux[i-CrossingStart] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm; // with crossing correction
				}

				if (CrossingStart > 0) {pReact->setCellFluxBottom(pReact->get_relative_pdtZPE());}
				else{pReact->setCellFluxBottom(pReact->get_relative_rctZPE());}

			}
			else if(!tunnelling){  // what to do if it's just crossing with no tunnelling

				const int BarrierHeight = int(pReact->get_ThresholdEnergy());

				if(BarrierHeight >= 0.0){  // if the barrier is positive
					for (int i = BarrierHeight ; i < MaximumCell ; ++i)     // Calculate k(E)s using RRKM expression.
						rxnFlux[i-BarrierHeight] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm;
					pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + BarrierHeight);
				}
				else{  // if the barrier is negative
					for (int i=0 ; i < MaximumCell ; ++i)     // Calculate k(E)s using RRKM expression.
						rxnFlux[i] = ConvolvedSumOfStates[i] * SpeedOfLight_in_cm;
					//double dummy=pReact->get_relative_rctZPE() + BarrierHeight;
					pReact->setCellFluxBottom(pReact->get_relative_rctZPE());
				}

			}
		}

		else{ // if there's neither crossing nor tunneling
			double SumOfStates = 0.0;
			for (int i = 0 ; i < MaximumCell ; ++i) {
				// Integrate transition state density of states.
				SumOfStates += TScellDOS[i];

				// Calculate microcanonical rate coefficients using RRKM expression.
				rxnFlux[i] = SumOfStates * SpeedOfLight_in_cm;
			}

			// the flux bottom energy is equal to the ZPE of the transition state
			pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + pReact->get_ThresholdEnergy());
		}

		return true;
	}

	bool SimpleRRKM::ParseData(PersistPtr)
	{
		if (!m_parent->get_TransitionState())
		{
			cerr << "Reaction " <<m_parent->getName() 
				<< " uses RRKM method, which should have transition state."<<endl;
			return false;
		}
		return true;
	}

}//namespace
