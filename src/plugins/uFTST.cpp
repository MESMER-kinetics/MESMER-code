//-------------------------------------------------------------------------------------------
//
// uFTST.cpp
//
// Author: Struan Robertson
// Date:   Feb/2019
//
// Calculates microcanonical rate coefficients for a barrilerless reaction coordinate.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../MesmerMath.h"
#include "../System.h"
#include "../gDensityOfStates.h"
#include "ParseForPlugin.h"
#include "PhaseIntergrals.h"
#include "FTSTPotential.h"

using namespace std;
using namespace Constants;
namespace mesmer
{

  class uFTST : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    uFTST(const char* id) : m_id(id),
      m_Frag1(NULL),
      m_Frag2(NULL),
      m_Adduct(NULL),
      m_nRxnCrdSteps(40),
      m_rxnCrdMin(0.5),
      m_rxnCrdMax(4.5),
      m_threshold(10000.0),
      m_R0(1.5),
      m_Beta0(1.0),
      m_addFrq(),
      m_frgFrq(),
      m_alpha(1.0),
      m_pPhaseIntegral(NULL),
      m_pFTSTPotential(NULL)
    {
      Register();
    }

    virtual ~uFTST() {}

    virtual const char* getID() { return m_id; }

    virtual uFTST* Clone() { return new uFTST(*this); }

    virtual bool ParseData(PersistPtr);

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual double get_ThresholdEnergy(Reaction* pReact);

  private:

    // This function returns the minimum energy potential for a given reaction coordinate value.
    double MinimumEnergyPotential(double rxnCrd) const;

    // Convolve the conserved modes into the sum of states.
    bool ConvolveConservedModes(double rxnCrd, vector<double> &wrk) const;

    // Determine the top combination,
    bool TopCombination();

    const char* m_id;

    Molecule* m_Frag1;
    Molecule* m_Frag2;
    Molecule* m_Adduct;

    size_t m_nRxnCrdSteps;
    double m_rxnCrdMin;
    double m_rxnCrdMax;

    // Reaction coordinate potential parameters.
    double m_threshold;
    double m_R0;
    double m_Beta0;

    vector<double> m_addFrq;
    vector<double> m_frgFrq;
    double m_alpha;

    PhaseIntegral *m_pPhaseIntegral;
    FTSTPotential *m_pFTSTPotential;

  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  uFTST theuFTST("uFTST");
  //************************************************************

  bool uFTST::ParseData(PersistPtr pp)
  {
    // First check that the reaction type is sensible.
    if (m_parent->getReactionType() == ISOMERIZATION || m_parent->getReactionType() == IRREVERSIBLE_EXCHANGE) {
      cerr << "Reaction " << m_parent->getName()
        << " cannot use uFTST method because it is not a barrierless reaction." << endl;
      return false;
    }

    // Reaction coordinate limits.
    double RxnCoordMin = pp->XmlReadDouble("me:RxnCoordMin", optional);
    if (!IsNan(RxnCoordMin))
      m_rxnCrdMin = RxnCoordMin;

    double RxnCoordMax = pp->XmlReadDouble("me:RxnCoordMax", optional);
    if (!IsNan(RxnCoordMax))
      m_rxnCrdMax = RxnCoordMax;

    int RxnCoordSteps = pp->XmlReadInteger("me:RxnCoordSteps", optional);
    if (!IsNan(RxnCoordSteps))
      m_nRxnCrdSteps = RxnCoordSteps;

    // Reaction coordinate potential (Morse Potential)
    PersistPtr ppRxnPtnl = pp->XmlMoveTo("me:MorseRxnCoordPotential");
    double threshold = ppRxnPtnl->XmlReadDouble("me:Threshold");
    if (!IsNan(threshold) && threshold > 0.0) {
      PersistPtr ppThreshold = ppRxnPtnl->XmlMoveTo("me:Threshold");
      const char* p = ppThreshold->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";
      m_threshold = getConvertedEnergy(units, threshold);
    }
    else {
      cerr << "Reaction " << m_parent->getName()
        << " Morse potential parameter Threshold should be greater than 0.0." << endl;
      return false;
    }

    double R0 = ppRxnPtnl->XmlReadDouble("me:EquilibriumLength");
    if (!IsNan(R0))
      m_R0 = R0;

    // Adjust minimum reaction coordinate value if necessay.
    m_rxnCrdMin = max(m_rxnCrdMin, R0);

    double Beta0 = ppRxnPtnl->XmlReadDouble("me:ReciprocalBeta");
    if (!IsNan(Beta0))
      m_Beta0 = Beta0;
    if (m_Beta0 > 0.0)
      m_Beta0 = 1.0 / m_Beta0;
    else {
      cerr << "Reaction " << m_parent->getName()
        << " Morse potential parameter ReciprocalBeta should be greater than 0.0." << endl;
      return false;
    }

    // Conserved modes: read in frequency mappings.
    PersistPtr ppConMod = pp->XmlMoveTo("me:ConservedModeMap");
    double alpha = ppConMod->XmlReadDouble("alpha", optional);
    if (!IsNan(alpha))
      m_alpha = alpha;

    while (ppConMod = ppConMod->XmlMoveTo("me:FrequencyMapping"))
    {
      double addFrq = ppConMod->XmlReadDouble("adduct", optional);
      if (!IsNan(addFrq))
        m_addFrq.push_back(addFrq);

      double frgFrq = ppConMod->XmlReadDouble("Fragments", optional);
      if (!IsNan(frgFrq))
        m_frgFrq.push_back(frgFrq);
    }

    // Sanity check: number of parent and child conserved mode frequncies should be the same.
    if (m_addFrq.size() != m_frgFrq.size()) {
      cerr << "Reaction " << m_parent->getName()
        << " cannot use uFTST method because the conserved modes are not properly defined." << endl;
      return false;
    }

    vector<Molecule*> molecules;
    int nRcts = m_parent->get_reactants(molecules);
    if (nRcts == 2) {
      m_Frag1 = molecules[0];
      m_Frag2 = molecules[1];
    }
    else {
      cerr << "Reaction " << m_parent->getName()
        << " cannot use uFTST method because the numnber of dissociation fragments is not two." << endl;
      return false;
    }
    molecules.clear();
    int nPdts = m_parent->get_products(molecules);
    if (nPdts == 1) {
      m_Adduct = molecules[0];
    }
    else {
      cerr << "Reaction " << m_parent->getName()
        << " cannot use uFTST method because there is no single adduct." << endl;
      return false;
    }

    // Sanity check: do we have the correct number of frequencies?
    size_t nCndMd = m_Frag1->getDOS().get_NoVibFreq() + m_Frag2->getDOS().get_NoVibFreq();
    if (nCndMd != m_addFrq.size()) {
      cerr << "Reaction " << m_parent->getName()
        << ": inconsistent number of conserved modes in uFTST method" << endl;
      return false;
    }

    m_pFTSTPotential = ParseForPlugin<FTSTPotential>(this, "me:FTSTPotential", pp);

    return true;
  }


  bool uFTST::calculateMicroCnlFlux(Reaction* pReact)
  {
    // get MaxCell from MesmerEnv structure via Reaction class
    const size_t MaximumCell = pReact->getEnv().MaxCell;
    const double cellSize = pReact->getEnv().CellSize;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);
    vector<double> OptRxnCrd(MaximumCell, m_rxnCrdMin);

    // Determine the combination tops and therefore the number of transitional modes.
    TopCombination();

    m_pPhaseIntegral->initialize(m_Frag1, m_Frag2, pReact);

    // Main loop over the reaction coordinate in Angrstroms.
    double rxnCrd(m_rxnCrdMin);
    double drxnCrd((m_rxnCrdMax - m_rxnCrdMin) / double(m_nRxnCrdSteps));
    for (size_t i(0); i < m_nRxnCrdSteps; i++) {

      // Some work space.
      vector<double> wrk(MaximumCell, 0.0);

      // Phase space integral of transitional modes.

      m_pPhaseIntegral->integrate(rxnCrd, wrk);

      // Convolve the conserved modes into the sum of states.

      ConvolveConservedModes(rxnCrd, wrk);

      // Adjust the sum of states to account for the reaction coordinate.

      double cellSize = m_parent->getEnv().CellSize;
      size_t imep = size_t((m_threshold - m_pFTSTPotential->MEPPotential(rxnCrd)) / cellSize);
      //size_t imep = size_t((m_threshold - MinimumEnergyPotential(rxnCrd)) / cellSize);

      if (i == 0) {
        for (size_t j(0), jj(imep); jj < wrk.size(); j++, jj++) {
          if (wrk[jj] > 0.0)
            rxnFlux[j] = wrk[jj];
          OptRxnCrd[j] = rxnCrd;
        }
      }
      else {
        for (size_t j(0), jj(imep); jj < wrk.size(); j++, jj++) {
          if (rxnFlux[j] > wrk[jj] && wrk[jj] > 0.0) {
            rxnFlux[j] = wrk[jj];
            OptRxnCrd[j] = rxnCrd;
          }
        }
      }

      // Move on to next reaction coordinate value.

      rxnCrd += drxnCrd;
    }

    ctest << endl;
    ctest << "  115 " << formatFloat(rxnFlux[154],  6, 14) << formatFloat(OptRxnCrd[154],  6, 14) << endl;
    ctest << "  248 " << formatFloat(rxnFlux[247],  6, 14) << formatFloat(OptRxnCrd[247],  6, 14) << endl;
    ctest << "  413 " << formatFloat(rxnFlux[412],  6, 14) << formatFloat(OptRxnCrd[412],  6, 14) << endl;
    ctest << "  624 " << formatFloat(rxnFlux[623],  6, 14) << formatFloat(OptRxnCrd[623],  6, 14) << endl;
    ctest << " 1040 " << formatFloat(rxnFlux[1039], 6, 14) << formatFloat(OptRxnCrd[1039], 6, 14) << endl;
    ctest << " 1248 " << formatFloat(rxnFlux[1247], 6, 14) << formatFloat(OptRxnCrd[1247], 6, 14) << endl;
    ctest << " 2007 " << formatFloat(rxnFlux[2006], 6, 14) << formatFloat(OptRxnCrd[2006], 6, 14) << endl;
    ctest << " 2080 " << formatFloat(rxnFlux[2079], 6, 14) << formatFloat(OptRxnCrd[2079], 6, 14) << endl;
    ctest << " 2600 " << formatFloat(rxnFlux[2599], 6, 14) << formatFloat(OptRxnCrd[2599], 6, 14) << endl;
    ctest << " 3120 " << formatFloat(rxnFlux[3119], 6, 14) << formatFloat(OptRxnCrd[3119], 6, 14) << endl;
    ctest << " 3419 " << formatFloat(rxnFlux[3418], 6, 14) << formatFloat(OptRxnCrd[3418], 6, 14) << endl;
    ctest << " 4015 " << formatFloat(rxnFlux[4014], 6, 14) << formatFloat(OptRxnCrd[4014], 6, 14) << endl;
    ctest << " 4445 " << formatFloat(rxnFlux[4444], 6, 14) << formatFloat(OptRxnCrd[4444], 6, 14) << endl;
    ctest << " 5554 " << formatFloat(rxnFlux[5553], 6, 14) << formatFloat(OptRxnCrd[5553], 6, 14) << endl;
    ctest << endl;


    // Calculate the flux.
    for (size_t i(0); i < MaximumCell; ++i) {
      rxnFlux[i] *= SpeedOfLight_in_cm;
    }

    // The flux bottom energy is equal to the ZPE of the transition state
    pReact->setCellFluxBottom(pReact->get_relative_rctZPE());

    return true;
  }

  // This function returns the minimum energy potential for a given reaction coordinate value.
  double uFTST::MinimumEnergyPotential(double rxnCrd) const {

    double potential(0.0);

    double A3    = -1.017234e-03;  // Hirst Surface parameters as fit by Hase et al.
    double A2    =  7.7738886e-02; // A3, A2, A1, AO parameters of the variable
    double A1    =  7.703640e-02;  // beta of the Modified Morse function.
    double A0    =  1.686690e0;    // 

    double delR  = rxnCrd - 1.0015;
    double delR2 = delR * delR ;
    double delR3 = delR * delR2;

    double bb  = A3 * delR3 + A2 * delR2 + A1 * delR + A0;
    double tmp = 1.0 - exp(-bb * delR);
    potential  = m_threshold * tmp*tmp;

    return potential;
  }


  // Convolve the conserved modes into the sum of states.
  bool uFTST::ConvolveConservedModes(double rxnCrd, vector<double> &wrk) const {

    //First calculate freqencies.
    vector<double> Freq;
    for (size_t j(0); j < m_addFrq.size(); j++) {
      double frq = m_frgFrq[j] + (m_addFrq[j] - m_frgFrq[j])*exp(-m_alpha * (rxnCrd - 1.09));
      Freq.push_back(frq);
    }

    // Implementation of the Beyer-Swinehart algorithm.
    double cellSize = m_parent->getEnv().CellSize;
    const size_t MaximumCell = m_parent->getEnv().MaxCell;
    for (size_t n(0); n < Freq.size(); ++n) {
      size_t freq = static_cast<size_t>(nint(Freq[n] / cellSize));
      if (freq > MaximumCell) {
        // This is to catch those occassional cases where the first excited 
        // vibrational state is above the cutoff, which can occur at low 
        // temperatures. 
        continue;
      }
      for (size_t m(0); m < MaximumCell - freq; ++m) {
        wrk[m + freq] += wrk[m];
      }
    }

    return true;
  }

  bool uFTST::TopCombination() {

    RotationalTop top1 = m_Frag1->getDOS().get_rotType();
    RotationalTop top2 = m_Frag2->getDOS().get_rotType();

    if (top1 == NONLINEAR && top2 == NONLINEAR) {
      m_pPhaseIntegral = new NLnrNLnrTops();
    }
    else if ((top1 == NONLINEAR && top2 == LINEAR) || (top1 == LINEAR && top2 == NONLINEAR)) {
      m_pPhaseIntegral = new NLnrLnrTops();
    }
    else if ((top1 == NONLINEAR && top2 == ATOMIC) || (top1 == ATOMIC && top2 == NONLINEAR)) {
      m_pPhaseIntegral = new NLnrAtmTops();
    }
    else if ((top1 == LINEAR && top2 == LINEAR)) {
      m_pPhaseIntegral = new LnrLnrTops();
    }
    else if ((top1 == LINEAR && top2 == ATOMIC) || (top1 == ATOMIC && top2 == LINEAR)) {
      m_pPhaseIntegral = new LnrAtmTops();
    }
    else if ((top1 == ATOMIC && top2 == ATOMIC)) {
      m_pPhaseIntegral = new AtmAtmTops();
    }

    return true;
  }

  double uFTST::get_ThresholdEnergy(Reaction* pReac) {

    return (pReac->get_relative_rctZPE() - pReac->get_relative_pdtZPE());

  }

}//namespace
