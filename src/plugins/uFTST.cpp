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
#include "../ParseForPlugin.h"
#include "PhaseIntegrals.h"
#include "../FTSTPotential.h"
#include "../MicroRate.h"

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
      m_pFTSTPotential(NULL),
      m_namedPhaseIntegral(),
      m_writeSOS(false),
      m_eneForSOS()
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

    // Testing parameters
    string m_namedPhaseIntegral;
    bool m_writeSOS;
    vector<double> m_eneForSOS;

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
    size_t nRcts = m_parent->get_reactants(molecules);
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
    size_t nPdts = m_parent->get_products(molecules);
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

    m_pFTSTPotential->Initialize();

    // The following items are for test purposes only and should not be exposed in documentation.

    // Check to see if there is an analytic phase intergral requested.
    const char* ptxt = pp->XmlReadValue("me:NamedPhaseIntegral", optional);
    if (ptxt)
      m_namedPhaseIntegral = string(ptxt);

    // Check to see if sum of states are to be printed and for what energies.
    PersistPtr ppWrtSOS = pp->XmlMoveTo("me:WriteSumsOfStates");
    if (ppWrtSOS) {
      m_writeSOS = true;
      const char* p = ppWrtSOS->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";

      p = pp->XmlReadValue("me:WriteSumsOfStates");
      if (p) {
        istringstream idata(p);
        double x;
        while (idata >> x)
          m_eneForSOS.push_back(getConvertedEnergy(units, x));
      }
    }

    return true;
  }

  bool uFTST::calculateMicroCnlFlux(Reaction* pReact)
  {
    // get MaxCell from MesmerEnv structure via Reaction class
    const size_t MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);
    vector<double> OptRxnCrd(MaximumCell, m_rxnCrdMin);

    // Determine the combination tops and therefore the number of transitional modes.
    TopCombination();

    m_pPhaseIntegral->initialize(m_Frag1, m_Frag2, m_pFTSTPotential, pReact);

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

    if (m_writeSOS) {
      ctest << endl;
      if (m_eneForSOS.size() > 0) {
        for (size_t i = 0; i < m_eneForSOS.size(); i++) {
          size_t ii = size_t(m_eneForSOS[i]);
          if (ii < rxnFlux.size()) {
            ctest << setw(6) << ii << formatFloat(rxnFlux[ii], 6, 14) << formatFloat(OptRxnCrd[ii], 6, 14) << endl;
          }
        }
      }
      else {
        for (size_t i = 0; i < rxnFlux.size(); i++) {
          ctest << setw(6) << i << formatFloat(rxnFlux[i], 6, 14) << formatFloat(OptRxnCrd[i], 6, 14) << endl;
        }
      }
      ctest << endl;
    }

    // Calculate the flux.
    for (size_t i(0); i < MaximumCell; ++i) {
      rxnFlux[i] *= SpeedOfLight_in_cm;
    }

    // The flux bottom energy is equal to the ZPE of the transition state
    pReact->setCellFluxBottom(pReact->get_relative_rctZPE());

    return true;
  }

  // Convolve the conserved modes into the sum of states.
  bool uFTST::ConvolveConservedModes(double rxnCrd, vector<double> &wrk) const {

    //First calculate freqencies.
    vector<double> Freq;
    for (size_t j(0); j < m_addFrq.size(); j++) {
      double frq = m_frgFrq[j] + (m_addFrq[j] - m_frgFrq[j])*exp(-m_alpha * (rxnCrd - m_R0));
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

    if (m_namedPhaseIntegral.length() > 0) {
      if (m_namedPhaseIntegral == "MethylPlusH_HW")
        m_pPhaseIntegral = new MethylPlusH_HW();
      else
        throw(std::runtime_error("__FUNCTION__: No named phase integral known"));
    }
    else if (top1 == NONLINEAR && top2 == NONLINEAR) {
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
