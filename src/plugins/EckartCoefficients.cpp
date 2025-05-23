//-------------------------------------------------------------------------------------------
//
// EckartCoefficients.h
//
// Author: Chi-Hsiu Liang
// Date:   _2008_03_19__10_02_55_
//
// Produces Eckart tunneling coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../System.h"
#include "../AssociationReaction.h"
#include "../IrreversibleExchangeReaction.h"
#include "../PseudoIsomerizationReaction.h"
#include "../gDensityOfStates.h"
#include "../Tunneling.h"

using namespace Constants;

namespace mesmer
{
  class EckartCoefficients : public TunnelingCalculator
  {
  public:

    ///Constructor which registers with the list of TunnelingCalculators in the base class
    EckartCoefficients(const char* id) : m_id(id) { Register(); }

    virtual ~EckartCoefficients() {}
    virtual const char* getID() { return m_id; }
    virtual EckartCoefficients* Clone() { return new EckartCoefficients(*this); }

    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability);
  private:
    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id
  EckartCoefficients theEckartCoefficients("Eckart");
  //************************************************************

  bool EckartCoefficients::calculateCellTunnelingCoeffs(Reaction* pReact, vector<double>& TunnelingProbability) {

    // First determine if V0 and V1 are measured from the classical well or zero point energy.
    PersistPtr ppEckart = pReact->get_PersistentPointer();
    ppEckart = ppEckart->XmlMoveTo("me:tunneling");
    bool useZPE = ppEckart->XmlReadBoolean("me:useZPE");

    double rctClassicalEnergy(pReact->get_reactant()->getDOS().getClassicalEnergy());
    if (pReact->getReactionType() == IRREVERSIBLE_EXCHANGE) {
      IrreversibleExchangeReaction* irex = static_cast<IrreversibleExchangeReaction*>(pReact);
      rctClassicalEnergy += irex->get_excessReactant()->getDOS().getClassicalEnergy();
    }
    else if (pReact->getReactionType() == ASSOCIATION || pReact->getReactionType() == PSEUDOISOMERIZATION) {
      AssociationReaction* asso = static_cast<AssociationReaction*>(pReact);
      rctClassicalEnergy += asso->get_excessReactant()->getDOS().getClassicalEnergy();
    }

    pReact->setUsesProductProperties();
    std::vector<Molecule *> pdt_temp;
    pReact->get_products(pdt_temp);
    const double pdtClassicalEnergy(pdt_temp.size() == 1
      ? pdt_temp[0]->getDOS().getClassicalEnergy()
      : pdt_temp[0]->getDOS().getClassicalEnergy() + pdt_temp[1]->getDOS().getClassicalEnergy());

    Molecule * p_TransitionState = pReact->get_TransitionState();

    //TC is the classical energy of the TS
    const double TC = p_TransitionState->getDOS().getClassicalEnergy();
    // Sanity check the barrier values and issue a warning if check fails
    if (rctClassicalEnergy > TC || pdtClassicalEnergy > TC) {
      throw std::runtime_error("Eckart tunnelling barrier is negative. Check classical energy correction.");
    }

    // TZ is the zpe of the TS
    const double TZ = pReact->get_relative_TSZPE();
    // barrier0 & barrier1 are the zpe corrected barrier heights in the forward/reverse directions.
    const int barrier0 = int(TZ - pReact->get_relative_rctZPE());
    const int barrier1 = int(TZ - pReact->get_relative_pdtZPE());

    double V0(0.0), V1(0.0);
    PersistPtr ppBarrier = ppEckart->XmlMoveTo("me:BarrierHeights");
    if (!ppBarrier) {
      if (useZPE) {
        // The zpe corrected barrier heights in the forward/reverse directions
        V0 = TZ - pReact->get_relative_rctZPE();
        V1 = TZ - pReact->get_relative_pdtZPE();
      }
      else {
        // V0 & V1 are the classical barrier heights in the forward/reverse directions.
        V0 = TC - rctClassicalEnergy;
        V1 = TC - pdtClassicalEnergy;
      }
    }
    else {
      const char* p = ppBarrier->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";
      V0 = ppBarrier->XmlReadDouble("V0", optional);
      V1 = ppBarrier->XmlReadDouble("V1", optional);

      V0 = getConvertedEnergy(units, V0);
      V1 = getConvertedEnergy(units, V1);
    }

    //imFreq is the imaginary frequency of the TS
    const double imFreq = pReact->get_TSImFreq();

    //get properties of vectors in which to include transmission coefficients
    const size_t MaximumCell = pReact->get_reactant()->getEnv().MaxCell;
    TunnelingProbability.clear();
    TunnelingProbability.resize(MaximumCell,0.0);

    // Set transmission coefficients to 0 where no tunneling is possible;
    // where tunneling may occur, the transmission coefficients are calculated
    // using a 1d eckart barrier as described by W.H. Miller, JACS, 101(23),
    // 1979. Note: the parameters a, b, and c defined below must be unitless.

    double tmp = (4.0 * M_PI / imFreq) / (1.0 / sqrt(V0) + 1.0 / sqrt(V1));
    double c = 2.0 * M_PI * sqrt(V0 * V1 / (imFreq*imFreq) - 1.0 / 16.0);
    double csh2c = cosh(c)*cosh(c);
    for (size_t i(0); i < MaximumCell; ++i) {
      int E = int(i) - barrier0;
      if ((E + barrier1) < 0) {
        TunnelingProbability[i] = 0.0;
      }
      else {
        double a = tmp * sqrt(double(E) + V0);
        double b = tmp * sqrt(double(E) + V1);
        TunnelingProbability[i] = (sinh(a) * sinh(b)) / (pow(sinh((a + b) / 2.0), 2.0) + csh2c);
        // following if statement to avoid nan at small values of E
        if (IsNan(TunnelingProbability[i])) 
          TunnelingProbability[i] = 0.0;
      }
    }

    if (pReact->getFlags().TunnellingCoeffEnabled) {
      ctest << "\nTunneling coefficients for: " << pReact->getName()
        << "\nV0 = " << V0 << ", V1 = " << V1
        << ", barrier0 = " << barrier0 << ", barrier1 = " << barrier1
        << ", imFreq = " << imFreq * SpeedOfLight_in_cm << " Hz\n{\n";
      for (size_t i(0); i < MaximumCell; ++i) {
        ctest << TunnelingProbability[i] << endl;
      }
      ctest << "}\n";
    }

    return true;
  }
}//namespace

