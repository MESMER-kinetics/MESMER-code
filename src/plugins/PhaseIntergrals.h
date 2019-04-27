//-------------------------------------------------------------------------------------------
// PhaseIntergrals.h
//
// Author: Struan Robertson
// Date:   Feb/2019
//
// This file contains the definitions of the various phase integrals needed for
// microcanonical (NOT J resolved) flexible transition state theory.
//-------------------------------------------------------------------------------------------

#ifndef GUARD_PhaseIntergrals_h
#define GUARD_PhaseIntergrals_h

#include "../Molecule.h"
#include "gStructure.h"

namespace mesmer
{
  class FTSTPotential {

  public:

    virtual void initialize(double rxnCrd) = 0;

    virtual double MEPPotential(double rxnCrd) = 0;

    virtual double HinderingPotential(double rxnCrd) = 0;

  };

  class HirstWardlawPotential : public FTSTPotential {

  public:

    HirstWardlawPotential() : m_V0(0.0),
      m_currentRxnCrd(-1.0),
      m_threshold(0.0) {}
    virtual ~HirstWardlawPotential() {}

    void initialize(double rxnCrd) {
      double ant = getConvertedEnergy("kcal/mol", 164.0);
      double bnt = 0.43;
      m_V0 = ant * exp(-bnt * rxnCrd*rxnCrd);
    }

    virtual double MEPPotential(double rxnCrd) {
      double potential(0.0);

      double A3 = -1.017234e-03; // Hirst Surface parameters as fit by Hase et al.
      double A2 = 7.7738886e-02; // A3, A2, A1, AO parameters of the variable
      double A1 = 7.703640e-02;  // beta of the Modified Morse function.
      double A0 = 1.686690e0;    // 

      double delR  = rxnCrd - 1.0015;
      double delR2 = delR * delR;
      double delR3 = delR * delR2;

      double bb = A3 * delR3 + A2 * delR2 + A1 * delR + A0;
      double tmp = 1.0 - exp(-bb * delR);
      potential = m_threshold * tmp*tmp;

      return potential;

    }

    virtual double HinderingPotential(double rxnCrd, double gamma) {
      if (rxnCrd != m_currentRxnCrd)
        initialize(rxnCrd);
      double tmp = sin(gamma);
      return  m_V0 * tmp * tmp;
    }

  private:

      double m_V0;
      double m_currentRxnCrd;
      double m_threshold;

  };

  class PhaseIntegral {

  public:

    PhaseIntegral() : m_Frag1(NULL),
      m_Frag2(NULL),
      m_MaxCell(50000),
      m_cellSize(1.0),
      m_MCPnts(10000),
      m_knmtcFctr(),
      m_potential(),
      m_Sym(1),
      m_top1(UNDEFINED_TOP),
      m_top2(UNDEFINED_TOP),
      m_mu(0.0) {}
    virtual ~PhaseIntegral() {}

    virtual void initialize(Molecule* Frag1, Molecule* Frag2, Reaction* pReact);

    virtual void integrate(double rxnCrd, vector<double> &wrk) = 0;

    void convolveExcessEnergy(size_t TDOF, vector<double> &cellSOS) const;

  protected:

    Molecule* m_Frag1;
    Molecule* m_Frag2;

    size_t m_MaxCell;
    double m_cellSize;

    size_t m_MCPnts;
    vector<double> m_knmtcFctr; // Configuration kinematic factor.
    vector<double> m_potential; // Torsion configuration potential.
    size_t m_Sym;

    RotationalTop m_top1;
    RotationalTop m_top2;

    double m_mu;
  };

  class NLnrNLnrTops : public PhaseIntegral {

  public:

    NLnrNLnrTops() : PhaseIntegral(), m_nIDOF(5) {}
    virtual ~NLnrNLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk);

  private:

    const int m_nIDOF;
  };

  class NLnrLnrTops : public PhaseIntegral {

  public:
    NLnrLnrTops() : PhaseIntegral(), m_nIDOF(4) {}
    virtual ~NLnrLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk);

  private:

    const int m_nIDOF;
  };

  class NLnrAtmTops : public PhaseIntegral {

  public:
    NLnrAtmTops() : PhaseIntegral(), m_nIDOF(2) {}
    virtual ~NLnrAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk) ;

  private:

    // Hirst potential
    double ptnl(double rxnCrd, double gamma, double phi) {
      //C and A the potential barrier Parameters;
      double tmp = sin(gamma);
      return  m_V0 * tmp * tmp;
    }

    double m_V0;
    const int m_nIDOF;
  };

  class LnrLnrTops : public PhaseIntegral {

  public:
    LnrLnrTops() : PhaseIntegral(), m_nIDOF(3) {}
    virtual ~LnrLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk);

  private:

    const int m_nIDOF;
  };

  class LnrAtmTops : public PhaseIntegral {

  public:
    LnrAtmTops() : PhaseIntegral(), m_nIDOF(1) {}
    virtual ~LnrAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk);

  private:

    const int m_nIDOF;
  };

  class AtmAtmTops : public PhaseIntegral {

  public:

    AtmAtmTops() : PhaseIntegral(), m_nIDOF(0) {}
    virtual ~AtmAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk) {}

  private:

    const int m_nIDOF;
  };

}//namespace

#endif // GUARD_PhaseIntergrals_h