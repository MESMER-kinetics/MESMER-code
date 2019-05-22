//-------------------------------------------------------------------------------------------
// PhaseIntegrals.h
//
// Author: Struan Robertson
// Date:   Feb/2019
//
// This file contains the definitions of the various phase integrals needed for
// microcanonical (NOT J resolved) flexible transition state theory.
//-------------------------------------------------------------------------------------------

#ifndef GUARD_PhaseIntegrals_h
#define GUARD_PhaseIntegrals_h

#include "../Molecule.h"
#include "gStructure.h"

namespace mesmer
{
  class FTSTPotential;

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
      m_pFTSTPotential(NULL),
      m_mu(0.0) {}
    virtual ~PhaseIntegral() {}

    virtual void initialize(Molecule* Frag1, Molecule* Frag2, FTSTPotential *pFTSTPotential, Reaction* pReact);

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

    FTSTPotential *m_pFTSTPotential;

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

#endif // GUARD_PhaseIntegrals_h