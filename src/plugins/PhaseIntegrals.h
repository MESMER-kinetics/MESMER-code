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
#include "../gStructure.h"

namespace mesmer
{
  class FTSTPotential;

  class PhaseIntegral {

  public:

    PhaseIntegral(size_t nIDOF) :
      m_nIDOF(nIDOF),
      m_Frag1(NULL),
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
      m_mu(0.0) {
    }
    virtual ~PhaseIntegral() {}

    virtual void initialize(Molecule* Frag1, Molecule* Frag2, FTSTPotential* pFTSTPotential, Reaction* pReact);

    virtual void integrate(double rxnCrd, vector<double>& wrk) = 0;

    void convolveExcessEnergy(vector<double>& cellSOS) const;

  protected:

    // Heavyside function integration.
    void HeavysideIntegration(vector<double>& cellSOS) const;

    const size_t m_nIDOF; // Number of transitional modes (excluding external rotation).

    Molecule* m_Frag1;
    Molecule* m_Frag2;

    size_t m_MaxCell;
    double m_cellSize;

    size_t m_MCPnts;            // Number of configurations to be used in Monte-Carlo integration.
    vector<double> m_knmtcFctr; // Configuration kinematic factor.
    vector<double> m_potential; // Torsion configuration potential.
    size_t m_Sym;               // Symmetry number.

    RotationalTop m_top1;
    RotationalTop m_top2;

    FTSTPotential* m_pFTSTPotential;

    double m_mu;
  };

  class NLnrNLnrTops : public PhaseIntegral {

  public:

    NLnrNLnrTops() : PhaseIntegral(5) {}
    virtual ~NLnrNLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk);

  private:

  };

  class NLnrLnrTops : public PhaseIntegral {

  public:
    NLnrLnrTops() : PhaseIntegral(4) {}
    virtual ~NLnrLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk);

  private:

  };

  class NLnrAtmTops : public PhaseIntegral {

  public:
    NLnrAtmTops() : PhaseIntegral(2) {}
    virtual ~NLnrAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk);

  private:

  };

  class LnrLnrTops : public PhaseIntegral {

  public:
    LnrLnrTops() : PhaseIntegral(3) {}
    virtual ~LnrLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk);

  private:

  };

  class LnrAtmTops : public PhaseIntegral {

  public:
    LnrAtmTops() : PhaseIntegral(1) {}
    virtual ~LnrAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk);

  private:

  };

  class AtmAtmTops : public PhaseIntegral {

  public:

    AtmAtmTops() : PhaseIntegral(0) {}
    virtual ~AtmAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk) {}

  private:

  };

  // The following phase integral is taken from Seakins et al JPC 9974, 101 (1997).

  class MethylPlusH_HW : public PhaseIntegral {

  public:
    MethylPlusH_HW() : PhaseIntegral(2) {}
    virtual ~MethylPlusH_HW() {}

    virtual void integrate(double rxnCrd, vector<double>& wrk);

  private:

  };

}//namespace

#endif // GUARD_PhaseIntegrals_h