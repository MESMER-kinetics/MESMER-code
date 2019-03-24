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
#include "../Sobol.h"

namespace mesmer
{

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

    virtual void initialize(Molecule* Frag1, Molecule* Frag2, Reaction* pReact) {
      m_Frag1 = Frag1;
      m_Frag2 = Frag2;

      m_MaxCell = pReact->getEnv().MaxCell;;
      m_cellSize = pReact->getEnv().CellSize;;

      // Calculate the reduced mass.
      m_Frag1->getStruc().ReadStructure();
      m_Frag2->getStruc().ReadStructure();
      double m1 = m_Frag1->getStruc().CalcMW();
      double m2 = m_Frag2->getStruc().CalcMW();
      m_mu = m1 * m2 / (m1 + m2);

      m_top1 = m_Frag1->getDOS().get_rotType();
      m_top2 = m_Frag2->getDOS().get_rotType();

      PersistPtr pp = pReact->get_PersistentPointer()->XmlMoveTo("me:MCRCMethod");
      m_Sym = pp->XmlReadInteger("me:SymmetryNumber", optional);
      int MCpnts = pp->XmlReadInteger("me:MCPoints", optional);
      if (MCpnts > 0)
        m_MCPnts = size_t(MCpnts);
      else
        cerr << "Number of of Monte-Carlo points for coupled rotors method must be greater than 0." << endl;

    };

    virtual void integrate(double rxnCrd, vector<double> &wrk) = 0;

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

    virtual void integrate(double rxnCrd, vector<double> &wrk) {}

  private:

    const int m_nIDOF;
  };

  class NLnrLnrTops : public PhaseIntegral {

  public:
    NLnrLnrTops() : PhaseIntegral(), m_nIDOF(4) {}
    virtual ~NLnrLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk) {}

  private:

    const int m_nIDOF;
  };

  class NLnrAtmTops : public PhaseIntegral {

  public:
    NLnrAtmTops() : PhaseIntegral(), m_nIDOF(2) {}
    virtual ~NLnrAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk) {

      Molecule *top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;

      const size_t MaximumCell = wrk.size();
      const double cellSize = m_cellSize;
      const size_t TDOF = m_nIDOF + 3;

      // Instantiate a random vector generator.
      Sobol sobol;

      // Configuration loop.
      long long seed(1);
      m_knmtcFctr.resize(m_MCPnts, 0.0);
      m_potential.resize(m_MCPnts, 0.0);
      for (size_t i(0); i < m_MCPnts; ++i) {

        // Select angular coordinates.
        vector<double> angles(m_nIDOF, 0.0);
        sobol.sobol(angles.size(), &seed, angles);
        double gamma = M_PI * angles[0];
        double phi   = 2.0 *  M_PI * angles[1];

        // Calculate the determinant of the Wilson G Matrix.
       m_knmtcFctr[i] = sin(gamma);

        // Calculate potential energy.
       m_potential[i] = ptnl(rxnCrd, gamma, phi); 
      }

      //// Restore rotatable bond IDs.
      //gs.setRotBondID(SavedRotBondIDs);

      // Heavyside function integration.
      vector<double> cellSOS(MaximumCell, 0.0);
      for (size_t i(0); i < m_MCPnts; ++i) {
        double kfctr = m_knmtcFctr[i];
        double ptnl = m_potential[i];
        size_t ll = size_t(ptnl / cellSize);
        for (size_t j(ll); j < cellSOS.size(); ++j) {
          cellSOS[j] += kfctr;
        }
      }

      // Conversion and symmetry number factor.
      vector<double> MntsInt;
      top->getDOS().get_rotConsts(MntsInt);
      double orbitalInertia = m_mu * rxnCrd*rxnCrd / conMntInt2RotCnt;
      double cnt = orbitalInertia / sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
      cnt *= 2.0 / double(3.0*m_MCPnts*m_Sym);
      for (size_t j(0); j < cellSOS.size(); ++j) {
        cellSOS[j] *= cnt;
      }

      // Convolve with remaining energy contributions.
      vector<double> tmpCellDOS(cellSOS.size(), 0.0);
      getCellEnergies(MaximumCell, cellSize, tmpCellDOS);
      double pwr = 0.5*double(TDOF) - 1.0;
      for (size_t j(0); j < tmpCellDOS.size(); ++j) {
        tmpCellDOS[j] += pow(tmpCellDOS[j], pwr);
      }
      FastLaplaceConvolution(cellSOS, tmpCellDOS, wrk);

    }

  private:

    // Hirst potential
    double ptnl(double rxnCrd, double gamma, double phi) {
      //C and A the potential barrier Parameters;
      double ant = getConvertedEnergy("kcal/mol", 164.0);
      double bnt = 0.43;
      double V0  = ant*exp(-bnt * rxnCrd*rxnCrd);
      double tmp = sin(gamma);
      return  V0 * tmp * tmp;
    }

    const int m_nIDOF;
  };

  class LnrLnrTops : public PhaseIntegral {

  public:
    LnrLnrTops() : PhaseIntegral(), m_nIDOF(3) {}
    virtual ~LnrLnrTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk) {}

  private:

    const int m_nIDOF;
  };

  class LnrAtmTops : public PhaseIntegral {

  public:
    LnrAtmTops() : PhaseIntegral(), m_nIDOF(1) {}
    virtual ~LnrAtmTops() {}

    virtual void integrate(double rxnCrd, vector<double> &wrk) {}

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