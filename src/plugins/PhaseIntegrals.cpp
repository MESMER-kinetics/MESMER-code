//-------------------------------------------------------------------------------------------
//
// PhaseIntergrals.cpp
//
// Author: Struan Robertson
// Date:   April/2019
//
// This file contains the definitions of the various phase integrals needed for
// microcanonical (NOT J resolved) flexible transition state theory.
//
//-------------------------------------------------------------------------------------------

#include "../gDensityOfStates.h"
#include "../Reaction.h"
#include "../Sobol.h"
#include "../FTSTPotential.h"
#include "PhaseIntegrals.h"

using namespace std;
using namespace Constants;

namespace mesmer
{

  void PhaseIntegral::initialize(Molecule* Frag1, Molecule* Frag2, FTSTPotential *pFTSTPotential, Reaction* pReact) {
    m_Frag1 = Frag1;
    m_Frag2 = Frag2;

    m_pFTSTPotential = pFTSTPotential;

    m_MaxCell  = pReact->getEnv().MaxCell;;
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

  void PhaseIntegral::convolveExcessEnergy(size_t TDOF, vector<double> &cellSOS) const {

    vector<double> tmpCellSOS(cellSOS);
    vector<double> ene(cellSOS.size(), 0.0);
    getCellEnergies(cellSOS.size(), m_cellSize, ene);
    double pwr = 0.5*double(TDOF) - 1.0;
    for (size_t j(0); j < ene.size(); ++j) {
      ene[j] = pow(ene[j], pwr);
    }
    FastLaplaceConvolution(tmpCellSOS, ene, cellSOS);

  }

  void NLnrNLnrTops::integrate(double rxnCrd, vector<double> &cellSOS) {

    const size_t TDOF = m_nIDOF + 3;
  
    // Convolve with remaining energy contributions.
    convolveExcessEnergy(TDOF, cellSOS);
  }

  void NLnrLnrTops::integrate(double rxnCrd, vector<double> &cellSOS) {

    const size_t TDOF = m_nIDOF + 3;

    // Convolve with remaining energy contributions.
    convolveExcessEnergy(TDOF, cellSOS);
  }

  void NLnrAtmTops::integrate(double rxnCrd, vector<double> &cellSOS) {

    m_pFTSTPotential->RxnCrdInitialize(rxnCrd);

    Molecule *top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;

    const double cellSize = m_cellSize;
    const size_t TDOF = m_nIDOF + 3;

    // Instantiate a random vector generator.
    Sobol sobol;

    // Configuration loop.
    long long seed(17);
    m_knmtcFctr.resize(m_MCPnts, 0.0);
    m_potential.resize(m_MCPnts, 0.0);
    for (size_t i(0); i < m_MCPnts; ++i) {

      // Select angular coordinates.
      vector<double> angles(m_nIDOF, 0.0);
      sobol.sobol(angles.size(), &seed, angles);
      angles[0] *= M_PI ;
      angles[1] *= 2.0 *  M_PI ;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sin(angles[0]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles) ;
    }

    // Heavyside function integration.
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
    double RotCnt = orbitalInertia / sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = 2.0*M_PI*RotCnt / double(3.0*m_MCPnts*m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

    // Convolve with remaining energy contributions.
    convolveExcessEnergy(TDOF, cellSOS);

  }

  void LnrLnrTops::integrate(double rxnCrd, vector<double> &cellSOS) {

    const size_t TDOF = m_nIDOF + 3;

    // Convolve with remaining energy contributions.
    convolveExcessEnergy(TDOF, cellSOS);
  }

  void LnrAtmTops::integrate(double rxnCrd, vector<double> &cellSOS) {

    const size_t TDOF = m_nIDOF + 3;

    // Convolve with remaining energy contributions.
    convolveExcessEnergy(TDOF, cellSOS);
  }

  void MethylPlusH_HW::integrate(double rxnCrd, vector<double>& cellSOS) {

    Molecule* top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;

    const double cellSize = m_cellSize;
    const size_t TDOF = m_nIDOF + 3;

    // Conversion and symmetry number factor.
    vector<double> MntsInt;
    top->getDOS().get_rotConsts(MntsInt);
    double orbitalInertia = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    double RotCnt = orbitalInertia / sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = 2.0 * RotCnt / (3.0 * m_Sym);

    // The following is a specific test for the CH3 + H system taken from JPC 101, 9974 (1997).

    // Hirst potential
    const double C = getConvertedEnergy("kcal/mol", 164.0);
    const double A = 0.43;
    const double V0 = C * exp(-A * rxnCrd * rxnCrd);

    vector<double> ene(cellSOS.size(), 0.0);
    getCellEnergies(cellSOS.size(), cellSize, ene);

    cnt *= 2.0;
    for (size_t j(0); j < cellSOS.size() ; ++j) {
      cellSOS[j] = (ene[j] < V0) ? cnt * ( 1 - (sqrt(1.0 - ene[j]/V0))) : cnt ;
    }

    // Convolve with remaining energy contributions.
    convolveExcessEnergy(TDOF, cellSOS);
  }

}
