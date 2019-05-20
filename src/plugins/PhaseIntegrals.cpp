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

#include "PhaseIntegrals.h"
#include "../gDensityOfStates.h"
#include "../Reaction.h"
#include "../Sobol.h"
#include "../FTSTPotential.h"

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

    const size_t MaximumCell = cellSOS.size();
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

    //// Restore rotatable bond IDs.
    //gs.setRotBondID(SavedRotBondIDs);

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

#ifdef _DEBUG

    // The following is a specific test for the CH3 + H system taken from JPC 101, 9974 (1997).

    vector<double> ene(cellSOS.size(), 0.0);
    vector<double> wrk = ene;
    vector<double> tmpCellSOS = ene;
    getCellEnergies(MaximumCell, cellSize, ene);
    size_t jj = size_t(m_V0);
    for (size_t j(0); j < jj; ++j) {
      wrk[j] = 1.0 / (2.0*sqrt(m_V0*m_V0 - m_V0 * ene[j]));
    }

    // Convolve with remaining energy contributions.
    double pwr = 0.5*double(TDOF);
    cnt = 4.0*RotCnt / (15.0*m_Sym);
    for (size_t j(0); j < ene.size(); ++j) {
      ene[j] += 2.0*cnt * pow(ene[j], pwr);
    }
    FastLaplaceConvolution(wrk, ene, tmpCellSOS);


    ctest << endl;
    ctest << "  115 " << formatFloat(wrk[154], 6, 14) << formatFloat(tmpCellSOS[154], 6, 14) << endl;
    ctest << "  248 " << formatFloat(wrk[247], 6, 14) << formatFloat(tmpCellSOS[247], 6, 14) << endl;
    ctest << "  413 " << formatFloat(wrk[412], 6, 14) << formatFloat(tmpCellSOS[412], 6, 14) << endl;
    ctest << "  624 " << formatFloat(wrk[623], 6, 14) << formatFloat(tmpCellSOS[623], 6, 14) << endl;
    ctest << " 1040 " << formatFloat(wrk[1039], 6, 14) << formatFloat(tmpCellSOS[1039], 6, 14) << endl;
    ctest << " 1248 " << formatFloat(wrk[1247], 6, 14) << formatFloat(tmpCellSOS[1247], 6, 14) << endl;
    ctest << " 2007 " << formatFloat(wrk[2006], 6, 14) << formatFloat(tmpCellSOS[2006], 6, 14) << endl;
    ctest << " 2080 " << formatFloat(wrk[2079], 6, 14) << formatFloat(tmpCellSOS[2079], 6, 14) << endl;
    ctest << " 2600 " << formatFloat(wrk[2599], 6, 14) << formatFloat(tmpCellSOS[2599], 6, 14) << endl;
    ctest << " 3120 " << formatFloat(wrk[3119], 6, 14) << formatFloat(tmpCellSOS[3119], 6, 14) << endl;
    ctest << " 3419 " << formatFloat(wrk[3418], 6, 14) << formatFloat(tmpCellSOS[3418], 6, 14) << endl;
    ctest << " 4015 " << formatFloat(wrk[4014], 6, 14) << formatFloat(tmpCellSOS[4014], 6, 14) << endl;
    ctest << " 4445 " << formatFloat(wrk[4444], 6, 14) << formatFloat(tmpCellSOS[4444], 6, 14) << endl;
    ctest << " 5554 " << formatFloat(wrk[5553], 6, 14) << formatFloat(tmpCellSOS[5553], 6, 14) << endl;
    ctest << endl;

    // wrk = tmpCellDOS;

#endif

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

}
