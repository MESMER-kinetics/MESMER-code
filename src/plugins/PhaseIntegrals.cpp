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
#include "../gStructure.h"
#include "../Reaction.h"
#include "../Sobol.h"
#include "../FTSTPotential.h"
#include "PhaseIntegrals.h"

using namespace std;
using namespace Constants;

namespace mesmer
{

  void PhaseIntegral::initialize(Molecule* Frag1, Molecule* Frag2, FTSTPotential* pFTSTPotential, Reaction* pReact) {
    m_Frag1 = Frag1;
    m_Frag2 = Frag2;

    m_pFTSTPotential = pFTSTPotential;

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

    // Get fragment coordinates and shit to cenre of mass.
    ReadCoordsAndShiftToCoM(m_Frag1, m_m1, m_x1, m_y1, m_z1);
    ReadCoordsAndShiftToCoM(m_Frag2, m_m2, m_x2, m_y2, m_z2);

    PersistPtr pp = pReact->get_PersistentPointer()->XmlMoveTo("me:MCRCMethod");
    m_Sym = pp->XmlReadInteger("me:SymmetryNumber", optional);
    int MCpnts = pp->XmlReadInteger("me:MCPoints", optional);
    if (MCpnts > 0)
      m_MCPnts = size_t(MCpnts);
    else
      cerr << "Number of of Monte-Carlo points for coupled rotors method must be greater than 0." << endl;

  };

  void PhaseIntegral::ReadCoordsAndShiftToCoM(Molecule* Frag, vector<double>& m, vector<double>& x, vector<double>& y, vector<double>& z) {

    // Get fragment coordinates.

    Frag->getStruc().getAtomicMasses(m);
    Frag->getStruc().getXCoords(x);
    Frag->getStruc().getYCoords(y);
    Frag->getStruc().getZCoords(z);

    // Find mass centres and centre fragments.

    double smm(0.0), smx(0.0), smy(0.0), smz(0.0);
    for (size_t i(0); i < m.size(); ++i) {
      smm += m[i];
      smx += m[i] * x[i];
      smy += m[i] * y[i];
      smz += m[i] * z[i];
    }

    smx /= smm; smy /= smm; smz /= smm;
    for (size_t i(0); i < m.size(); ++i) {
      x[i] -= smx;
      y[i] -= smy;
      z[i] -= smz;
    }

  }

  void PhaseIntegral::convolveExcessEnergy(vector<double>& cellSOS) const {

    vector<double> tmpCellSOS(cellSOS);
    vector<double> ene(cellSOS.size(), 0.0);
    getCellEnergies(cellSOS.size(), m_cellSize, ene);
    double pwr = 0.5 * double(m_nIDOF + 3) - 1.0;
    for (size_t j(0); j < ene.size(); ++j) {
      ene[j] = pow(ene[j], pwr);
    }
    FastLaplaceConvolution(tmpCellSOS, ene, cellSOS);

  }

  void PhaseIntegral::convolveExcessEnergy(vector<double>& cellSOS, size_t J) const {

    if (J == 0) {
      for (size_t i(0); i < cellSOS.size(); ++i) {
        cellSOS[i] = 1.0;
      }
    }
    else {
      vector<double> tmpCellSOS(cellSOS);
      vector<double> ene(cellSOS.size(), 0.0);
      getCellEnergies(cellSOS.size(), m_cellSize, ene);
      double pwr = 0.5 * double(m_nIDOF) - 1.0;
      for (size_t j(0); j < ene.size(); ++j) {
        ene[j] = pow(ene[j], pwr);
      }
      FastLaplaceConvolution(tmpCellSOS, ene, cellSOS);
    }

  }

  // Heaviside function integration.
  void PhaseIntegral::HeavisideIntegration(vector<double>& cellSOS) const {
    for (size_t i(0); i < m_MCPnts; ++i) {
      double kfctr = m_knmtcFctr[i];
      double ptnl = m_potential[i];
      size_t ll = size_t(ptnl / m_cellSize);
      for (size_t j(ll); j < cellSOS.size(); ++j) {
        cellSOS[j] += kfctr;
      }
    }
  }

  // Instantaneous moments of inertia.
  void PhaseIntegral::InstMoI(double rxnCrd, double& Ba, double& Bb, double& Bc) const {

    // Rotate fragments.

    // Shift fragments relative to each other along z-axis.
    double m1 = m_Frag1->getStruc().CalcMW();
    double m2 = m_Frag2->getStruc().CalcMW();
    double r1 = m1 * rxnCrd / (m1 + m2);
    double r2 = m2 * rxnCrd / (m1 + m2);

    size_t natoms = m_m1.size() + m_m2.size();
    vector<double> m(natoms, 0.0), x(natoms, 0.0), y(natoms, 0.0), z(natoms, 0.0);
    vector<double> r(3, 0.0);
    size_t ll(0);
    for (size_t i(0); i < m_m1.size(); i++, ll++) {
      r[0] = m_x1[i]; r[1] = m_y1[i]; r[2] = m_z1[i];
      r *= m_rot1;
      m[ll] = m_m1[i];
      x[ll] = r[0];
      y[ll] = r[1];
      z[ll] = r[2] - r2;
    }
    for (size_t i(0); i < m_m2.size(); i++, ll++) {
      r[0] = m_x2[i]; r[1] = m_y2[i]; r[2] = m_z2[i];
      r *= m_rot2;
      m[ll] = m_m2[i];
      x[ll] = r[0];
      y[ll] = r[1];
      z[ll] = r[2] + r1;
    }

    // Calculate moments of interia of the ensembly.

    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
    for (size_t i(0); i < m.size() ; i++)
    {
      sxx += m[i] * x[i] * x[i];
      syy += m[i] * y[i] * y[i];
      szz += m[i] * z[i] * z[i];
      sxy += m[i] * x[i] * y[i];
      sxz += m[i] * x[i] * z[i];
      syz += m[i] * y[i] * z[i];
    }

    dMatrix MI(3);
    MI[0][0] = syy + szz;
    MI[1][1] = sxx + szz;
    MI[2][2] = sxx + syy;
    MI[0][1] = MI[1][0] = -sxy;
    MI[0][2] = MI[2][0] = -sxz;
    MI[1][2] = MI[2][1] = -syz;

    vector<double> PrincipalMI(3, 0.0); // amuAng2
    MI.diagonalize(&PrincipalMI[0]);

    Ba = conMntInt2RotCnt/PrincipalMI[2];
    Bb = conMntInt2RotCnt/PrincipalMI[1];
    Bc = conMntInt2RotCnt/PrincipalMI[0];

  }

  void LnrAtmTops::integrate(double rxnCrd, vector<double>& cellSOS) {

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
      angles[0] *= M_PI;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sin(angles[0]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == LINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    double cnt = M_PI * RotCnt / double(2.0 * m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void LnrAtmTops::integrate(double rxnCrd, vector<double>& cellSOS, size_t J) {

    // Instantiate a random vector generator.
    Sobol sobol;

    // Configuration loop.

    double B0 = conMntInt2RotCnt / (m_mu * rxnCrd * rxnCrd);
    double Ba(0.0), Bb(0.0), Bc(0.0);
    m_knmtcFctr.resize(m_MCPnts, 0.0);
    m_potential.resize(m_MCPnts, 0.0);
    vector<double> angles(m_nIDOF, 0.0);
    vector<double> tmp(m_nIDOF + 2, 0.0);
    long long seed(17);
    for (size_t i(0); i < m_MCPnts; ++i) {

      // Select angular coordinates.
      vector<double> angles(m_nIDOF, 0.0);
      sobol.sobol(tmp.size(), &seed, tmp);
      angles[0] = M_PI * tmp[0];
      double nu = 2.0 * M_PI * tmp[1];
      double gamma = 2.0 * tmp[2] - 1.0;

      // Calculate the instantaneous moments of inertia.
      InstMoI(rxnCrd, Ba, Bb, Bc);

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sqrt(Ba * Bb * Bc) * sin(angles[0]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);

      // Add instantaneous rotational energy.
      m_potential[i] += B0 * J * J;
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = 1.0 / B0;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == LINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    double cnt = 2.0 * J * J * RotCnt / double(m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void NLnrAtmTops::integrate(double rxnCrd, vector<double>& cellSOS) {

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
      angles[0] *= M_PI;
      angles[1] *= 2.0 * M_PI;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sin(angles[0]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = 2.0 * M_PI * RotCnt / double(3.0 * m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void NLnrAtmTops::integrate(double rxnCrd, vector<double>& cellSOS, size_t J) {

    if (J == 0)
      return;

    // Instantiate a random vector generator.
    Sobol sobol;

    // Configuration loop.
    double B0 = conMntInt2RotCnt / (m_mu * rxnCrd * rxnCrd);
    double Ba(0.0), Bb(0.0), Bc(0.0);

    m_rot2[0][0] = m_rot2[1][1] = m_rot2[2][2] = 1.0;

    m_knmtcFctr.resize(m_MCPnts, 0.0);
    m_potential.resize(m_MCPnts, 0.0);
    vector<double> angles(m_nIDOF, 0.0);
    vector<double> tmp(m_nIDOF + 2, 0.0);
    long long seed(17);
    for (size_t i(0); i < m_MCPnts; ++i) {

      // Select angular coordinates.
      sobol.sobol(tmp.size(), &seed, tmp);
      angles[0] = 2.0 * M_PI * tmp[0];
      angles[1] = M_PI * tmp[1];
      double nu = 2.0 * M_PI * tmp[2];
      double gamma = 2.0 * tmp[3] - 1.0;

      dMatrix rotZ(3, 0.0);
      rotZ[0][0] = rotZ[1][1] = cos(angles[0]);
      rotZ[2][2] = 1.0;
      rotZ[0][1] = sin(angles[0]);
      rotZ[1][0] = -rotZ[0][1];

      dMatrix rotY(3, 0.0);
      rotY[0][0] = rotY[2][2] = cos(angles[1]);
      rotY[1][1] = 1.0;
      rotY[0][2] = -sin(angles[1]);
      rotY[2][0] = -rotY[0][2];

      m_rot1 = rotY * rotZ;

      // Calculate the instantaneous moments of inertia.
      InstMoI(rxnCrd, Ba, Bb, Bc);

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sqrt(Ba * Bb * Bc) * sin(angles[0]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);

      // Add instantaneous rotational energy.
      double g2 = gamma * gamma;
      double S1 = sin(nu);
      double S2 = S1 * S1;
      m_potential[i] += J * J * ((Ba*S2 + Bb*(1-S2)) * (1.0 - g2) + Bc*g2);
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = 1.0 / B0;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = 2.0 * M_PI * J * J * RotCnt / double(m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void LnrLnrTops::integrate(double rxnCrd, vector<double>& cellSOS) {

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
      angles[0] *= M_PI;
      angles[1] *= M_PI;
      angles[2] *= 2.0 * M_PI;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sin(angles[0]) * sin(angles[1]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    vector<double> MntsInt;
    m_Frag1->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    m_Frag2->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    double cnt = M_PI * M_PI * RotCnt / double(8.0 * m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void LnrLnrTops::integrate(double rxnCrd, vector<double>& cellSOS, size_t J) {

    if (J == 0)
      return;

    // Instantiate a random vector generator.
    Sobol sobol;

    // Configuration loop.
    double B0 = conMntInt2RotCnt / (m_mu * rxnCrd * rxnCrd);
    m_knmtcFctr.resize(m_MCPnts, 0.0);
    m_potential.resize(m_MCPnts, 0.0);
    vector<double> angles(m_nIDOF, 0.0);
    vector<double> tmp(m_nIDOF + 2, 0.0);
    long long seed(17);
    for (size_t i(0); i < m_MCPnts; ++i) {

      // Select angular coordinates.
      sobol.sobol(tmp.size(), &seed, tmp);
      angles[0] = M_PI * tmp[0];
      angles[1] = M_PI * tmp[1];
      angles[2] = 2.0 * M_PI * tmp[2];
      double nu = 2.0 * M_PI * tmp[3];
      double gamma = 2.0 * tmp[4] - 1.0;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sqrt(B0 * B0 * B0) * sin(angles[0]) * sin(angles[1]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);

      // Add instantaneous rotational energy.
      m_potential[i] += B0 * J * J;
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = 1.0 / B0;
    vector<double> MntsInt;
    m_Frag1->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    m_Frag2->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    double cnt = 2.0 * M_PI * J * J * RotCnt / double(m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void NLnrLnrTops::integrate(double rxnCrd, vector<double>& cellSOS) {

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
      angles[0] *= M_PI;
      angles[1] *= 2.0 * M_PI;
      angles[2] *= M_PI;
      angles[3] *= 2.0 * M_PI;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sin(angles[0]) * sin(angles[2]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    top = (m_top1 == LINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    double cnt = 2.0 * M_PI * M_PI * RotCnt / double(15.0 * m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }
  }

  void NLnrLnrTops::integrate(double rxnCrd, vector<double>& cellSOS, size_t J) {

    if (J == 0)
      return;
    // Instantiate a random vector generator.
    Sobol sobol;

    // Configuration loop.
    double B0 = conMntInt2RotCnt / (m_mu * rxnCrd * rxnCrd);
    m_knmtcFctr.resize(m_MCPnts, 0.0);
    m_potential.resize(m_MCPnts, 0.0);
    vector<double> angles(m_nIDOF, 0.0);
    vector<double> tmp(m_nIDOF + 2, 0.0);
    long long seed(17);
    for (size_t i(0); i < m_MCPnts; ++i) {

      // Select angular coordinates.
      sobol.sobol(tmp.size(), &seed, tmp);
      angles[0] = M_PI * tmp[0];
      angles[1] = 2.0 * M_PI * tmp[1];
      angles[2] = M_PI * tmp[2];
      angles[3] = 2.0 * M_PI * tmp[3];
      double nu = 2.0 * M_PI * tmp[4];
      double gamma = 2.0 * tmp[5] - 1.0;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sqrt(B0 * B0 * B0) * sin(angles[0]) * sin(angles[2]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);

      // Add instantaneous rotational energy.
      m_potential[i] += B0 * J * J;
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = 1.0 / B0;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    top = (m_top1 == LINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= MntsInt[0];
    double cnt = M_PI * M_PI * J * J * RotCnt / double(m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }
  }

  void NLnrNLnrTops::integrate(double rxnCrd, vector<double>& cellSOS) {

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
      angles[0] *= M_PI;
      angles[1] *= 2.0 * M_PI;
      angles[2] *= M_PI;
      angles[3] *= 2.0 * M_PI;
      angles[4] *= 2.0 * M_PI;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sin(angles[0]) * sin(angles[2]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    vector<double> MntsInt;
    m_Frag1->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    m_Frag2->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = M_PI * M_PI * M_PI * RotCnt / double(24.0 * m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }
  }

  void NLnrNLnrTops::integrate(double rxnCrd, vector<double>& cellSOS, size_t J) {

    if (J == 0)
      return;
    // Instantiate a random vector generator.
    Sobol sobol;

    // Configuration loop.
    double B0 = conMntInt2RotCnt / (m_mu * rxnCrd * rxnCrd);
    m_knmtcFctr.resize(m_MCPnts, 0.0);
    m_potential.resize(m_MCPnts, 0.0);
    vector<double> angles(m_nIDOF, 0.0);
    vector<double> tmp(m_nIDOF + 2, 0.0);
    long long seed(17);
    for (size_t i(0); i < m_MCPnts; ++i) {

      // Select angular coordinates.
      sobol.sobol(tmp.size(), &seed, tmp);
      angles[0] = M_PI * tmp[0];
      angles[1] = 2.0 * M_PI * tmp[1];
      angles[2] = M_PI * tmp[2];
      angles[3] = 2.0 * M_PI * tmp[3];
      angles[4] = 2.0 * M_PI * tmp[4];
      double nu = 2.0 * M_PI * tmp[5];
      double gamma = 2.0 * tmp[6] - 1.0;

      // Calculate the determinant of the Wilson G Matrix.
      m_knmtcFctr[i] = sqrt(B0 * B0 * B0) * sin(angles[0]) * sin(angles[2]);

      // Calculate potential energy.
      m_potential[i] = m_pFTSTPotential->HinderingPotential(rxnCrd, angles);

      // Add instantaneous rotational energy.
      m_potential[i] += B0 * J * J;
    }

    // Heaviside function integration.
    HeavisideIntegration(cellSOS);

    // Conversion and symmetry number factor.
    double RotCnt = 1.0 / B0;
    vector<double> MntsInt;
    m_Frag1->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    m_Frag2->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = 4.0 * M_PI * M_PI * RotCnt / double(3.0 * m_MCPnts * m_Sym);
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] *= cnt;
    }

  }

  void MethylPlusH_HW::integrate(double rxnCrd, vector<double>& cellSOS) {

    // Conversion and symmetry number factor.
    double RotCnt = m_mu * rxnCrd * rxnCrd / conMntInt2RotCnt;
    vector<double> MntsInt;
    Molecule* top = (m_top1 == NONLINEAR) ? m_Frag1 : m_Frag2;
    top->getDOS().get_rotConsts(MntsInt);
    RotCnt /= sqrt(MntsInt[0] * MntsInt[1] * MntsInt[2]);
    double cnt = 2.0 * RotCnt / (3.0 * m_Sym);

    // The following is a specific test for the CH3 + H system taken from JPC 101, 9974 (1997).

    // Hirst potential
    const double C = getConvertedEnergy("kcal/mol", 164.0);
    const double A = 0.43;
    const double V0 = C * exp(-A * rxnCrd * rxnCrd);

    vector<double> ene(cellSOS.size(), 0.0);
    getCellEnergies(cellSOS.size(), m_cellSize, ene);

    cnt *= 2.0;
    for (size_t j(0); j < cellSOS.size(); ++j) {
      cellSOS[j] = (ene[j] < V0) ? cnt * (1 - (sqrt(1.0 - ene[j] / V0))) : cnt;
    }

  }

}
