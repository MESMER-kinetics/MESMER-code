//-------------------------------------------------------------------------------------------
//
// ExchangeReaction.cpp
//
// Author: Struan Robertson
// Date:   Aug/2019
//
// This file contains the implementation of the ExchangeReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "ExchangeReaction.h"
#include "gStructure.h"
#include "gWellProperties.h"
#include "Fragmentation.h"
#include "ParseForPlugin.h"

using namespace Constants;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  // Read the Molecular data from input stream.
  bool ExchangeReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if (!ppReactantList)
      ppReactantList = ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1 = ppReactantList->XmlMoveTo("reactant");      // Read reactant details.
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if (!pMol1) {
      cerr << "Cannot find 1st reactant molecule definition for association reaction " << getName() << ".";
      return false;
    }
    PersistPtr ppReactant2 = ppReactant1->XmlMoveTo("reactant");
    Molecule* pMol2 = GetMolRef(ppReactant2);
    if (!pMol2)
    {
      cerr << "Cannot find 2nd reactant molecule definition for association reaction " << getName() << ".";
      return false;
    }

    // if deficientReactantLocation=true, then pMol1 (the first rct
    // in the XML input) is the deficient reactant (m_rct1)

    Molecule* tmp_rct1 = pMol1;
    Molecule* tmp_rct2 = pMol2;

    if (deficientReactantLocation) {
      m_rct1 = tmp_rct1;
      m_rct2 = tmp_rct2;
    }
    else {
      m_rct1 = tmp_rct2;
      m_rct2 = tmp_rct1;
    }

    if (!m_rct1) {
      cerr << "the deficient reactant in the association reaction is undefined" << endl;
      return false;
    }
    if (!m_rct2) {
      cerr << "the excess reactant in the association reaction is undefined" << endl;
      return false;
    }

    PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
    if (!ppProductList)
      ppProductList = ppReac; //Be forgiving; we can get by without a productList element

    PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");     // Read product details. Save them as type Molecule
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);
      if (pMol1) {
        m_pdt1 = pMol1;
      }
      else {
        cerr << "Exchange reaction" << getName() << " has no products defined.";
      }

      PersistPtr ppProduct2 = ppProduct1->XmlMoveTo("product");
      if (ppProduct2) {
        pMol2 = GetMolRef(ppProduct2);
        if (pMol2) {
          m_pdt2 = pMol2;
        }
        else {
          cerr << "Exchange reaction " << getName() << " has only one product defined.";
        }
      }
    }

    // Determine fragmentation model.
    PersistPtr ppDstbn = ppReac->XmlMoveTo("me:FragmentDist");
    m_fragDist = ParseForPlugin<FragDist>(this, "me:FragmentDist", ppReac);
    if (!m_fragDist)
      return false;

    m_fragDist->ReadParameters(ppDstbn, this->getName());

    // Read heat of reaction and rate parameters.
    if (!ReadRateCoeffParameters(ppReac))
      return false;

    // Check the transistion state is defined. 
    if (!m_TransitionState)
      return false;

    return true;
  }

  double ExchangeReaction::calcEquilibriumConstant() {   // Calculate reaction equilibrium constant.
    // equilibrium constant:
    double Keq(0.0);
    const double beta = getEnv().beta;

    vector<double> cellEne;
    getCellEnergies(m_rct1->getEnv().MaxCell, m_rct1->getEnv().CellSize, cellEne);

    // Rovibronic partition function for products/reactants.
    vector<double> tmpDOS;
    m_rct1->getDOS().getCellDensityOfStates(tmpDOS);
    double Qrcts = canonicalPartitionFunction(tmpDOS, cellEne, beta);
    tmpDOS.clear();
    m_rct2->getDOS().getCellDensityOfStates(tmpDOS);
    Qrcts *= canonicalPartitionFunction(tmpDOS, cellEne, beta);

    m_pdt1->getDOS().getCellDensityOfStates(tmpDOS);
    double Qpdts = canonicalPartitionFunction(tmpDOS, cellEne, beta);
    tmpDOS.clear();
    m_pdt2->getDOS().getCellDensityOfStates(tmpDOS);
    Qpdts *= canonicalPartitionFunction(tmpDOS, cellEne, beta);

    // Rovibronic partition function for reactants/products multiplied by translation contribution.
    Qrcts *= translationalContribution(m_rct1->getStruc().getMass(), m_rct2->getStruc().getMass(), beta);
    Qpdts *= translationalContribution(m_pdt1->getStruc().getMass(), m_pdt2->getStruc().getMass(), beta);

    Keq = Qpdts / Qrcts;

    // Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells.
    const double HeatOfReaction = getHeatOfReaction();
    Keq *= exp(-beta * HeatOfReaction);

    return Keq;
  }

  // Is reaction equilibrating and therefore contributes
  // to the calculation of equilibrium fractions.
  bool ExchangeReaction::isEquilibratingReaction(double& Keq, Molecule** rct, Molecule** pdt) {

    Keq = calcEquilibriumConstant();

    *rct = m_rct1;
    *pdt = m_pdt1;

    return true;
  }

  // Add exchange reaction terms to collision matrix.
  void ExchangeReaction::AddReactionTerms(qdMatrix* CollOptr, molMapType& isomermap, const double rMeanOmega)
  {
    // Locate isomers in system matrix.
    const size_t rctLocation = (*m_sourceMap)[m_rct1];
    const size_t pdtLocation = isomermap[m_pdt1];

    // Get densities of states for detailed balance.
    vector<double> pdtEne;
    vector<double> pdtBltz;
    m_pdt1->getDOS().getGrainEnergies(pdtEne);
    m_pdt1->getColl().normalizedGrnBoltzmannDistribution(pdtBltz);

    // Need to know the number of grouped grains in product.
    const size_t pShiftedGrains(m_pdt1->getColl().reservoirShift());

    const size_t pColloptrsize = m_pdt1->getColl().get_colloptrsize();

    if (pColloptrsize != pdtBltz.size()) {
      string msg = __FUNCTION__ + string(": Error: collison operator and distribution have different sizes.\n");
      throw(runtime_error(msg));
    }

    // Get TS Boltzmann distribution for partition of rate.
    vector<double> TStPopFrac; // Population fraction of the transition state.
    m_TransitionState->getColl().set_colloptrsize(pColloptrsize);
    m_TransitionState->getColl().normalizedGrnBoltzmannDistribution(TStPopFrac);

    // Calculate the weighted sum of the energy partition distribution.

    // 1. First need to find the isoenergetic grains for the product and TS.

    calcEffGrnThresholds();
    size_t nRvrThresh = get_EffGrnRvsThreshold();

    // 2. Form the weighted sum of distributions.
    m_fragDist->initialize(this);
    vector<double> totalFragDist(pdtBltz.size(), 0.0);
    for (size_t i(nRvrThresh), idx(0); i < pdtBltz.size(); ++i, ++idx) {

      vector<double> fragDist(pdtBltz.size(), 0.0);
      m_fragDist->calculate(pdtEne[i], fragDist);

      double weight = TStPopFrac[idx];
      for (size_t j(0); j < pdtBltz.size(); ++j) {
        totalFragDist[j] += weight * fragDist[j];
      }
    }

    // 3. Re-normalize the fragmentation distribution.
    normalizeDist(totalFragDist);

    // Get equilibrium constant.
    const qd_real Keq = qd_real(calcEquilibriumConstant());

    // Call high pressure rate method because exchange reaction rate coefficients
    // do not depend on pressure. The extra logic is to control output. 
    bool writeStatus = getFlags().testRateConstantEnabled;
    getFlags().testRateConstantEnabled = false;
    HighPresRateCoeffs(NULL);
    getFlags().testRateConstantEnabled = writeStatus;

    qd_real AssocRateCoeff(m_MtxGrnKf[0]), qdMeanOmega(rMeanOmega);

    qd_real diag(qdMeanOmega * AssocRateCoeff / Keq);
    qd_real offDiag(qdMeanOmega * AssocRateCoeff / sqrt(Keq));

    // Note: reverseThreshE will always be greater than pShiftedGrains here.

    const size_t jj = rctLocation;
    (*CollOptr)[jj][jj] -= qdMeanOmega * AssocRateCoeff;
    for (size_t i(0), j(0); i < pColloptrsize; ++i, ++j) {
      size_t ii(pdtLocation + i);
      qd_real dist(totalFragDist[i]), bltz(pdtBltz[i]);
      (*CollOptr)[ii][ii] -= diag * dist / bltz;          // Loss of the adduct to the source
      (*CollOptr)[jj][ii] = offDiag * dist / sqrt(bltz);  // Reactive gain of the source
      (*CollOptr)[ii][jj] = (*CollOptr)[jj][ii];         // Reactive gain (symmetrization)
    }

    // Return any resources used.
    m_fragDist->clear();

  }

  // Calculate high pressure rate coefficients at current T.
  void ExchangeReaction::HighPresRateCoeffs(vector<double>* pCoeffs) {

    double k_forward(0.0);
    const double beta = getEnv().beta;

    string mrcID = string(m_pMicroRateCalculator->getID());

    if (mrcID == "CanonicalRateCoefficient") {
      m_pMicroRateCalculator->calculateMicroCnlFlux(this);
      k_forward = m_CellFlux[0];
    }
    else {
      throw(std::runtime_error(string("The irreversible exchange reaction " + getName() + " is defined without either a transition state or Arrhenius parameters.\n")));
    }

    // Save high pressure rate coefficient for use in the construction of the collision operator.
    m_MtxGrnKf.clear();
    m_MtxGrnKf.push_back(k_forward * m_ERConc);

    if (pCoeffs) {
      pCoeffs->push_back(k_forward);
      const double Keq = calcEquilibriumConstant();
      if (!IsNan(Keq)) {
        pCoeffs->push_back(k_forward / Keq);
        pCoeffs->push_back(Keq);
      }
    }
    else if (getFlags().testRateConstantEnabled) {
      const double temperature = 1. / (boltzmann_RCpK * beta);
      stest << endl << "Canonical pseudo first order rate constant of irreversible reaction "
        << getName() << " = " << k_forward * m_ERConc << " s-1 (" << temperature << " K)" << endl;
      stest << "Canonical bimolecular rate constant of irreversible reaction "
        << getName() << " = " << k_forward << " cm^3/molec/s (" << temperature << " K)" << endl;
    }
  }

  const int ExchangeReaction::get_rctsGrnZPE() {
    double grnZpe = (m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin) / getEnv().GrainSize; //convert to grains
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void ExchangeReaction::calcEffGrnThresholds(void) {

    const double rctsEne = m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe();
    const double pdtsEne = m_pdt1->getDOS().get_zpe() + m_pdt2->getDOS().get_zpe();
    const double tsEne = m_TransitionState->getDOS().get_zpe();

    if (tsEne < rctsEne) {
      string msg = "Reaction " + getName() + " cannot have a negative threshold.";
      throw(std::runtime_error(msg));
    }
    else {
      set_EffGrnFwdThreshold(int(tsEne - rctsEne) / getEnv().GrainSize);
      set_EffGrnRvsThreshold(int(tsEne - pdtsEne) / getEnv().GrainSize);
    }
  }

  void ExchangeReaction::calcFluxFirstNonZeroIdx(void) {
    double thresh = get_ThresholdEnergy();
    if (thresh < 0.0) { m_GrnFluxFirstNonZeroIdx = int(-thresh / getEnv().GrainSize); }
    else { m_GrnFluxFirstNonZeroIdx = 0; }
  }

}//namespace
