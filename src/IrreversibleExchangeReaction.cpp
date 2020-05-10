//-------------------------------------------------------------------------------------------
//
// IrreversibleExchangeReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IrreversibleExchangeReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IrreversibleExchangeReaction.h"
#include "gStructure.h"

using namespace Constants;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  // Read the Molecular data from input stream.
  bool IrreversibleExchangeReaction::InitializeReaction(PersistPtr ppReac)
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

    // Read rate coefficient parameters.
    if (!ReadRateCoeffParameters(ppReac)) return false;

    if (m_ERConc == 0.0 || IsNan(m_ERConc))
    {
      // If not already read in the MicroRateCalculator
      cinfo << "Not a unimolecular reaction: look for excess reactant concentration." << endl;
      if (!ReadExcessReactantConcentration(ppReac)) return false;
    }

    return true;
  }

  double IrreversibleExchangeReaction::calcEquilibriumConstant() {   // Calculate reaction equilibrium constant.
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

  // Add exchange reaction terms to collision matrix.
  void IrreversibleExchangeReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega)
  {
    // Call high pressure rate method because exchange reaction rate coefficients
    // do not depend on pressure. The extra logic is to control output. 
    bool writeStatus = getFlags().testRateConstantEnabled;
    getFlags().testRateConstantEnabled = false;

    HighPresRateCoeffs(NULL);
    const int jj = (*m_sourceMap)[get_pseudoIsomer()];
    (*CollOptr)[jj][jj] -= qd_real(rMeanOmega) * qd_real(m_MtxGrnKf[0]);

    getFlags().testRateConstantEnabled = writeStatus;
  }

  // Calculate high pressure rate coefficients at current T.
  void IrreversibleExchangeReaction::HighPresRateCoeffs(vector<double> *pCoeffs) {

    double k_forward(0.0);
    const double beta = getEnv().beta;

    string mrcID = string(m_pMicroRateCalculator->getID());

    if (mrcID == "RRKM" || mrcID == "SimpleRRKM") {
      // Check to see if there is a transition state defined.
      if (m_TransitionState) {

        // Energies of cells. 
        std::vector<double> cellEne;
        getCellEnergies(m_rct1->getEnv().MaxCell, m_rct1->getEnv().CellSize, cellEne);

        // Calculate ro-vibrational canonical partition function of reactants.
        vector<double> tmpDOS;
        m_rct1->getDOS().getCellDensityOfStates(tmpDOS);
        double prtfn = canonicalPartitionFunction(tmpDOS, cellEne, beta);
        tmpDOS.clear();
        m_rct2->getDOS().getCellDensityOfStates(tmpDOS);
        prtfn *= canonicalPartitionFunction(tmpDOS, cellEne, beta);

        // Calculate ro-vibrational canonical partition function of transitions state.
        tmpDOS.clear();
        m_TransitionState->getDOS().getCellDensityOfStates(tmpDOS);
        k_forward = canonicalPartitionFunction(tmpDOS, cellEne, beta);
        k_forward /= prtfn;

        // Calculate the translational partition function ratio and thence the rate coefficient.
        const double trans = translationalContribution(m_rct1->getStruc().getMass(), m_rct2->getStruc().getMass(), beta);
        k_forward /= trans;
        k_forward *= exp(-beta*get_ThresholdEnergy());
        k_forward *= SpeedOfLight_in_cm / beta;
      }
      else {
        throw(std::runtime_error(string("The irreversible exchange reaction" + getName() + "is defined without a transition state.\n")));
      }
    }
    else if (mrcID == "CanonicalRateCoefficient") {
      m_pMicroRateCalculator->calculateMicroCnlFlux(this);
      k_forward = m_CellFlux[0];
    }
    else {
      throw(std::runtime_error(string("The irreversible exchange reaction " + getName() + " is defined without either a transition state or Arrhenius parameters.\n")));
    }

    // Save high pressure rate coefficient for use in the construction of the collision operator.
    m_MtxGrnKf.clear();
    m_MtxGrnKf.push_back(k_forward*m_ERConc);

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
        << getName() << " = " << k_forward*m_ERConc << " s-1 (" << temperature << " K)" << endl;
      stest << "Canonical bimolecular rate constant of irreversible reaction "
        << getName() << " = " << k_forward << " cm^3/molec/s (" << temperature << " K)" << endl;
    }
  }

  const int IrreversibleExchangeReaction::get_rctsGrnZPE() {
    double grnZpe = (m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin) / getEnv().GrainSize; //convert to grains
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void IrreversibleExchangeReaction::calcEffGrnThresholds(void) {
    int TS_en = get_fluxGrnZPE();// see the comments in calcEffGrnThresholds under AssociationReaction.cpp  
    int rct_en = get_rctsGrnZPE();

    double thresh = get_ThresholdEnergy();
    if (thresh < 0.0) {
      set_EffGrnFwdThreshold(0 / getEnv().GrainSize);
    }
    else {
      set_EffGrnFwdThreshold((TS_en - rct_en) / getEnv().GrainSize);
    }
  }

  void IrreversibleExchangeReaction::calcFluxFirstNonZeroIdx(void) {
    double thresh = get_ThresholdEnergy();
    if (thresh < 0.0) { m_GrnFluxFirstNonZeroIdx = int(-thresh / getEnv().GrainSize); }
    else { m_GrnFluxFirstNonZeroIdx = 0; }
  }

}//namespace
