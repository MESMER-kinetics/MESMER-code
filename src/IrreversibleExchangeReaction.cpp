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

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  // Read the Molecular data from input stream.
  bool IrreversibleExchangeReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");      // Read reactant details.
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      cerr << "Cannot find 1st reactant molecule definition for association reaction " << getName() << ".";
      return false;
    }
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    Molecule* pMol2 = GetMolRef(ppReactant2);
    if(!pMol2)
    {
      cerr << "Cannot find 2nd reactant molecule definition for association reaction " << getName() << ".";
      return false;
    }

    // if deficientReactantLocation=true, then pMol1 (the first rct
    // in the XML input) is the deficient reactant (m_rct1)

    ModelledMolecule* tmp_rct1 = dynamic_cast<ModelledMolecule*>(pMol1);
    ModelledMolecule* tmp_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);

    if(deficientReactantLocation){
      m_rct1 = tmp_rct1;
      m_rct2 = tmp_rct2;
    }
    else {
      m_rct1 = tmp_rct2;
      m_rct2 = tmp_rct1;
    }

    if(!m_rct1){
      cerr << "the deficient reactant in the association reaction is undefined" << endl;
      return false;
    }
    if(!m_rct2){
      cerr << "the excess reactant in the association reaction is undefined" << endl;
      return false;
    }

    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");     // Read product details. Save them as type Molecule
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);
      if (pMol1) {
        m_pdt1 = dynamic_cast<ModelledMolecule*>(pMol1) ;
      } else {
        cerr << "Exchange reaction" << getName() << " has no products defined.";
      }

      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      if (ppProduct2) {
        pMol2 = GetMolRef(ppProduct2);
        if (pMol2) {
          m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2) ;
        } else {
          cerr << "Exchange reaction " << getName() << " has only one product defined.";
        }
      }
    }

    return ReadRateCoeffParameters(ppReac) ;       // Read heat of reaction and rate parameters.

  }

  double IrreversibleExchangeReaction::calcEquilibriumConstant() {   // Calculate reaction equilibrium constant.
    double Keq(0.0) ;
    return Keq ;
  }


  void IrreversibleExchangeReaction::AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, 
    const double rMeanOmega)              // Add exchange reaction terms to collision matrix.
  {
    const int jj     = (*m_sourceMap)[get_pseudoIsomer()] ;
    (*CollOptr)[jj][jj] -= qd_real(rMeanOmega * get_fwdGrnCanonicalRate());
  }

  bool IrreversibleExchangeReaction::calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)    // Calculate rovibrational reactant DOS
  {
    std::vector<double> rctsCellDOS;
    getRctsCellDensityOfStates(rctsCellDOS);

    std::vector<double> shiftedCellDOS;
    std::vector<double> shiftedCellEne;
    const int MaximumCell = getEnv().MaxCell;
    const int cellOffset = get_pseudoIsomer()->get_cellOffset();
    std::vector<double> rctsCellEne;
    getCellEnergies(MaximumCell, rctsCellEne);
    shiftCells(MaximumCell, cellOffset, rctsCellDOS, rctsCellEne, shiftedCellDOS, shiftedCellEne);

    string catName = m_rct1->getName() + " + " + m_rct2->getName();
    calcGrainAverages(getEnv().MaxGrn, getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, grainDOS, grainEne, catName);

    return true;
  }

  // Get reactants cell density of states.
  void IrreversibleExchangeReaction::getRctsCellDensityOfStates(vector<double> &cellDOS) {
    get_rctsDensityOfStatesCalculator()->countDimerCellDOS(m_rct1, m_rct2, cellDOS);
  }

  // Calculate grained forward and reverse k(E)s from transition state flux

  void IrreversibleExchangeReaction::calcGrainRateCoeffs(){

    vector<double> rctGrainDOS;
    vector<double> rctGrainEne;
    calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

    calcEffGrnThresholds();
    const int forwardTE = get_EffGrnFwdThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    const int MaximumGrain = (getEnv().MaxGrn - fluxStartIdx);
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain , 0.0);

    for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i];
    }

    if (getFlags().kfEGrainsEnabled){               // printing of the forward k(E)s
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }

    testRateConstant();  // always execute this routine for irreversible exchange reactions in order to derive
                         // the irreversible collision matrix element
  }

  void IrreversibleExchangeReaction::testRateConstant() {    // Test k(T)

    const int MaximumGrain = (getEnv().MaxGrn-get_fluxFirstNonZeroIdx());
    const double beta = getEnv().beta;
    vector<double> rctsGrainDOS;
    vector<double> rctsGrainEne;
    double k_forward(0.0);

    calcRctsGrainDensityOfStates(rctsGrainDOS, rctsGrainEne);
    for(int i(0); i < MaximumGrain; ++i){
      k_forward += m_GrainKfmc[i] * exp( log(rctsGrainDOS[i]) - beta * rctsGrainEne[i]); 
    }
    const double prtfn = canonicalPartitionFunction(rctsGrainDOS, rctsGrainEne, beta);
    const double trans = translationalContribution(m_rct1->getMass(), m_rct2->getMass(), beta);
    k_forward /= prtfn;
    k_forward /= trans;
    k_forward *= m_ERConc;
    set_fwdGrnCanonicalRate(k_forward);

    if (getFlags().testRateConstantEnabled){
      const double temperature = 1. / (boltzmann_RCpK * beta);
      ctest << endl << "Canonical pseudo first order rate constant of irreversible reaction " 
        << getName() << " = " << get_fwdGrnCanonicalRate() << " s-1 (" << temperature << " K)" << endl;
      ctest << "Canonical bimolecular rate constant of irreversible reaction " 
        << getName() << " = " << get_fwdGrnCanonicalRate()/m_ERConc << " cm^3/mol/s (" << temperature << " K)" << endl;
    }
  }

  const int IrreversibleExchangeReaction::get_rctsGrnZPE(){
    double grnZpe = (m_rct1->get_zpe()+m_rct2->get_zpe()-getEnv().EMin) / getEnv().GrainSize ; //convert to grains
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void IrreversibleExchangeReaction::calcEffGrnThresholds(void){           // see the comments in
    double thresh = get_ThresholdEnergy();  // calcEffGrnThresholds under AssociationReaction.cpp
    int TS_en = get_fluxGrnZPE();
    int rct_en = get_rctsGrnZPE();
    if(thresh<0.0){
      set_EffGrnFwdThreshold(0);
    }
    else{
      set_EffGrnFwdThreshold(TS_en-rct_en);
    }
  }

}//namespace
