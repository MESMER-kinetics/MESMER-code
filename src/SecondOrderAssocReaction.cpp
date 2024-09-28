//-------------------------------------------------------------------------------------------
//
// SecondOrderAssocReaction.cpp
//
// Author: Struan Robertson
// Date:   24/Aug/2014
//
// This file contains the implementation of the SecondOrderAssocReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "SecondOrderAssocReaction.h"
#include "gWellProperties.h"
#include <math.h>

using namespace Constants;
using namespace std;

namespace mesmer
{
  //
  // Read molecular data for second order association reaction from input stream.
  // Note: the convention adopted here is that there are two (one repeated) reactants
  // and one product (adduct).
  // 
  bool SecondOrderAssocReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if (!ppReactantList)
      ppReactantList = ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1 = ppReactantList->XmlMoveTo("reactant");
    m_rct1 = GetMolRef(ppReactant1);
    m_rct2 = m_rct1;

    if (!m_rct1) {
      cerr << "The deficient reactant in the association reaction " << getName() << " is undefined." << endl;
      return false;
    }

    //Read product details.
    PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
    if (!ppProductList)
      ppProductList = ppReac; //Be forgiving; we can get by without a productList element

    PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
    m_pdt1 = GetMolRef(ppProduct1);
    if (!m_pdt1) {
      cerr << "Cannot find product molecule definition for association reaction " << getName() << ".";
      return false;
    }

    // Adjust the zero point energy to reflect second order.
    const double pseudoIsomerZPE = get_pseudoIsomer()->getDOS().get_zpe();
    double sourceTermZPE = 2.0 * pseudoIsomerZPE;
    get_pseudoIsomer()->getDOS().set_zpe(sourceTermZPE);

    // Read rate coefficient parameters.
    if (!ReadRateCoeffParameters(ppReac)) return false;

    if (!(get_concExcessReactant() > 0.0))
    {
      // If not already read in the MicroRateCalculator
      if (!ReadExcessReactantConcentration(ppReac)) return false;
    }

    return true;
  }

  // Reset zero point energy locations of the reactants such that
  // location of the pair is entirely on the pseudoisomer - in the
  // second order case this redundant and so we simply return the 
  // ZPE of the reactant.
  double SecondOrderAssocReaction::resetZPEofReactants() {
    return get_pseudoIsomer()->getDOS().get_zpe();
  }

  void SecondOrderAssocReaction::AddReactionTerms(qdMatrix* CollOptr,
    molMapType& isomermap,
    const double rMeanOmega)
  {
    // Get densities of states of the adduct for detailed balance.
    vector<double> pdtDOS;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS);

    // Locate isomers in system matrix.
    const int pdtLoc = isomermap[m_pdt1];
    const int jj = (*m_sourceMap)[get_pseudoIsomer()];

    // Get equilibrium constant.
    const qd_real Keq = qd_real(calcEquilibriumConstant());

    // Get Boltzmann distribution for detailed balance.
    vector<double> adductPopFrac; // Population fraction of the adduct
    const size_t pShiftedGrains(m_pdt1->getColl().reservoirShift());
    m_pdt1->getColl().normalizedGrnBoltzmannDistribution(adductPopFrac);

    qd_real DissRateCoeff(0.0), qdMeanOmega(rMeanOmega);

    const size_t pdtRxnOptPos(pdtLoc - pShiftedGrains);
    const size_t colloptrsize = m_pdt1->getColl().get_colloptrsize() + pShiftedGrains;
    const int reverseThreshE = get_EffGrnRvsThreshold();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    // Note: reverseThreshE will always be greater than pShiftedGrains here.

  // In following factors 2.0 and 4.0 appear. These arise from the the Taylor
  // expansion of the non-linear term about the the equilibrium point. 

    for (size_t i = reverseThreshE, j = fluxStartIdx; i < colloptrsize; ++i, ++j) {
      size_t ii(pdtRxnOptPos + i);
      size_t kk(i - pShiftedGrains);
      qd_real Flux(m_GrainFlux[j]), dos(pdtDOS[i]), addPop(adductPopFrac[kk]);
      (*CollOptr)[ii][ii] -= qdMeanOmega * Flux / dos;                                     // Loss of the adduct to the source
      (*CollOptr)[jj][ii] += qdMeanOmega * Flux * qd_real(2.0) * sqrt(Keq * addPop) / dos; // Reactive gain of the source
      (*CollOptr)[ii][jj]  = (*CollOptr)[jj][ii];                                          // Reactive gain (symmetrization)
      DissRateCoeff += Flux * addPop / dos;
    }
    (*CollOptr)[jj][jj] -= qd_real(4.0) * qdMeanOmega * DissRateCoeff * Keq;           // Loss of the source from detailed balance.
  }

  //
  // Add isomer reaction terms to contracted basis reaction matrix.
  //
  void SecondOrderAssocReaction::AddContractedBasisReactionTerms(qdMatrix* CollOptr, molMapType& isomermap)
  {
    // Get densities of states for detailed balance.
    vector<double> pdtDOS;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS);

    // Get equilibrium constant.
    const double Keq = calcEquilibriumConstant();

    // Get Boltzmann distribution for detailed balance.
    vector<double> adductPopFrac; // Population fraction of the adduct
    m_pdt1->getColl().normalizedGrnBoltzmannDistribution(adductPopFrac);

    const size_t pdtColloptrsize = m_pdt1->getColl().get_colloptrsize();
    const int reverseThreshE = get_EffGrnRvsThreshold();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    double DissRateCoeff(0.0);

    vector<double> RvsMicroRateCoef(pdtColloptrsize, 0.0);
    vector<double> CrsMicroRateCoef(pdtColloptrsize, 0.0);
    for (int i = fluxStartIdx, j = reverseThreshE, k = 0; j < pdtColloptrsize; ++i, ++j, ++k) {
      int mm = k + reverseThreshE;
      RvsMicroRateCoef[mm] = m_GrainFlux[i] / pdtDOS[mm];                          // Backward loss reaction. 
      CrsMicroRateCoef[mm] = RvsMicroRateCoef[mm] * sqrt(adductPopFrac[mm] * Keq); // Reactive gain from detailed balance.
      DissRateCoeff += RvsMicroRateCoef[mm] * adductPopFrac[mm];
    }

    // Calculate the elements of the product block.

    const size_t pdtLocation = isomermap[m_pdt1];
    const size_t pdtBasisSize = m_pdt1->getColl().get_nbasis();
    for (size_t i = 0, ii(pdtLocation), egvI(pdtColloptrsize - 1); i < pdtBasisSize; i++, ii++, --egvI) {
      (*CollOptr)[ii][ii] -= m_pdt1->getColl().matrixElement(egvI, egvI, RvsMicroRateCoef);
      for (size_t j = i + 1, jj(pdtLocation + j), egvJ(pdtColloptrsize - j - 1); j < pdtBasisSize; j++, jj++, --egvJ) {
        qd_real tmp = m_pdt1->getColl().matrixElement(egvI, egvJ, RvsMicroRateCoef);
        (*CollOptr)[ii][jj] -= tmp;
        (*CollOptr)[jj][ii] -= tmp;
      }
    }

    // Calculate the elements of the reactant block.

    const int jj = (*m_sourceMap)[get_pseudoIsomer()];
    (*CollOptr)[jj][jj] -= qd_real(DissRateCoeff * Keq);       // Loss of the source from detailed balance.

    // Calculate the elements of the cross blocks.

    vector<double> pdtBasisVector(pdtColloptrsize, 0.0);
    for (size_t i = 0, pdtEgv(pdtColloptrsize - 1); i < pdtBasisSize; i++, --pdtEgv) {
      size_t ii(pdtLocation + i);
      qd_real tmp(0.0);

      if (i == 0) {

        // Special case for equilibrium eigenvectors which obey a detailed balance relation.
        // SHR, 8/Mar/2009: are there other relations like this I wonder.

        qd_real elmti = (*CollOptr)[ii][ii];
        qd_real elmtj = (*CollOptr)[jj][jj];
        tmp = sqrt(elmti * elmtj);

      }
      else {

        // General case.

        m_pdt1->getColl().eigenVector(pdtEgv, pdtBasisVector);
        double sum = 0.0;
        for (size_t k(0), n(pdtColloptrsize - 1); k < (pdtColloptrsize - reverseThreshE); --n, k++) {
          sum += pdtBasisVector[n] * CrsMicroRateCoef[n];
        }

        tmp = qd_real(sum);
      }
      (*CollOptr)[ii][jj] += tmp;
      (*CollOptr)[jj][ii] += tmp;
    }

  }

  const int SecondOrderAssocReaction::get_rctsGrnZPE() {
    double grnZpe = (m_rct1->getDOS().get_zpe() - getEnv().EMin) / getEnv().GrainSize; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

}//namespace
