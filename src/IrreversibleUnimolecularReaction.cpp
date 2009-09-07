//-------------------------------------------------------------------------------------------
//
// IrreversibleUnimolecularReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IrreversibleUnimolecularReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IrreversibleUnimolecularReaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool IrreversibleUnimolecularReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if(!ppReactantList)
      ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      cerr << "Cannot find reactant molecule definition for irreversible reaction " << getName() << ".";
      return false;
    }

    // Save reactant as Molecule.

    Molecule* pColMol = pMol1;
    if(pColMol){
      m_rct1 = pColMol;
    } else {
      cerr << "Reactant of a irreversible reaction must be a colliding molecule";
      return false;
    }

    // Read product details. The detail of products may be absent or, may be needed
    // to calculate the microcanonical rates. If there are products, save them as
    // type Molecule.
    PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
    if(!ppProductList)
      ppProductList=ppReac; //Be forgiving; we can get by without a productList element

    PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);
      if (pMol1){
        m_pdt1 = pMol1;
      }
      else {
        cerr << "Irreversible reaction" << getName() << " has no product defined." << endl;
        return false;
      }

      Molecule* pMol2 = NULL ;
      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      if (ppProduct2) {
        pMol2 = GetMolRef(ppProduct2);
        if (pMol2){
          m_pdt2 = pMol2;
        }
        else {
          cinfo << "Irreversible reaction " << getName() << " has only one product defined." << endl;
        }
      }
    }

    // Read heat of reaction and rate parameters.
    return ReadRateCoeffParameters(ppReac) ;

  }

  //
  // Calculate reaction equilibrium constant.
  //
  double IrreversibleUnimolecularReaction::calcEquilibriumConstant() {

    // equilibrium constant:
    double Keq(0.0) ;
    const double beta = getEnv().beta ;

    // partition function for each products
    double Qpdts = pdtsRovibronicGrnCanPrtnFn();

    // rovibronic partition function for products multiplied by translation contribution
    if (m_pdt2){
      Qpdts *= translationalContribution(m_pdt1->getStruc().getMass(), m_pdt2->getStruc().getMass(), beta);
    }

    // rovibronic partition function for reactant
    const double Qrct1 = m_rct1->getDOS().rovibronicGrnCanPrtnFn() ;

    Keq = Qpdts / Qrct1;

    // Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells
    const double HeatOfReaction = getHeatOfReaction() ;
    const double _expon = -beta * HeatOfReaction;
    Keq *= exp(_expon) ;

    return Keq ;
  }

  //
  // Add dissociation reaction terms to collision matrix.
  //
  void IrreversibleUnimolecularReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) {
    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    int rNGG(m_rct1->getColl().getNumberOfGroupedGrains());
    const bool rIsCemetery(m_rct1->getColl().isCemetery());
    const int rShiftedGrains(rNGG == 0 ? 0 : (rIsCemetery ? rNGG : rNGG - 1));

    const int colloptrsize   = m_rct1->getColl().get_colloptrsize();
    const int forwardThreshE = get_EffGrnFwdThreshold();
    const int fluxStartIdx   = get_fluxFirstNonZeroIdx();

    for ( int i=fluxStartIdx, j = forwardThreshE, k=0; j < colloptrsize; ++i, ++j, ++k) {
      int ll = k + forwardThreshE;
      int ii(rctLocation + ll - rShiftedGrains) ;
      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * m_GrainFlux[i] / rctDOS[ll]);                     // Forward loss reaction.
    }
  }

  //
  // Add contracted basis set reaction terms to the reaction matrix.
  //
  void IrreversibleUnimolecularReaction::AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) {

    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;

    const int rctColloptrsize = m_rct1->getColl().get_colloptrsize();
    const int forwardThreshE  = get_EffGrnFwdThreshold();
    const int fluxStartIdx    = get_fluxFirstNonZeroIdx();

    vector<double> fwdMicroRateCoef(rctColloptrsize, 0.0) ;
    for ( int i=fluxStartIdx, j = forwardThreshE, k=0; j < rctColloptrsize; ++i, ++j, ++k) {
      int ll = k + forwardThreshE;
      fwdMicroRateCoef[ll] = m_GrainFlux[i] / rctDOS[ll] ;    // Forward loss reaction.
    }

    // Locate isomers in system matrix.
    const int rctLocation = isomermap[m_rct1] ;

    // Calculate the elements of the reactant block.
    const int rctBasisSize = static_cast<int>(m_rct1->getColl().get_nbasis()) ;

    for (int i=0, ii(rctLocation), egvI(rctColloptrsize-1) ; i < rctBasisSize ; i++, ii++, --egvI) {
      (*CollOptr)[ii][ii] -= m_rct1->getColl().matrixElement(egvI, egvI, fwdMicroRateCoef) ;
      for (int j=i+1, jj(rctLocation + j), egvJ(rctColloptrsize-j-1) ; j < rctBasisSize ; j++, jj++, --egvJ) {
        qd_real tmp = m_rct1->getColl().matrixElement(egvI, egvJ, fwdMicroRateCoef) ;
        (*CollOptr)[ii][jj] -= tmp ;
        (*CollOptr)[jj][ii] -= tmp ;
      }
    }
  }

  //
  // Calculate grained forward k(E)s from transition state flux
  //
  void IrreversibleUnimolecularReaction::calcGrainRateCoeffs(){
    vector<double> rctGrainDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS) ;

    calcEffGrnThresholds();
    const int forwardTE = get_EffGrnFwdThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    const int MaximumGrain = (getEnv().MaxGrn-fluxStartIdx);
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain , 0.0);

    // calculate forward k(E)s from flux
    for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i];
    }

    // the code that follows is for printing the forward k(E)s
    if (getFlags().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().testRateConstantEnabled)
      testRateConstant();
  }

  void IrreversibleUnimolecularReaction::calcFluxFirstNonZeroIdx(void) {
    //Use base class function(which references product properties) only with reverse ILT
    if(isReverseReactionILT_Ea())
      Reaction::calcFluxFirstNonZeroIdx();
    else
      m_GrnFluxFirstNonZeroIdx = 0;
  }

  // Test k(T)
  void IrreversibleUnimolecularReaction::testRateConstant() {

    double k_forward(0.0);
    vector<double> rctGrainDOS, rctGrainEne;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS);
    m_rct1->getDOS().getGrainEnergies(rctGrainEne);
    const int MaximumGrain = (getEnv().MaxGrn-get_fluxFirstNonZeroIdx());
    const double beta = getEnv().beta;
    const double temperature = 1. / (boltzmann_RCpK * beta);

    for(int i(0); i < MaximumGrain; ++i)
      k_forward += m_GrainKfmc[i] * exp( log(rctGrainDOS[i]) - beta * rctGrainEne[i]);

    const double rctprtfn = canonicalPartitionFunction(rctGrainDOS, rctGrainEne, beta);
    k_forward /= rctprtfn;
    set_fwdGrnCanonicalRate(k_forward);

    ctest << endl << "Canonical pseudo first order forward rate constant of irreversible reaction "
      << getName() << " = " << get_fwdGrnCanonicalRate() << " s-1 (" << temperature << " K)" << endl;
  }

  void IrreversibleUnimolecularReaction::calcEffGrnThresholds(void){       
    int TS_en = get_fluxGrnZPE();
    int rct_en = m_rct1->getColl().get_grnZPE();

    //This function uses product properties, although the reaction is irreversible
    //This is only appropriate for reverse ILT.
    if(!isReverseReactionILT_Ea())
    {
      set_EffGrnFwdThreshold(TS_en-rct_en);
      return;
    }

    double RxnHeat   = getHeatOfReaction();    // see the comments in
    double threshold = get_ThresholdEnergy();  // calcEffGrnThresholds under AssociationReaction.cpp

    int pdtsGrnZPE = get_pdtsGrnZPE();
    int rctGrnZPE  = m_rct1->getColl().get_grnZPE();
    int GrainedRxnHeat  = pdtsGrnZPE - rctGrnZPE;

    if(threshold<0.0){
      set_EffGrnFwdThreshold(0);
    }
    else if(threshold>0.0 && threshold<RxnHeat){// if the reverse threshold energy is negative
      set_EffGrnFwdThreshold( GrainedRxnHeat);  // forward grained flux threshold energy = heat of reaction
    }
    else{
      set_EffGrnFwdThreshold(TS_en-rct_en);
    }
  }

  //
  // Get products cell density of states.
  //
  void IrreversibleUnimolecularReaction::getPdtsCellDensityOfStates(vector<double> &cellDOS) {
    countDimerCellDOS(m_pdt1->getDOS(), m_pdt2->getDOS(), cellDOS);
  }

  const int IrreversibleUnimolecularReaction::get_pdtsGrnZPE(){
    double zpe = m_pdt1->getDOS().get_zpe() - getEnv().EMin;
    if (m_pdt2) zpe += m_pdt2->getDOS().get_zpe();
    double grnZpe = zpe / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cinfo << "Grain zero point energy is negative in " << getName() << ".";

    return int(grnZpe);
  }

  //
  // Calculate the rovibrational density of states of products.
  //
  bool IrreversibleUnimolecularReaction::calcPdtsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)
  {
    std::vector<double> pdtsCellDOS;
    getPdtsCellDensityOfStates(pdtsCellDOS);

    std::vector<double> shiftedCellDOS;
    std::vector<double> shiftedCellEne;
    const int MaximumCell = getEnv().MaxCell;
    const int cellOffset = m_pdt1->getDOS().get_cellOffset(); // ** temporary statement to get cellOffset from one of the molecules.
    std::vector<double> pdtsCellEne;
    getCellEnergies(MaximumCell, pdtsCellEne);
    shiftCells(MaximumCell, cellOffset, pdtsCellDOS, pdtsCellEne, shiftedCellDOS, shiftedCellEne);

    const string catName = m_pdt1->getName() + " + " + m_pdt2->getName();

    if (getFlags().cyclePrintCellDOS){
      ctest << endl << "Cell rovibronic density of states of " << catName << endl << "{" << endl;
      for (int i = 0; i < MaximumCell; ++i){
        formatFloat(ctest, pdtsCellEne[i],  6,  15) ;
        formatFloat(ctest, pdtsCellDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
      getFlags().cyclePrintCellDOS = false;
    }

    calcGrainAverages(getEnv().MaxGrn, getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, grainDOS, grainEne);

    if (getFlags().cyclePrintGrainDOS){
      ctest << endl << "Grain rovibronic density of states of " << catName << endl << "{" << endl;
      for (int i = 0; i < getEnv().MaxGrn; ++i){
        formatFloat(ctest, grainEne[i],  6,  15) ;
        formatFloat(ctest, grainDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
      getFlags().cyclePrintGrainDOS = false;
    }

    return true;
  }

  //
  // Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
  //
  double IrreversibleUnimolecularReaction::rctsRovibronicGrnCanPrtnFn() { return m_rct1->getDOS().rovibronicGrnCanPrtnFn();}
  double IrreversibleUnimolecularReaction::pdtsRovibronicGrnCanPrtnFn() {
    if (!m_pdt2){ // Irreversible isomerization
      return m_pdt1->getDOS().rovibronicGrnCanPrtnFn();
    }
    else{         // Irreversible dissociation
      vector<double> pdtGrainDOS;
      vector<double> pdtGrainEne;
      calcPdtsGrainDensityOfStates(pdtGrainDOS, pdtGrainEne);

      // Calculate the rovibronic partition function based on the grain DOS
      // The following catches the case where the molecule is a single atom
      double CanPrtnFn = max(canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, getEnv().beta), 1.0) ;
      if (CanPrtnFn == 1.0){
        // Electronic partition function for atom is accounted here.
        CanPrtnFn = double(m_pdt1->getDOS().getSpinMultiplicity() * m_pdt2->getDOS().getSpinMultiplicity()) ;
      }
      return CanPrtnFn ;
    }
  }

}//namespace
