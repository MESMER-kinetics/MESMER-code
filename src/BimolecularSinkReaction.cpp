//-------------------------------------------------------------------------------------------
//
// BimolecularSinkReaction.cpp
//
// Author: Dave Glowacki
// Date:   13 Dec 2011
//
// This file contains the implementation of the BimolecularSinkReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "BimolecularSinkReaction.h"

using namespace Constants ;
using namespace std;

namespace mesmer
{

  // Read the Molecular data from input stream.
  bool BimolecularSinkReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if(!ppReactantList)
      ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");      // Read reactant details.
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      cerr << "Cannot find 1st reactant molecule definition for bimolecular sink reaction " << getName() << ".";
      return false;
    }
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    Molecule* pMol2 = GetMolRef(ppReactant2);
    if(!pMol2)
    {
      cerr << "Cannot find 2nd reactant molecule definition for bimolecular sink reaction " << getName() << ".";
      return false;
    }

    // if excessReactantLocation=false, then pMol1 (the first rct
    // in the XML input) is the deficient reactant (m_rct1)

    Molecule* tmp_rct1 = pMol1;
    Molecule* tmp_rct2 = pMol2;

    if(!excessReactantLocation){
      m_rct1 = tmp_rct1;
      m_rct2 = tmp_rct2;
    } else {
      m_rct1 = tmp_rct2;
      m_rct2 = tmp_rct1;
    }

    if(!m_rct1){
      cerr << "the modelled molecule in the bimolecular sink reaction is undefined" << endl;
      return false;
    }
    if(!m_rct2){
      cerr << "the excess reactant in the bimolecular sink reaction is undefined" << endl;
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
      } else {
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

    return ReadRateCoeffParameters(ppReac) ;       // Read heat of reaction and rate parameters.
	}

  
	double BimolecularSinkReaction::rctsRovibronicGrnCanPrtnFn() { return m_rct1->getDOS().rovibronicGrnCanPrtnFn();}

  void BimolecularSinkReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) {

		cinfo << "adding reaction terms to the collision matrix " << endl;

    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    const int rShiftedGrains(m_rct1->getColl().reservoirShift());

    const int colloptrsize   = m_rct1->getColl().get_colloptrsize();
    const int forwardThreshE = get_EffGrnFwdThreshold();
    const int fluxStartIdx   = get_fluxFirstNonZeroIdx();

	m_MtxGrnKf.clear();
    m_MtxGrnKf.resize(colloptrsize , 0.0);

    for ( int i=fluxStartIdx, j = forwardThreshE; j < colloptrsize + rShiftedGrains; ++i, ++j) {
      int ii(rctLocation + j - rShiftedGrains) ;
	  double rtcf = m_GrainFlux[i] / rctDOS[j] ;
      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega * rtcf);  // Forward loss reaction.
	  m_MtxGrnKf[j - rShiftedGrains] = rtcf ;
    }
  }

  void BimolecularSinkReaction::AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) {
		cinfo << "adding reaction terms to the contracted basis matrix " << endl;
  }


  void BimolecularSinkReaction::calcGrainRateCoeffs(){

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

  
	void BimolecularSinkReaction::calcFluxFirstNonZeroIdx(void) {
		  double thresh = get_ThresholdEnergy();
			if(thresh<0.0){m_GrnFluxFirstNonZeroIdx = int(-thresh/getEnv().GrainSize);}
			else{m_GrnFluxFirstNonZeroIdx = 0;}
  }

  // Test k(T)
  void BimolecularSinkReaction::testRateConstant() {

		cinfo << "testing rate coefficient for bimolecular sink reaction " << endl;

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

  const int BimolecularSinkReaction::get_pdtsGrnZPE(){
    double zpe = m_pdt1->getDOS().get_zpe() - getEnv().EMin;
    if (m_pdt2) zpe += m_pdt2->getDOS().get_zpe();
    double grnZpe = zpe / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cinfo << "Grain zero point energy is negative in " << getName() << ".";

    return int(grnZpe);
  }

// calculate Keq, copied from IrreversibleExchangeReaction
	double BimolecularSinkReaction::calcEquilibriumConstant() {   // Calculate reaction equilibrium constant.
    // equilibrium constant:
    double Keq(0.0) ;
    const double beta = getEnv().beta ;

    // rovibronic partition function for products/reactants
    double Qrcts = rctsRovibronicGrnCanPrtnFn();
    double Qpdts = pdtsRovibronicGrnCanPrtnFn();

    // rovibronic partition function for reactants/products multiplied by translation contribution
    Qrcts *= translationalContribution(m_rct1->getStruc().getMass(), m_rct2->getStruc().getMass(), beta);
    Qpdts *= translationalContribution(m_pdt1->getStruc().getMass(), m_pdt2->getStruc().getMass(), beta);

    Keq = Qpdts / Qrcts;

    // Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells
    const double HeatOfReaction = getHeatOfReaction() ;
    const double _expon = -beta * HeatOfReaction;
    Keq *= exp(_expon);

    return Keq ;
	}

// calculate grain threshold, presently set simply to zero
  void BimolecularSinkReaction::calcEffGrnThresholds(void){           
		set_EffGrnFwdThreshold(0);
  }

}//namespace

