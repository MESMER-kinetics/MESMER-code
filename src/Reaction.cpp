//-------------------------------------------------------------------------------------------
//
// Reaction.cpp
//
// Author: Struan Robertson
// Date:   23/Feb/2003
//
// This file contains the implementation of the Reaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "Reaction.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  Reaction::Reaction(MoleculeManager *pMoleculeManager):
    m_Name(),
    m_pMoleculeManager(pMoleculeManager),
    m_srct(NULL),
    m_rct1(NULL),
    m_rct2(NULL),
    m_pdt1(NULL),
    m_pdt2(NULL),
    m_TransitionState(NULL),
    m_reactiontype(ERROR_REACTION),
    m_HeatOfReaction(0.0),
    m_kfwd(0.0),
    m_CellKfmc(),
    m_GrainKfmc(),
    m_ActEne(0.0),
    m_PreExp(0.0),
    m_ppPersist(),
    m_pMicroRateCalculator(NULL)
  {}

  Reaction::~Reaction(){}

  /*
  Reaction::Reaction(const Reaction& reaction) {
  // Copy constructor - define later SHR 23/Feb/2003
  }

  Reaction& Reaction::operator=(const Reaction& reaction) {
  // Assignment operator - define later SHR 23/Feb/2003

  return *this ;
  }
  */

  //
  // Read the Molecular data from inout stream.
  //
  bool Reaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;
    stringstream errorMsg;

    //Read reaction ID
    const char* id = ppReac->XmlReadValue("id");
    if(id) m_Name = id; //Continues if reaction id not found
    else{
      errorMsg << "Reaction ID not found\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    Molecule* pMol1(NULL) ;
    Molecule* pMol2(NULL) ;

    //Read reactants
    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
    pMol1 = GetMolRef(ppReactant1);
    if(!pMol1) return false;
    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    if(ppReactant2)
    {
      pMol2 = GetMolRef(ppReactant2);
      if(!pMol2) return false;
    }

    //Put the CollidingMolecule into m_rct1, even if it is second in datafile
    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol)
    {
      m_rct1 = pColMol;
      m_rct2 = dynamic_cast<ModelledMolecule*>(pMol2);
    }
    else{
      pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
      if(!pColMol){
        errorMsg << "Either " << pMol1->getName() << " or "
            << pMol2->getName() <<" has to be a colliding molecule";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
      m_rct1 = pColMol;
      m_rct2 = dynamic_cast<ModelledMolecule*>(pMol1);
    }

    if (pMol1 && pMol2){
      // Create a SuperMolecule from reactants
      SuperMolecule* pSMol(NULL);
      pSMol = m_pMoleculeManager->addSuperMol(); // add a blank SuperMolecule pointer to molmap

      m_srct = pSMol;
      m_srct->setMembers(m_rct1, m_rct2);
    }

    //Read products ... if any.
    pMol1=NULL; pMol2=NULL;
    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);

      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      pMol2 = GetMolRef(ppProduct2);

      //Put the Colliding Molecule into m_pdt1, even if it is second in datafile
      pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
      if(pColMol)
      {
        m_pdt1 = pColMol;
        m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol2);
      }
      else
      {
        pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
        if(!pColMol)
        {
          stringstream errorMsg;
          errorMsg << "Either " << pMol1->getName() << " or "
            << pMol2->getName() <<" has to be a modelled molecule";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
          return false;
        }
        m_pdt1 = pColMol;
        m_pdt2 = dynamic_cast<ModelledMolecule*>(pMol1);
      }
    }

    // Read the transition state (if present)
    PersistPtr ppTransitionState = ppReac->XmlMoveTo("me:transitionState") ;
    if (ppTransitionState)
    {
      PersistPtr ppmol = ppTransitionState->XmlMoveTo("molecule");
      if(ppmol)
      {
        const char* pRef = ppmol->XmlReadValue("ref");
        if(!pRef)
            return false;
        m_TransitionState = dynamic_cast<TransitionState*>(m_pMoleculeManager->find(pRef));
      }

      /* It would be better to use the ZPEs rather than threshold
      if(m_TransitionState)
      {
      const char* pthreshtxt = ppReac->XmlReadValue("me:threshold",false);
      if(pthreshtxt)
      {
      stringstream ss(pthreshtxt);
      ss >> m_E0;
      }
      */
    }

    // Read the heat of reaction (if present).
    const char* pHeatRxntxt = ppReac->XmlReadValue("me:HeatOfReaction",false);
    if (pHeatRxntxt){
      stringstream s1(pHeatRxntxt);
      s1 >> m_HeatOfReaction ;
    }
    else{ // calculate HeatOfReaction
      double zpe_pdt1 = m_pdt1 ? m_pdt1->get_zpe() : 0.;
      double zpe_pdt2 = m_pdt2 ? m_pdt2->get_zpe() : 0.;
      double zpe_rct1 = m_rct1 ? m_rct1->get_zpe() : 0.;
      double zpe_rct2 = m_rct2 ? m_rct2->get_zpe() : 0.;
      m_HeatOfReaction = zpe_pdt1 + zpe_pdt2 - zpe_rct1 - zpe_rct2;

    }

    //Read in Kinetic rate parameters, if present
    m_ActEne = std::numeric_limits<double>::quiet_NaN();//means not set

    const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy",false);
    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential",false);
    if (pActEnetxt && pPreExptxt)
    {
      stringstream s1(pActEnetxt); s1 >> m_ActEne ;
      stringstream s2(pPreExptxt); s2 >> m_PreExp ;
    }

    // Classify the reaction.

    if (m_rct1 && m_pdt1 && !m_rct2 && !m_pdt2)
      m_reactiontype = ISOMERIZATION ;
    else if (m_rct1 && m_pdt1 && m_rct2 && !m_pdt2)
      m_reactiontype = ASSOCIATION ;
    else if (m_rct1)
      m_reactiontype = DISSOCIATION ;
    else {
      m_reactiontype = ERROR_REACTION ;
      stringstream errorMsg; errorMsg << "Unknown combination of reactants and products";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    // Determine the method of MC rate coefficient calculation.
    const char* pMCRCMethodtxt = ppReac->XmlReadValue("me:MCRCMethod") ;
    if(pMCRCMethodtxt)
    {
      m_pMicroRateCalculator = MicroRateCalculator::Find(pMCRCMethodtxt);
      if(!m_pMicroRateCalculator)
      {
        stringstream errorMsg;
        errorMsg << "Unknown method " << pMCRCMethodtxt
          << " for the determination of Microcanonical rate coefficients in reaction "
          << m_Name;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
    }


    return true;
  }

  //
  // Locate molecule in molecular map.
  //
  Molecule* Reaction::GetMolRef(PersistPtr pp)
  {
    Molecule* pMol;
    if(!pp)
        return NULL;
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol){
      return NULL;
    }

    const char* pRef = ppmol->XmlReadValue("ref");
    if(pRef)
      pMol = m_pMoleculeManager->find(pRef);

    if(!pMol)
    {
      stringstream errorMsg;
      errorMsg << "Unknown molecule: " << pRef;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return NULL;
    }

    return pMol;
  }

  //
  // Returns the unimolecular species in each reaction, i.e. for association
  // (source term) or dissociation (sink term) reaction one species is returned,
  // for an isomerization reaction two species are returned.
  //
  int Reaction::get_unimolecularspecies(vector<CollidingMolecule *> &unimolecularspecies) const
  {
    switch(m_reactiontype) {
    case ISOMERIZATION :
      unimolecularspecies.push_back(m_rct1) ;
      unimolecularspecies.push_back(m_pdt1) ;
      return 3;

    case ASSOCIATION :
      unimolecularspecies.push_back(m_pdt1) ;
      return 2;

    case DISSOCIATION :
      unimolecularspecies.push_back(m_rct1) ;
      return 1;

    default :
      return 0;
    }
    return 0;
  }

  // Returns the bi-molecular speices (reactants) for association reaction.
  int Reaction::get_bi_molecularspecies(SuperMolecule* bi_mol) const
  {
    switch(m_reactiontype) {
    case ISOMERIZATION :
      return 3;

    case ASSOCIATION :
      bi_mol = m_srct ;
      return 2;

    case DISSOCIATION :
      return 1;

    default :
      return 0;
    }
    return 0;
  }

  //
  // Get the prinincipal source reactant (i.e. reactant not in excess) if it exists.
  // (Not sure if this is a good idea, may be better to pass a Map in.)
  //
  CollidingMolecule *Reaction::get_pseudoIsomer() const
  {
    CollidingMolecule *pseudoIsomer = NULL ;
    if(m_reactiontype == ASSOCIATION) pseudoIsomer = m_pdt1 ;
    return pseudoIsomer ;
  }

  //
  // Calculate reaction equilibrium constant.
  //
  double Reaction::calcEquilibriumConstant(const MesmerEnv &mEnv) {

    double Keq(0.0) ;

    // Get Canonical partition functions.

    double Qrct1 = m_rct1->grnCanPrtnFn(mEnv) ;
    double Qpdt1 = m_pdt1->grnCanPrtnFn(mEnv) ;

    double Qrct2 = (m_rct2)? m_rct2->grnCanPrtnFn(mEnv) : 0.0 ;
    double Qpdt2 = (m_pdt2)? m_pdt2->grnCanPrtnFn(mEnv) : 0.0 ;

    // Calculate the equilibrium constant.

    double beta = mEnv.beta ;

    if      ( m_rct2 &&  m_pdt2 ) Keq = (Qpdt1*Qpdt2)/(Qrct1*Qrct2) ;     // Exchange reaction.
    else if ( m_rct2 && !m_pdt2 ) Keq =  Qpdt1/(Qrct1*Qrct2) ;            // Assocations reaction.
    else if (!m_rct2 &&  m_pdt2 ) Keq = (Qpdt1*Qpdt2)/Qrct1 ;             // Dissociation reaction.
    else if (!m_rct2 && !m_pdt2 ) Keq =  Qpdt1/Qrct1 ;                    // Isomerization reaction.
    else { //no chance to get here
      stringstream errorMsg;
      errorMsg << "Error: calculating equilibrium constant for reaction," << endl
               << "unknown reactant/product combination" << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    Keq *= exp(-beta*m_HeatOfReaction) ;

    return Keq ;
  }

  //
  // Access microcanoincal rate coefficients.
  //
  void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc, const MesmerEnv &mEnv) {
    calcGrnAvrgMicroRateCoeffs(mEnv);
    kmc = m_GrainKfmc ;
  }

  //
  // Calculate grain averaged microcanonical rate coefficients.
  //
  bool Reaction::calcGrnAvrgMicroRateCoeffs(const MesmerEnv &mEnv) {

    // Calculate microcanonical rate coefficients.
    if (m_CellKfmc.size()==0)
    {
      if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this, m_CellKfmc, mEnv) ||
        (mEnv.microRateEnabled && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_CellKfmc, m_ppPersist, mEnv)))
        return false;
    }

    // Calculate Grain averages of microcanonical rate coefficients.
    if (m_GrainKfmc.size()==0)
        return grnAvrgMicroRateCoeffs(mEnv) ;
    return true;
  }

  //
  // Access microcanonical rate coefficients - cell values are averaged
  // to give grain values. This code is similar to that in Molecule.cpp
  // and this averaging should be done there. SHR 19/Sep/2004.
  //
  bool Reaction::grnAvrgMicroRateCoeffs(const MesmerEnv &mEnv) {

    int MaximumGrain = mEnv.MaxGrn;
    double currentGrainSize = mEnv.GrainSize;
    m_GrainKfmc.resize(MaximumGrain);

    // Extract density of states of equilibrium molecule.

    vector<double> cellDOS(mEnv.MaxCell,0.0) ;
    m_rct1->getCellDensityOfStates(cellDOS, mEnv) ;

    // Check that there are enough cells.

    if (currentGrainSize < 1) {
      stringstream errorMsg;
      errorMsg << "Not enought Cells to produce requested number of Grains.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i = 0; i < MaximumGrain ; ++i ) {

      int idx3 = idx1 ;

      // Calculate the number of states in a grain.

      double smt = .0 ;
      for (int j = 0 ; j < currentGrainSize ; ++j, ++idx1 )
        smt += cellDOS[idx1] ;

      // Calculate average energy of the grain if it contains sum states.

      if ( smt > 0.0 ) {

        double smat = .0;
        for (int j= 0 ; j < currentGrainSize ; ++j, ++idx3 )
            smat += m_CellKfmc[idx3] * cellDOS[idx3] ;

        m_GrainKfmc[idx2] = smat/smt ;
        idx2++ ;
      }
    }

    // Issue warning if number of grains produced is less that requested.
    if ( idx2 < MaximumGrain ) {
      stringstream errorMsg;
      errorMsg << "Number of grains produced is less than requested" << endl
        << "Number of grains requested: " << MaximumGrain << endl
        << "Number of grains produced : " << idx2 << ".";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    return true;
  }

  //
  // Add microcanonical terms to collision operator
  //
  void Reaction::AddMicroRates(dMatrix *CollOptr,
      isomerMap &isomermap,
      sourceMap &sourcemap,
      const double rMeanOmega,
      const MesmerEnv &mEnv)
  {

    // Calculate Microcanonical rate coefficients.

    calcGrnAvrgMicroRateCoeffs(mEnv) ;

    // Add microcanonical rates to the collision operator.

    switch(m_reactiontype) {
    case ISOMERIZATION :
      AddIsomerReactionTerms(CollOptr, isomermap, rMeanOmega, mEnv) ;
      break;

    case ASSOCIATION :
      AddAssocReactionTerms(CollOptr, isomermap, sourcemap, rMeanOmega, mEnv) ;
      break;

    case DISSOCIATION :
      AddDissocReactionTerms(CollOptr, isomermap, rMeanOmega) ;
      break;

    default :
      stringstream errorMsg;
      errorMsg << "Unknown reaction type";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);

    }
  }

  //
  // Add isomer reaction terms to collision matrix.
  //
  void Reaction::AddIsomerReactionTerms(dMatrix         *CollOptr,
                                        isomerMap       &isomermap,
                                        const double    rMeanOmega,
                                        const MesmerEnv &mEnv)
  {
    // Locate isomers in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    const int pdtLocation = isomermap[m_pdt1] ;

    // Get densities of states for detailed balance.
    const int MaximumGrain = mEnv.MaxGrn;
    vector<double> rctDOS(MaximumGrain, 0.0) ;
    vector<double> pdtDOS(MaximumGrain, 0.0) ;

    m_rct1->grnDensityOfStates(rctDOS, mEnv) ;
    m_pdt1->grnDensityOfStates(pdtDOS, mEnv) ;

    const int idx = m_pdt1->get_grnZpe() - m_rct1->get_grnZpe() ;
    for ( int i = max(0,-idx) ; i < min(MaximumGrain,(MaximumGrain-idx)) ; ++i ) {
      int ll = i + idx ;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation + i) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[ll] ;                            // Forward loss reaction.
      (*CollOptr)[jj][jj] -= rMeanOmega * m_GrainKfmc[ll]*rctDOS[ll]/pdtDOS[i] ;       // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(rctDOS[ll]/pdtDOS[i]) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                     // Reactive gain.
    }
  }

  //
  // Add (REVERSIBLE) association reaction terms to collision matrix.
  //
  void Reaction::AddAssocReactionTerms(dMatrix         *CollOptr,
                                       isomerMap       &isomermap,
                                       sourceMap       &sourcemap,
                                       const double    rMeanOmega,
                                       const MesmerEnv &mEnv)
  {
    // Locate isomers in system matrix.
    const int pdtLoc  = isomermap[m_pdt1] ;
    const int srcLoc  = sourcemap[m_srct] ;

    // Get equilibrium constant.
    double Keq = calcEquilibriumConstant(mEnv) ;

    // Get Boltzmann distribution for detailed balance.
    const int MaximumGrain = mEnv.MaxGrn ;
    vector<double> srcBoltz(MaximumGrain, 0.0) ;
    m_srct->grnBoltzDist(srcBoltz, mEnv) ;

    int sL(srcLoc) ;
    double DissRateCoeff(0.0) ;

    const int idx = int(m_HeatOfReaction / mEnv.GrainSize); //any possible leak? need to test.
    //the same number of atoms on reactants and product
    for ( int i = max(0,-idx) ; i < min(MaximumGrain,(MaximumGrain-idx)) ; ++i ) {
      int sG = m_srct->get_grnZpe() + max(0,idx);

      /*for well-dependent grain vector one should count from grnZpe (related with the minimum of all wells) of
      that specific isomer*/
      int ll = i + idx ; int pL(pdtLoc + ll) ;

      (*CollOptr)[pL][pL] -= rMeanOmega * m_GrainKfmc[sG] ;                        // Forward loss reaction.
      (*CollOptr)[pL][sL]  = rMeanOmega * m_GrainKfmc[sG]*sqrt(srcBoltz[sG]/Keq) ; // Reactive gain.
      (*CollOptr)[sL][pL]  = (*CollOptr)[pL][sL] ;                                 // Reactive gain.
      DissRateCoeff       += m_GrainKfmc[sG]*srcBoltz[sG] ;
    }
    (*CollOptr)[sL][sL] -= DissRateCoeff/Keq ;       // Backward loss reaction from detailed balance.

  }

  //
  // Add dissociation reaction terms to collision matrix.
  //
  void Reaction::AddDissocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[m_rct1] ;

    for ( int i = 0 ; i < m_rct1->get_colloptrsize() ; ++i ) {
        int ii(rctLocation + i) ;
        (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[i] ;                            // Forward loss reaction.
    }
  }

}//namespace
