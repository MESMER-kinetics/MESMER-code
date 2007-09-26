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
#include "System.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  Reaction::Reaction(MoleculeManager *pMoleculeManager):
            m_pMoleculeManager(pMoleculeManager),
            m_Reactant(NULL),
            m_Reactant2(NULL),
            m_Product(NULL),
            m_Product2(NULL),
            m_TransitionState(NULL),
            m_kfwd(0.0),
            m_CellKfmc(),
            m_GrainKfmc()
            {}

  Reaction::~Reaction()
  {
    if(m_pMicroRateCalculator) delete m_pMicroRateCalculator;
  }

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
  bool Reaction::InitializeReaction(System* pSys, PersistPtr ppReac)
  {
    m_pSys = pSys;
    m_ppPersist = ppReac;

    //Read reaction ID
    const char* id = ppReac->XmlReadValue("id");
    if(id)
        m_Name = id; //Continues if reaction id not found

    Molecule* pMol1 = NULL;
    Molecule* pMol2 = NULL;

    //Read reactants
    PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
    pMol1 = GetMolRef(ppReactant1);
    if(!pMol1) return false;

    PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
    if(ppReactant2)
    {
      pMol2 = GetMolRef(ppReactant2);
      if(!pMol2) return false;
      m_Reactant2 = pMol2;
    }
    //Put the Colliding Molecule into m_Reactant, even if it is second in datafile
    CollidingMolecule* pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
    if(pColMol)
    {
      m_Reactant = pColMol;
      m_Reactant2 = pMol2;
    }
    else
    {
      pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
      if(!pColMol)
      {
        stringstream errorMsg;
        errorMsg << "Either " << pMol1->getName() << " or "
                 << pMol2->getName() <<" has to be a modelled molecule";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obAuditMsg);
        return false;
      }
      m_Reactant = pColMol;
      m_Reactant2 = pMol1;
    }

    //Read products ... if any.
    pMol2=NULL;
    PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);

      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      pMol2 = GetMolRef(ppProduct2);

      //Put the Colliding Molecule into m_Product, even if it is second in datafile
      pColMol = dynamic_cast<CollidingMolecule*>(pMol1);
      if(pColMol)
      {
        m_Product = pColMol;
        m_Product2 = pMol2;
      }
      else
      {
        pColMol = dynamic_cast<CollidingMolecule*>(pMol2);
        if(!pColMol)
        {
          stringstream errorMsg;
          errorMsg << "Either " << pMol1->getName() << " or "
                   << pMol2->getName() <<" has to be a modelled molecule";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obAuditMsg);
          return false;
        }
        m_Product = pColMol;
        m_Product2 = pMol1;
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

    //Read in Kinetic rate parameters, if present
    m_ActEne = std::numeric_limits<double>::quiet_NaN();//means not set

    const char* pActEnetxt = ppReac->XmlReadValue("me:activationEnergy",false);
    const char* pPreExptxt = ppReac->XmlReadValue("me:preExponential",false);
    if (pActEnetxt && pPreExptxt)
    {
      stringstream s1(pActEnetxt);
      s1 >> m_ActEne ;
      stringstream s2(pPreExptxt);
      s2 >> m_PreExp ;
    }

    // Classify the reaction.

    if (m_Reactant && m_Product && !m_Reactant2 && !m_Product2)
        m_reactiontype = ISOMERIZATION ;
    else if (m_Reactant && m_Product && m_Reactant2 && !m_Product2)
        m_reactiontype = ASSOCIATION ;
    else if (m_Reactant)
        m_reactiontype = DISSOCIATION ;
    else {
      m_reactiontype = ERROR_REACTION ;
      stringstream errorMsg;
      errorMsg << "Unknown combination of reactants and products";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obAuditMsg);
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
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        return false;
      }
    }
    return true;
  }


  Molecule* Reaction::GetMolRef(PersistPtr pp)
  {
    Molecule* pMol;
    if(!pp)
        return NULL;
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) return false;
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
  void Reaction::get_unimolecularspecies(vector<CollidingMolecule *> &unimolecularspecies) const
  {
    if(m_Reactant2 == NULL){                // Possible dissociation or isomerization.
        unimolecularspecies.push_back(m_Reactant) ;
    }

    if(m_Product && m_Product2 == NULL){  // Possible association or isomerization.
        unimolecularspecies.push_back(m_Product) ;
    }
  }

  //
  // Get the prinincipal source reactant (i.e. reactant not in excess) if it exists.
  // (Not sure if this is a good idea, may be better to pass a Map in.)
  //
  CollidingMolecule *Reaction::get_pseudoIsomer() const
  {
    CollidingMolecule *pseudoIsomer = NULL ;
    if(m_reactiontype == ASSOCIATION) pseudoIsomer = m_Product ;
    return pseudoIsomer ;
  }

  //
  // Access microcanoincal rate coefficients.
  //
  void Reaction::get_MicroRateCoeffs(std::vector<double> &kmc) {
    calcGrnAvrgMicroRateCoeffs();
    kmc = m_GrainKfmc ;
  }

  //
  // Calculate grain averaged microcanonical rate coefficients.
  //
  bool Reaction::calcGrnAvrgMicroRateCoeffs() {

    // Calculate microcanonical rate coefficients.
    if (m_CellKfmc.size()==0)
    {
      if(!m_pMicroRateCalculator->calculateMicroRateCoeffs(this, m_CellKfmc) ||
        (GetSys()->TestMicroRatesEnabled() && !m_pMicroRateCalculator->testMicroRateCoeffs(this, m_CellKfmc, m_ppPersist)))
        return false;
    }
    // Calculate Grain averages of microcanonical rate coefficients.
    if (m_GrainKfmc.size()==0)
        return grnAvrgMicroRateCoeffs() ;
    return true;
  }

  //
  // Access microcanonical rate coefficients - cell values are averaged
  // to give grain values. This code is similar to that in Molecule.cpp
  // and this averaging should be done there. SHR 19/Sep/2004.
  //
  bool Reaction::grnAvrgMicroRateCoeffs() {

    int ngrn = GetSys()->MAXGrn();
    int currentGrainSize = GetSys()->getGrainSize();
    m_GrainKfmc.resize(ngrn);

    // Extract density of states of equilibrium molecule.

    vector<double> cellDOS(GetSys()->MAXCell(),0.0) ;
    m_Reactant->cellDensityOfStates(cellDOS) ;

    // Check that there are enough cells.

    if (currentGrainSize < 1) {
      stringstream errorMsg;
      errorMsg << "Not enought Cells to produce requested number of Grains.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i = 0; i < ngrn ; ++i ) {

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
    if ( idx2 < ngrn ) {
      stringstream errorMsg;
      errorMsg << "Number of grains produced is less than requested" << endl
               << "Number of grains requested: " << ngrn << endl
               << "Number of grains produced : " << idx2 << ".";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    return true;
  }

  //
  // Add microcanonical terms to collision operator
  //
  void Reaction::AddMicroRates(dMatrix *CollOptr, isomerMap &isomermap, sourceMap &sourcemap, const double rMeanOmega) {

    // Calculate Microcanonical rate coefficients.

    calcGrnAvrgMicroRateCoeffs() ;

    // Add microcanonical rates to the collision operator.

    switch(m_reactiontype) {
      case ISOMERIZATION :
        AddIsomerReactionTerms(CollOptr, isomermap, rMeanOmega) ;
        break;

      case ASSOCIATION :
        AddAssocReactionTerms(CollOptr, isomermap, sourcemap, rMeanOmega) ;
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
  void Reaction::AddIsomerReactionTerms(dMatrix *CollOptr,
                                        isomerMap &isomermap,
                                        const double rMeanOmega)
  {
    // Locate isomers in system matrix.

    const int rctLocation = isomermap[m_Reactant] ;
    const int pdtLocation = isomermap[m_Product] ;

    // Get densities of states for detailed balance.

    const int ngrn = GetSys()->MAXGrn();
    vector<double> rctDOS(ngrn, 0.0) ;
    vector<double> pdtDOS(ngrn, 0.0) ;

    m_Reactant->grnDensityOfStates(rctDOS) ;
    m_Product->grnDensityOfStates(pdtDOS) ;

    const int idx = m_Product->get_grnZpe() - m_Reactant->get_grnZpe() ;
    for ( int i = max(0,-idx) ; i < min(ngrn,(ngrn-idx)) ; ++i ) {
      int ll = i + idx ;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation + i) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[ll] ;                            // Forward loss reaction.
      (*CollOptr)[jj][jj] -= rMeanOmega * m_GrainKfmc[ll]*rctDOS[ll]/pdtDOS[i] ;       // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(rctDOS[ll]/pdtDOS[i]) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                 // Reactive gain.
    }
  }

  //
  // Add (reversible) association reaction terms to collision matrix.
  //
  void Reaction::AddAssocReactionTerms(dMatrix      *CollOptr,
                                       isomerMap    &isomermap,
                                       sourceMap    &sourcemap,
                                       const double rMeanOmega)
  {
    // Locate isomers in system matrix.

    const int rctLocation = isomermap[m_Reactant] ;
    const int pdtLocation = sourcemap[m_Product] ;

    // Get densities of states for detailed balance.

    const int ngrn = GetSys()->MAXGrn();
    vector<double> rctDOS(ngrn, 0.0) ;

    m_Reactant->grnDensityOfStates(rctDOS) ;

    const int idx = m_Product->get_grnZpe() - m_Reactant->get_grnZpe() ;
    for ( int i = max(0,-idx) ; i < min(ngrn,(ngrn-idx)) ; ++i ) {
      int ll = i + idx ;
      int ii(rctLocation + ll) ;
      int jj(pdtLocation) ;
      (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[ll] ;                            // Forward loss reaction.
      //(*CollOptr)[jj][jj] -= rMeanOmega * m_GrainKfmc[ll]*rctDOS[ll]/pdtDOS[i] ;       // Backward loss reaction from detailed balance.
      //(*CollOptr)[ii][jj]  = rMeanOmega * m_GrainKfmc[ll]*sqrt(rctDOS[ll]/pdtDOS[i]) ; // Reactive gain.
      //(*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                                 // Reactive gain.
    }
  }

  //
  // Add dissociation reaction terms to collision matrix.
  //
  void Reaction::AddDissocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) {

      // Locate reactant in system matrix.

      const int rctLocation = isomermap[m_Reactant] ;

      for ( int i = 0 ; i < m_Reactant->get_colloptrsize() ; ++i ) {
          int ii(rctLocation + i) ;
          (*CollOptr)[ii][ii] -= rMeanOmega * m_GrainKfmc[i] ;                            // Forward loss reaction.
      }
  }

}//namespace
