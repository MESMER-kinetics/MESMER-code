//-------------------------------------------------------------------------------------------
//
// ReactionManager.cpp
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This file contains the implementation of the ReactionManager class.
//
//-------------------------------------------------------------------------------------------
#include "ReactionManager.h"

using namespace Constants ;
using namespace std ;

namespace mesmer
{
  ReactionManager::ReactionManager(MoleculeManager *pMoleculeManager)
    :m_reactions(),
    m_pMoleculeManager(pMoleculeManager),
    m_reactionOperator(0),
    m_eigenvectors(0),
    m_eigenvalues(),
    m_basisMatrix(0),
    m_reducedBasisMatrix(0),
    m_reducedEigenvectors(0),
    m_reducedEigenvalues(),
    m_eqVector(),
    m_reducedEqVector(),
    m_locSizeMap(),
    m_divMap(),
    m_isomers(),
    m_sources(),
    m_sinkRxns(),
    m_SpeciesSequence(),
    m_meanOmega(0.0),
    punchSymbolGathered(false)
  {};

  //
  // Add a new reaction to the map.
  //
  bool ReactionManager::addreactions(PersistPtr ppReacList, const MesmerEnv& mEnv, MesmerFlags& mFlags)
  {
    bool readStatus(true);
    PersistPtr ppReac = ppReacList;
    while(ppReac = ppReac->XmlMoveTo("reaction"))
    {
      //ignore reactions with attribute active="false"
      const char* active = ppReac->XmlReadValue("active", optional);
      if(active && !strcmp(active, "false"))
        continue;

      //Read reaction ID
      const char* id = ppReac->XmlReadValue("id", false);
      if(!id){
        cinfo << "Reaction ID not found.\n";
        return false;
      }
      ErrorContext c(id);
      cinfo << "Parsing reaction..." << endl;

      // Read reactant and product types.

      string rct1Name, rct1Type, rct2Name, rct2Type ;
      string pdt1Name, pdt1Type, pdt2Name, pdt2Type ;
      bool bRct2(false), bPdt1(false), bPdt2(false) ;

      PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
      if(!ppReactantList)
        ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

      PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
      if(ppReactant1) {
        readStatus = GetMoleculeInfo(ppReactant1, rct1Name, rct1Type) ;

        PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
        if(ppReactant2) {
          readStatus = (readStatus && GetMoleculeInfo(ppReactant2, rct2Name, rct2Type)) ;
          bRct2 = true;
        }
      }

      PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
      if(!ppProductList)
        ppProductList=ppReac; //Be forgiving; we can get by without a productList element

      PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
      if (ppProduct1) {
        readStatus = (readStatus && GetMoleculeInfo(ppProduct1, pdt1Name, pdt1Type)) ;
        bPdt1 = true ;

        PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
        if (ppProduct2){
          readStatus = (readStatus && GetMoleculeInfo(ppProduct2, pdt2Name, pdt2Type)) ;
          bPdt2 = true ;
        }
      }

      if (!readStatus)
        return false ;

      //
      // Create a new Reaction.  For association & exchange reactions, if rct1Type == reactant,
      // bool is true and rct1 is the pseudoisomer.  if not, bool is false, and rct1 is the excess
      //
      Reaction *preaction ;
      if     (!bRct2 && bPdt1 && pdt1Type == "modelled" && !bPdt2)
        preaction = new IsomerizationReaction(m_pMoleculeManager, mEnv, mFlags, id) ;
      else if( bRct2 && bPdt1 && !bPdt2)
        preaction = new AssociationReaction(m_pMoleculeManager, mEnv, mFlags, id, (rct1Type == "deficientReactant")) ;
      else if(!bRct2 && bPdt1 && (pdt1Type == "sink" || pdt2Type == "sink"))
        preaction = new IrreversibleUnimolecularReaction(m_pMoleculeManager, mEnv, mFlags, id) ;
      else if( bRct2 && bPdt1 && (pdt1Type == "sink" || pdt2Type == "sink"))
        preaction = new IrreversibleExchangeReaction(m_pMoleculeManager, mEnv, mFlags, id, (rct1Type == "deficientReactant")) ;
      else {
        cinfo << "Unknown reaction type.\n";
        return false ;
      }

      // The information of the products of a dissociation reaction is necessary, as in
      // the xml output, Mesmer needs to know the products to draw the potential energy
      // surface. In addition, for dissociation reaction with QM tunneling, Mesmer also
      // needs to know the barrier height on the products side.

      //
      // Initialize Reaction from input stream.
      //
      readStatus = preaction->InitializeReaction(ppReac);
      if(!readStatus){
        //delete preaction;
        //return false;
        cerr << "UNSATISFACTORY REACTION\n"; //but keep parsing
      }

      //
      // Add reaction to map.
      //

      //need to check if there is duplicate reaction name/species: CHL

      m_reactions.push_back(preaction) ;
    }

    return readStatus;
  }

  void ReactionManager::resetCalcFlags(){
    for (size_t i(0) ; i < size() ; ++i) {
      m_reactions[i]->resetCalcFlag();
    }
  }

  bool ReactionManager::SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne)
  {
    //  Grain size and number of grain:
    //
    //  - Either grain size or number of grains can be specified, but not both.
    //
    //  - Uses the value of grain size in the datafile, if specified.
    //
    //  - If grain size is not specified but number of grains is, use a grain size to fit the energy range.
    //  If neither is specified, the grain size is set to 100cm-1 and the number of grains set so that
    //  the energy range is sufficient.
    //
    //  Energy Range:
    //
    //  - The required total energy domain extends from the lowest zero point energy of the lowest molecule
    //  to 10 k_B T above the highest.

    mEnv.EMin = minEne;
    mEnv.EMax = maxEne;

    /*For testing purposes, set the maxGrn based on the highest temperature we use in all calculations.*/
    const double MaximumTemperature = mEnv.MaximumTemperature;

    /*EAboveHill: Max energy above the highest hill. The temperature refers to the current condition.*/
    if (mFlags.useTheSameCellNumber){
      mEnv.EMax += mEnv.EAboveHill * MaximumTemperature * boltzmann_RCpK;
    }
    else{
      mEnv.EMax += mEnv.EAboveHill / mEnv.beta;
    }

    if(mEnv.GrainSize <= 0.0){
      mEnv.GrainSize = 100; //default 100cm-1
      cerr << "Grain size was invalid. Reset grain size to default: 100";
    }

    mEnv.MaxGrn = (int)((mEnv.EMax-mEnv.EMin)/mEnv.GrainSize + 0.5);
    mEnv.MaxCell = mEnv.GrainSize * mEnv.MaxGrn;

    //clog << "Cell number = " << mEnv.MaxCell << ", Grain number = " << mEnv.MaxGrn << endl;
    cinfo << "Cell number = " << mEnv.MaxCell << ", Grain number = " << mEnv.MaxGrn << endl;

    return true;
  }

  bool ReactionManager::BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags)
  {
    // reset the DOS calculation flags before building the reaction operator
    resetCalcFlags();

    const double SUPREMUM =  9e23 ;
    const double INFIMUM  = -SUPREMUM ;
    //
    // Find all the unique wells and lowest zero point energy.
    //
    m_isomers.clear();

    double minEnergy(SUPREMUM) ; // The minimum & maximum ZPE amongst all wells, set artificially large and small
    double maxEnergy(INFIMUM) ;  // to guarantee that each is overwritten in setting minEnergy and maxEnergy.
    Molecule *pBathGasMolecule = m_pMoleculeManager->get_BathGasMolecule();

    // populate molMapType with unimolecular species and determine minimum/maximum energy on the PES
    for (size_t i(0) ; i < size() ; ++i) {
      double TS_ZPE(INFIMUM);

      // Transition State
      // third check for the transition state in this reaction
      Molecule *pTransitionState = m_reactions[i]->get_TransitionState();
      if (pTransitionState){
        TS_ZPE = pTransitionState->getDOS().get_zpe();
        maxEnergy = max(maxEnergy, TS_ZPE) ;
      }

      // unimolecular species
      vector<Molecule *> unimolecules ;
      m_reactions[i]->get_unimolecularspecies(unimolecules) ;
      // populate molMapType with unimolecular species
      for (size_t j(0) ; j < unimolecules.size() ; ++j) {
        // wells
        Molecule *pCollidingMolecule = unimolecules[j] ;
        const double collidingMolZPE(pCollidingMolecule->getDOS().get_zpe());
        if(pCollidingMolecule && m_isomers.find(pCollidingMolecule) == m_isomers.end()){ // New isomer
          m_isomers[pCollidingMolecule] = 0 ; //initialize to a trivial location

          minEnergy = min(minEnergy, collidingMolZPE) ;
          maxEnergy = max(maxEnergy, collidingMolZPE) ;
        }

        //calculate the lowest barrier associated with this well(species)
        if (TS_ZPE != INFIMUM){
          const double barrierHeight = TS_ZPE - collidingMolZPE;
          if (barrierHeight < pCollidingMolecule->getColl().getLowestBarrier()){
            pCollidingMolecule->getColl().setLowestBarrier(barrierHeight);
          }
        }
      }

      //
      // For Association reactions determine zero point energy location of the
      // associating pair.
      //
      AssociationReaction *pReaction = dynamic_cast<AssociationReaction*>(m_reactions[i]) ;
      if (pReaction) {
        double pseudoIsomerZPE = pReaction->get_pseudoIsomer()->getDOS().get_zpe();
        double excessReactantZPE = pReaction->get_excessReactant()->getDOS().get_zpe();
        double sourceTermZPE = pseudoIsomerZPE + excessReactantZPE;
        minEnergy = min(minEnergy, sourceTermZPE) ;
        maxEnergy = max(maxEnergy, sourceTermZPE) ;

        // Calculate the lowest barrier associated with this well(species)
        // For association reaction, it is assumed that the barrier height is close to the source term energy
        // and in a sense, it is preferable to set this variable to the source term energy even there is an explicit
        // transition state.
        double adductZPE = unimolecules[0]->getDOS().get_zpe();
        double barrierHeight = sourceTermZPE - adductZPE;
        if (barrierHeight < unimolecules[0]->getColl().getLowestBarrier()){
          unimolecules[0]->getColl().setLowestBarrier(barrierHeight);
        }
      }

      //
      // For dissociation reactions determine zero point energy location of the
      // dissociation pair.
      //
      IrreversibleUnimolecularReaction *pDissnRtn = dynamic_cast<IrreversibleUnimolecularReaction*>(m_reactions[i]) ;
      if (pDissnRtn) {
        double thresholdEstimate = pDissnRtn->get_ThresholdEnergy();
        minEnergy = min(minEnergy, thresholdEstimate) ;
        maxEnergy = max(maxEnergy, thresholdEstimate) ;

        // Calculate the lowest barrier associated with this well(species).
        if (thresholdEstimate < unimolecules[0]->getColl().getLowestBarrier()){
          unimolecules[0]->getColl().setLowestBarrier(thresholdEstimate);
        }
      }

    }

    // set grain parameters for the current Temperature/pressure condition
    if(!SetGrainParams(mEnv, mFlags, minEnergy, maxEnergy))
      return false;

    // Calculate flux and k(E)s
    for (size_t i(0) ; i < size() ; ++i) {
      m_reactions[i]->calcGrnAvrgMicroRateCoeffs() ;
    }

    if (!mFlags.rateCoefficientsOnly){
      //
      // Shift all wells to the same origin, calculate the size of the reaction operator,
      // calculate the mean collision frequency and initialize all collision operators.
      //
      int msize(0) ; // size of the collision matrix
      m_meanOmega = 0.0;

      Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
      for (; isomeritr != m_isomers.end() ; ++isomeritr) {

        Molecule *isomer = isomeritr->first ;
        isomeritr->second = msize ; //set location

        int grnZpe = isomer->getColl().get_grnZPE() ; //set grain ZPE (with respect to the minimum of all wells)

        int colloptrsize = mEnv.MaxGrn - grnZpe ;
        isomer->getColl().set_colloptrsize(colloptrsize) ;
        msize += colloptrsize ;

        if(!isomer->getColl().initCollisionOperator(mEnv.beta, pBathGasMolecule)){
          cerr << "Failed initializing collision operator for " << isomer->getName();
          return false;
        }

        // update the size of the collision operator if it is different.
        int nGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
        if (nGroupedGrains != 0){
          msize -= (nGroupedGrains - 1);
        }

        m_meanOmega += isomer->getColl().get_collisionFrequency() ;
      }
      m_meanOmega /= double(m_isomers.size());

      //
      // Find all source terms.
      // Note: 1. A source term contains the deficient reactant.  It is possible for
      //          there to be more than one source term.
      m_sources.clear(); // Maps the location of source in the system matrix.
      for (size_t i(0) ; i < size() ; ++i) {
        AssociationReaction *pReaction = dynamic_cast<AssociationReaction*>(m_reactions[i]) ;
        if (pReaction) {  // if the reaction is an association reaction
          Molecule *pPseudoIsomer = pReaction->get_pseudoIsomer();
          if (pPseudoIsomer && m_sources.find(pPseudoIsomer) == m_sources.end()){ // reaction includes a new
            m_sources[pPseudoIsomer] = msize ;                                    // pseudoisomer
            ++msize ;
          }
          pReaction->putSourceMap(&m_sources) ;
        }
        IrreversibleExchangeReaction *IREreaction = dynamic_cast<IrreversibleExchangeReaction*>(m_reactions[i]);
        if(IREreaction){       // if the reaction is an irreversible exchange reaction
          Molecule *pPseudoIsomer = IREreaction->get_pseudoIsomer();
          if(pPseudoIsomer && m_sources.find(pPseudoIsomer) == m_sources.end()){  // reaction includes a new
            m_sources[pPseudoIsomer] = msize ;                                    // pseudoisomer
            ++msize;
          }
          IREreaction->putSourceMap(&m_sources);
        }
      }

      // Build reaction operator.
      //
      // One of two methods for building the reaction operator are available:
      // the conventional energy grained master equation method which is based
      // on energy grains and a contracted basis set method in which a basis
      // set is generated from the individual collision operators and a
      // representation of the reaction operator build upon this basis.

      if (!mFlags.doBasisSetMethod) {

        // Full energy grained reaction operator.

        constructGrainMatrix(msize);

      } else {

        // Contracted basis set reaction operator.

        constructBasisMatrix();

      }
    }

    return true;
  }

  void ReactionManager::printReactionOperator(const MesmerFlags& mFlags)
  {
    const int smsize = int(m_reactionOperator->size()) ;

    switch (mFlags.printReactionOperatorNum)
    {
    case -1:
      ctest << "Printing all (" << smsize << ") columns/rows of the Reaction Operator:\n";
      (*m_reactionOperator).showFinalBits(smsize, mFlags.print_TabbedMatrices);
      break;
    case -2:
      ctest << "Printing final 1/2 (" << smsize/2 << ") columns/rows of the Reaction Operator:\n";
      (*m_reactionOperator).showFinalBits(smsize/2, mFlags.print_TabbedMatrices);
      break;
    case -3:
      ctest << "Printing final 1/3 (" << smsize/3 << ") columns/rows of the Reaction Operator:\n";
      (*m_reactionOperator).showFinalBits(smsize/3, mFlags.print_TabbedMatrices);
      break;
    default: // the number is either smaller than -3 or positive
      if (abs(mFlags.printReactionOperatorNum) > smsize){
        ctest << "Printing all (" << smsize << ") columns/rows of the Reaction Operator:\n";
        (*m_reactionOperator).showFinalBits(smsize, mFlags.print_TabbedMatrices);
      }
      else{
        ctest << "Printing final " << abs(mFlags.printReactionOperatorNum) << " columns/rows of the Reaction Operator:\n";
        (*m_reactionOperator).showFinalBits(abs(mFlags.printReactionOperatorNum), mFlags.print_TabbedMatrices);
      }
    }
  }

  void ReactionManager::printEigenvectors(const MesmerFlags& mFlags)
  {
    const int smsize = int(m_eigenvectors->size()) ;

    switch (mFlags.printReactionOperatorNum)
    {
    case -1:
      ctest << "Printing all (" << smsize << ") columns/rows of the eigenvectors:\n";
      (*m_eigenvectors).showFinalBits(smsize, mFlags.print_TabbedMatrices);
      break;
    case -2:
      ctest << "Printing final 1/2 (" << smsize/2 << ") columns/rows of the eigenvectors:\n";
      (*m_eigenvectors).showFinalBits(smsize/2, mFlags.print_TabbedMatrices);
      break;
    case -3:
      ctest << "Printing final 1/3 (" << smsize/3 << ") columns/rows of the eigenvectors:\n";
      (*m_eigenvectors).showFinalBits(smsize/3, mFlags.print_TabbedMatrices);
      break;
    default: // the number is either smaller than -3 or positive
      if (abs(mFlags.printReactionOperatorNum) > smsize){
        ctest << "Printing all (" << smsize << ") columns/rows of the eigenvectors:\n";
        (*m_eigenvectors).showFinalBits(smsize, mFlags.print_TabbedMatrices);
      }
      else{
        ctest << "Printing final " << abs(mFlags.printReactionOperatorNum) << " columns/rows of the eigenvectors:\n";
        (*m_eigenvectors).showFinalBits(abs(mFlags.printReactionOperatorNum), mFlags.print_TabbedMatrices);
      }
    }
  }

  void ReactionManager::diagReactionOperator(const MesmerFlags &mFlags, const int precision)
  {
    // Allocate space for eigenvalues.
    const int smsize = int(m_reactionOperator->size()) ;
    m_eigenvalues.clear();
    m_eigenvalues.resize(smsize, 0.0);
    if (m_eigenvectors) delete m_eigenvectors;
    m_eigenvectors = new qdMatrix(smsize, 0.0) ;

    // This block prints Reaction Operator before diagonalization
    if (mFlags.printReactionOperatorNum){
      ctest << "Reaction operator --- ";
      printReactionOperator(mFlags);
    }

    //-------------------------------------------------------------
    // diagonalize the whole matrix
    switch (precision){
      case 0: // diagonalize in double
        {
          dMatrix dDiagM(smsize);
          for ( int i = 0 ; i < smsize ; ++i )
            for ( int j = 0 ; j < smsize ; ++j )
              dDiagM[i][j] = to_double((*m_reactionOperator)[i][j]) ;
          vector<double>  dEigenValue(smsize, 0.0);
          dDiagM.diagonalize(&dEigenValue[0]) ;
          for ( int i = 0 ; i < smsize ; ++i )
            m_eigenvalues[i] = dEigenValue[i];
          for ( int i = 0 ; i < smsize ; ++i )
            for ( int j = 0 ; j < smsize ; ++j )
              (*m_eigenvectors)[i][j] = dDiagM[i][j] ;
          break;
        }
      case 1: // diagonalize in double double
        {
          ddMatrix ddDiagM(smsize);
          for ( int i = 0 ; i < smsize ; ++i )
            for ( int j = 0 ; j < smsize ; ++j )
              ddDiagM[i][j] = to_dd_real((*m_reactionOperator)[i][j]) ;
          vector<dd_real> ddEigenValue(smsize, 0.0);
          ddDiagM.diagonalize(&ddEigenValue[0]) ;
          for ( int i = 0 ; i < smsize ; ++i )
            m_eigenvalues[i] = ddEigenValue[i];
          for ( int i = 0 ; i < smsize ; ++i )
            for ( int j = 0 ; j < smsize ; ++j )
              (*m_eigenvectors)[i][j] = ddDiagM[i][j] ;
          break;
        }
      default: // diagonalize in quad double
        {
          qdMatrix qdDiagM(smsize);
          for ( int i = 0 ; i < smsize ; ++i )
            for ( int j = 0 ; j < smsize ; ++j )
              qdDiagM[i][j] = (*m_reactionOperator)[i][j] ;
          qdDiagM.diagonalize(&m_eigenvalues[0]) ;
          for ( int i = 0 ; i < smsize ; ++i )
            for ( int j = 0 ; j < smsize ; ++j )
              (*m_eigenvectors)[i][j] = qdDiagM[i][j] ;
        }
    }
    // diagonalize the whole matrix
    //-------------------------------------------------------------

    // This block prints Eigenvectors
    if (mFlags.printReactionOperatorNum){
      ctest << "Eigenvectors --- ";
      printEigenvectors(mFlags);
    }

    int numberStarted = 0;
    int numberPrinted = smsize; // Default prints all of the eigenvalues
    if (mFlags.printEigenValuesNum > 0 && mFlags.printEigenValuesNum <= smsize){ //at least prints 1 eigenvalue
      numberPrinted = mFlags.printEigenValuesNum;
      numberStarted = smsize - mFlags.printEigenValuesNum;
    }

    ctest << "\nTotal number of eigenvalues = " << smsize << endl;
    ctest << "Eigenvalues\n{\n";
    for (int i = numberStarted ; i < smsize; ++i) {
      qd_real tmp = (mFlags.doBasisSetMethod)? m_eigenvalues[i] : m_eigenvalues[i] * m_meanOmega ;
      formatFloat(ctest, tmp, 6, 15) ;
      ctest << endl ;
    }
    ctest << "}\n";
  }

  //
  // Extract molecule information from XML stream.
  //
  bool ReactionManager::GetMoleculeInfo(PersistPtr pp, string& MolName, string& MolType)
  {
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) {
      cerr << "Ill formed molecule tag." << endl;
      return false;
    }
    const char *pmolname, *pmoltype;
    pmolname = ppmol->XmlReadValue("ref", false);
    if(pmolname) {
      ErrorContext c(pmolname);
      MolName = pmolname;
      pmoltype = ppmol->XmlReadValue("me:type");
      if(pmoltype)
        MolType = pmoltype;
    }
    return pmolname && pmoltype;
  }

  bool ReactionManager::calculateEquilibriumFractions(const double beta)
  { /* Consider a three well system: e.g., A <-> B <-> C where A <-> B has Keq = K1 & B <-> C has Keq = K2.
    This routine uses the fact that the normalized equilibrated system may be described
    by a 3x3 matrix and a vector which satisfy the following:
    |-K1  1   0| |A|   |0|
    | 0  -K2  1| |B| = |0|
    | 1   1   1| |C|   |1|
    The equilibrium fraction of each isomer (or pseudo isomer, in the case of a source term) may be
    obtained by inverting the matrix shown above, and taking the elements in the final column of the inverse.
    Any system, with an arbitrary number of wells and connections, may be described by such a Matrix */

    // determine the total number of isomers + sources from the m_isomers and m_sources maps
    int eqMatrixSize = int(m_isomers.size() + m_sources.size());

    // intialize the matrix which holds the system of equations that describe the equilibrium distribution
    dMatrix  eqMatrix(eqMatrixSize);

    // initialize a map of equilibrium fractions
    m_SpeciesSequence.clear();

    // loop over the number of reactions in order to assign elements to the m_SpeciesSequence map
    // and then update the corresponding matrix elements in eqMatrix

    int counter(0);   //counter keeps track of how may elements are in the m_SpeciesSequence map
    for (size_t i(0) ; i < size() ; ++i) {  //iterate through m_reactions

      Molecule* rct;
      Molecule* pdt;
      double Keq(0.0);

      //only need eq fracs for species in isom & assoc rxns
      if (m_reactions[i]->isEquilibratingReaction(Keq, &rct, &pdt)){

        int ploc, rloc ;

        Reaction::molMapType::iterator rctitr = m_SpeciesSequence.find(rct);   //check if the reactant is in the map
        bool rval = (rctitr != m_SpeciesSequence.end()) ;       //if the reactant isnt in the map
        if (rval)
          rloc = rctitr->second ;        //if the reactant is in the map, get the location

        Reaction::molMapType::iterator pdtitr = m_SpeciesSequence.find(pdt);   //check if the product is in the map
        bool pval = (pdtitr != m_SpeciesSequence.end()) ;       //if the product isnt in the map
        if (pval)
          ploc = pdtitr->second;        //if the product is in the map, get the location

        if(!rval && !pval){             // if neither reactant nor product are in the m_SpeciesSequence map
          m_SpeciesSequence[rct] = counter;            // update the eqMatrix elements
          counter++ ;
          m_SpeciesSequence[pdt] = counter;
          eqMatrix[counter-1][counter-1] -= Keq;
          eqMatrix[counter-1][counter] += 1.0;
          counter++ ;
        }
        else if(!rval && pval){        // if reactant isnt in m_SpeciesSequence map & product is
          m_SpeciesSequence[rct] = counter;            // update the eqMatrix matrix elements
          eqMatrix[counter-1][ploc] += 1.0;
          eqMatrix[counter-1][counter] -= Keq;
          counter++ ;
        }
        else if(rval && !pval){        // if reactant is in m_SpeciesSequence map & product isnt
          m_SpeciesSequence[pdt] = counter;            // update the eqMatrix matrix elements
          eqMatrix[counter-1][rloc] -= Keq;
          eqMatrix[counter-1][counter] += 1.0 ;
          counter++ ;
        }
        else if(rval && pval){        // if both reactant & product are in m_SpeciesSequence map

          double pdtRowSum(0.0), rctRowSum(0.0);

          for(int j(0);j<counter;++j){           // calculate pdt & rct rowSums of EqMatrix to see if the rxn is redundant
            pdtRowSum += eqMatrix[ploc][j];
            rctRowSum += eqMatrix[rloc][j];
          }

          if(pdtRowSum!=0.0 && rctRowSum!=0.0){ // connection is redundant
            eqMatrix[counter-1][ploc] += 1.0 ;
            eqMatrix[counter-1][rloc] -= Keq ;
          }
          else if(rctRowSum==0.0){              // connection is not redundant, pdts lack specification
            eqMatrix[rloc][ploc] += 1.0 ;
            eqMatrix[rloc][rloc] -= Keq ;
          }
          else if(pdtRowSum==0.0){
            eqMatrix[ploc][ploc] += 1.0 ;        // connection is not redundant, rcts lack specification
            eqMatrix[ploc][rloc] -= Keq ;
          }
        }
      }
    }

    // if counter==0 after the for loop above, then there are no equilibrating reactions (i.e., all the reactions
    // are irreversible).  In that case, the lone isomer has an equilibrium fraction of 1.  Thus, we increment 
    // counter so that the 1 is added to the eqMatrix in the for loop immediately following
    if (counter==0){
      Molecule* rct=(m_isomers.begin())->first;
      m_SpeciesSequence[rct] = counter;
      ++counter;
    }

    for(int i=0; i < counter; ++i){         // add ones to the final row of the matrix
      eqMatrix[counter-1][i]= 1.0;
    }

    //    ctest << "matrix elements for calculating isomer equilibrium fractions:" << endl;
    //    eqMatrix.showFinalBits(counter);

    dMatrix backup(eqMatrix);  //backup EqMatrix for error reporting

    ctest << endl << "Eq fraction matrix:" << endl;
    backup.showFinalBits(counter);

    if(eqMatrix.invertGaussianJordan()){
      cerr << "Inversion of matrix for calculating Eq fractions failed.  Matrix before inversion is: ";
      backup.showFinalBits(counter);
    }

    ctest << "inverse of Eq fraction matrix:" << endl;
    eqMatrix.showFinalBits(counter);

    Reaction::molMapType::iterator itr1;

    for(itr1= m_SpeciesSequence.begin(); itr1!=m_SpeciesSequence.end(); ++itr1){  //assign Eq fraction to appropriate Molecule
      int seqMatrixLoc = itr1->second;                          //in the Eq frac map
      Molecule* key = itr1->first;
      key->getPop().setEqFraction(eqMatrix[seqMatrixLoc][counter-1]);    //set Eq fraction to last column in eqMatrix
      ctest << "Equilibrium Fraction for " << key->getName() << " = " << key->getPop().getEqFraction() << endl;
    }
    return true;
  }

  bool ReactionManager::produceInitialPopulationVector(vector<double>& n_0){

    double populationSum = 0.0;

    Reaction::molMapType::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
      Molecule* isomer = ipos->first;                        // to get isomer initial populations
      populationSum += isomer->getPop().getInitPopulation();
    }

    Reaction::molMapType::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map to get
      Molecule* source = spos->first;                         // source initial populations
      populationSum += source->getPop().getInitPopulation();
    }

    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      Molecule* isomer = ipos->first;                        // get initial population of each isomer
      double initFrac = isomer->getPop().getInitPopulation();
      if (initFrac != 0.0){                                           // if isomer initial populations are nonzero
        initFrac /= populationSum;                                    // normalize initial pop fraction
        int rxnMatrixLoc = ipos->second;
        const int colloptrsize = isomer->getColl().get_colloptrsize();
        const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
        vector<double> boltzFrac;
        isomer->getColl().normalizedGrnBoltzmannDistribution(boltzFrac, colloptrsize, numberGrouped);
        if (numberGrouped == 0){
          for (int i = 0; i < colloptrsize; ++i){
            n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
          }
        }
        else{
          for (int i = 0; i < colloptrsize - colloptrsize + 1; ++i){
            n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
          }
        }
      }
    }

    // if there is no source term and the populationSum is still zero, set population = 1.0 for the first isomer
    int sizeSource = static_cast<int>(m_sources.size());
    if (populationSum == 0. && sizeSource == 0){
      ipos = m_isomers.begin();
      Molecule* isomer = ipos->first;
      isomer->getPop().setInitPopulation(1.0); // set initial population for the first isomer
      double initFrac = isomer->getPop().getInitPopulation();
      cinfo << "No population was assigned, and there is no source term."  << endl
        << "Initialize a Boltzmann distribution in the first isomer." << endl;
      int rxnMatrixLoc = ipos->second;
      const int colloptrsize = isomer->getColl().get_colloptrsize();
      const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
      vector<double> boltzFrac;
      isomer->getColl().normalizedInitialDistribution(boltzFrac, colloptrsize, numberGrouped);
      if (numberGrouped == 0){
        for (int i = 0; i < colloptrsize; ++i){
          n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
        }
      }
      else{
        for (int i = 0; i < colloptrsize - colloptrsize + 1; ++i){
          n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
        }
      }
    }

    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      Molecule* source = spos->first;
      int rxnMatrixLoc = spos->second;
      if (populationSum == 0. && spos == m_sources.begin()){
        cinfo << "No population was assigned. Initialize the first source term to 1.0." << endl;
        n_0[rxnMatrixLoc] = 1.0;
      }else{
        double initFrac = source->getPop().getInitPopulation() / populationSum;
        n_0[rxnMatrixLoc] = initFrac;
      }
    }

    return true;
  }

  //
  // The vector produced by this function contains the sqrt of
  // the normalized equilibrium distribution.
  //
  bool ReactionManager::produceEquilibriumVector()
  {

    m_eqVector.clear();
    m_eqVector.resize(m_reactionOperator->size());

    Reaction::molMapType::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map to get
      Molecule* source = spos->first;                            // eq Fractions
      int rxnMatrixLoc = spos->second;
      double eqFrac = source->getPop().getEqFraction();
      m_eqVector[rxnMatrixLoc] = sqrt(eqFrac);
    }

    Reaction::molMapType::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
      Molecule* isomer = ipos->first;                        // to get eq Fractions
      int rxnMatrixLoc = ipos->second;
      double eqFrac = isomer->getPop().getEqFraction();
      const int colloptrsize = isomer->getColl().get_colloptrsize();
      const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
      vector<double> boltzFrac;
      isomer->getColl().normalizedGrnBoltzmannDistribution(boltzFrac, colloptrsize, numberGrouped);
      if (numberGrouped == 0) {
        for(int i(0);i<colloptrsize;++i){
          m_eqVector[rxnMatrixLoc + i]= sqrt(eqFrac * boltzFrac[i]);
        }
      }
      else{
        for(int i(0);i<colloptrsize - numberGrouped + 1;++i){
          m_eqVector[rxnMatrixLoc + i]= sqrt(eqFrac * boltzFrac[i]);
        }
      }
    }
    return true;
  }

  bool ReactionManager::timeEvolution(MesmerFlags& mFlags)
  {
    int smsize = int(m_eigenvectors->size());

    if (!produceEquilibriumVector()){
      cerr << "Calculation of equilibrium vector failed.";
      return false;
    }

    vector<double> n_0(smsize, 0.); // initial distribution
    if (!produceInitialPopulationVector(n_0)){
      cerr << "Calculation of initial conditions vector failed.";
      return false;
    }

    // |n_0> = F^(-1)*|n_0>
    for (int j = 0; j < smsize; ++j) {
      n_0[j] /= m_eqVector[j];
    }

    // |n_0> is the initial populations of the grains for all species
    // |n_t> = U exp(Lamda t) U^-1 |n_0>
    // |r_0> = U^-1 |n_0>
    vector<double> r_0(smsize, 0.);

    double maxEvoTime = 0.;
    // set the default maximum evolution time
    if (mFlags.maxEvolutionTime <= 0. || mFlags.maxEvolutionTime > 1.0e8)
      maxEvoTime = 1.0e8;
    else
      maxEvoTime = mFlags.maxEvolutionTime;

    // calculate the time points
    vector<double> timePoints;
    for (int i = 0; i <= 160; ++i){
      double time = pow(10., static_cast<double>(i) / 10. - 11.);
      if (time > maxEvoTime)
        break;
      timePoints.push_back(time);
    }

    //initialize dt vector for calculating product yields
    vector<double> dt(timePoints.size()-1,0.0);
    dt[0] = timePoints[0];
    for (int i = 1; i < int(dt.size()); ++i){
      dt[i] = timePoints[i] - timePoints[i-1];
    }


    dMatrix totalEigenVecs(smsize); // copy full eigenvectors of the system
    for ( int i = 0 ; i < smsize ; ++i )
      for ( int j = 0 ; j < smsize ; ++j )
        totalEigenVecs[i][j] = to_double((*m_eigenvectors)[i][j]);


    for (int i = 0; i < smsize; ++i) {
      double sum = 0.;
      for (int j = 0; j < smsize; ++j) {
        sum += n_0[j] * totalEigenVecs[j][i];
      }
      r_0[i] = sum;  // now |r_0> = V^(T)*|init> = U^(-1)*|n_0>
    }

    for (int i = 0; i < smsize; ++i) {
      double tmp = m_eqVector[i];
      for (int j = 0; j < smsize; ++j) {
        totalEigenVecs[i][j] *= tmp;
      }
    }

    const int maxTimeStep = int(dt.size());
    db2D grnProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
    vector<double> work2(smsize, 0.);

    for (int timestep = 0; timestep < maxTimeStep; ++timestep){
      double numColl = m_meanOmega * timePoints[timestep];
      for (int j = 0; j < smsize; ++j) {
        work2[j] = r_0[j] * exp(to_double(m_eigenvalues[j]) * numColl);
      } // now |wk2> = exp(Dt)*V^(T)*|init> = exp(Dt)*U^(-1)*|n_0>
      for (int j = 0; j < smsize; ++j) {
        double sum = 0.;
        for (int l = 0; l < smsize; ++l) {
          sum += work2[l] * totalEigenVecs[j][l];
        }
        grnProfile[j][timestep] = sum;
      } // now |grnProfile(t)> = |grnProfile(i)> = F*V*exp(Dt)*V^(T)*|init> = U*exp(Dt)*U^(-1)*|n_0>
    }

    //------------------------------
    // print grained species profile
    if (mFlags.grainedProfileEnabled) {
      ctest << "\nGrained species profile (the first row is time points in unit of second):\n{\n";
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        formatFloat(ctest, timePoints[timestep], 6,  15);
      }
      ctest << endl;
      for (int j = 0; j < smsize; ++j) {
        for (int timestep = 0; timestep < maxTimeStep; ++timestep){
          formatFloat(ctest, grnProfile[j][timestep], 6,  15);
        }
        ctest << endl;
      }
      ctest << "}\n";
    }
    //------------------------------

    ctest<<"mean collision frequency = " << m_meanOmega << "/s" << endl;

    vector<double> totalIsomerPop(maxTimeStep, 0.);
    vector<double> totalPdtPop(maxTimeStep, 0.);

    for(int timestep(0); timestep<maxTimeStep; ++timestep){
      for(int j(0);j<smsize;++j){
        totalIsomerPop[timestep] += grnProfile[j][timestep];
      }
      double popTime = totalIsomerPop[timestep];
      if (popTime > 1.0){
        popTime = 1.0; // correct some numerical error
        //totalIsomerPop[timestep] = 1.0; // Not very sure if we should cover up this numerical error entirely!!?
      }
      else if (popTime < 0.0){
        popTime = 0.0;
        //totalIsomerPop[timestep] = 0.0; // Not very sure if we should cover up this numerical error entirely!!?
      }
      totalPdtPop[timestep] = 1.0 - popTime;
    }

    int sinkpos(0);
    m_sinkRxns.clear();                      // Populate map: m_sinkRxns[Rxns] = location of rct
    for (size_t i(0) ; i < size() ; ++i) {
      bool Irreversible = (m_reactions[i])->isIrreversible() ;
      if (Irreversible && m_sinkRxns.find(m_reactions[i]) == m_sinkRxns.end()) {   // add an irreversible rxn to the map
        Reaction* reaction = m_reactions[i];
        Molecule* source = reaction->get_reactant();
        Molecule* isomer = source;
        if(isomer){
          int rxnMatrixLoc = m_isomers[isomer];
          m_sinkRxns[reaction] = rxnMatrixLoc;
          m_SinkSequence[reaction] = sinkpos;               // populate SinkSequence map with Irreversible Rxns
          ++sinkpos;
        }
        else if(source){
          int rxnMatrixLoc = m_sources[source];
          m_sinkRxns[reaction] = rxnMatrixLoc;
          m_SinkSequence[reaction] = sinkpos;
          ++sinkpos;
        }
      }
    }

    if (mFlags.speciesProfileEnabled){
      ctest << endl << "Print time dependent species and product profiles" << endl << "{" << endl;
      int numberOfSpecies = static_cast<int>(m_isomers.size() + m_sources.size() + m_sinkRxns.size());
      db2D speciesProfile(numberOfSpecies, maxTimeStep);
      int speciesProfileidx(0);

      ctest << setw(16) << "Timestep/s";

      Reaction::molMapType::iterator spos;
      for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map
        Molecule* source = spos->first ;                        // to get source profile vs t
        ctest << setw(16) << source->getName();
        int rxnMatrixLoc = spos->second;
        for (int timestep = 0; timestep < maxTimeStep; ++timestep){
          double gPf = grnProfile[rxnMatrixLoc][timestep];
          speciesProfile[speciesProfileidx][timestep] = gPf;
        }
        ++speciesProfileidx;
      }

      Reaction::molMapType::iterator ipos;
      for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
        Molecule* isomer = ipos->first;                        // to get isomer profile vs t
        ctest << setw(16) << isomer->getName();
        int rxnMatrixLoc = ipos->second;
        const int colloptrsize = isomer->getColl().get_colloptrsize();
        const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
        if (numberGrouped == 0){
          for (int timestep = 0; timestep < maxTimeStep; ++timestep){
            for(int i = 0; i < colloptrsize; ++i){
              speciesProfile[speciesProfileidx][timestep] += grnProfile[i+rxnMatrixLoc][timestep];
            }
          }
        }
        else{
          for (int timestep = 0; timestep < maxTimeStep; ++timestep){
            for(int i = 0; i < colloptrsize - numberGrouped + 1; ++i){
              speciesProfile[speciesProfileidx][timestep] += grnProfile[i+rxnMatrixLoc][timestep];
            }
          }
        }
        ++speciesProfileidx;
      }

      sinkMap::iterator pos;      // iterate through sink map to get product profile vs t
      int pdtProfileStartIdx = speciesProfileidx;
      for (pos = m_sinkRxns.begin(); pos != m_sinkRxns.end(); ++pos){
        vector<double> KofEs;                             // vector to hold sink k(E)s
        Reaction* sinkReaction = pos->first;
        const int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
        vector<Molecule*> pdts;                               // in the sink reaction
        sinkReaction->get_products(pdts);
        if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
          KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
          ctest << setw(11) << pdts[0]->getName()<< setw(5) << "(bim)";
        }
        else{   // if the collision operator size is >1, there are k(E)s for the irreversible loss
          KofEs = sinkReaction->get_GrainKfmc();          // assign sink k(E)s, the vector size == maxgrn
          ctest << setw(16) << pdts[0]->getName();
        }
        int rxnMatrixLoc = pos->second;                       // get sink location
        double TimeIntegratedProductPop(0.0);
        for (int timestep = 0; timestep < maxTimeStep; ++timestep){
          for(int i = 0; i < colloptrsize; ++i){
            speciesProfile[speciesProfileidx][timestep] += KofEs[i]*grnProfile[i+rxnMatrixLoc][timestep]*dt[timestep];
          }
          TimeIntegratedProductPop += speciesProfile[speciesProfileidx][timestep];
          speciesProfile[speciesProfileidx][timestep]= TimeIntegratedProductPop;
        }
        ++speciesProfileidx;
        KofEs.clear();
      }

      if (pdtProfileStartIdx < speciesProfileidx){
        for(int timestep = 0; timestep < maxTimeStep; ++timestep){    // normalize product profile to account for small
          double normConst(0.0);                          // numerical errors in TimeIntegratedProductPop
          double pdtYield(0.0);
          for(int i(pdtProfileStartIdx); i<speciesProfileidx; ++i){   // calculate normalization constant
            pdtYield += speciesProfile[i][timestep];
          }
          normConst = totalPdtPop[timestep] / pdtYield;
          for(int i(pdtProfileStartIdx); i<speciesProfileidx; ++i){   // apply normalization constant
            speciesProfile[i][timestep] *= normConst;
          }
        }
      }

      ctest << setw(16)<< "totalIsomerPop" << setw(16)<< "totalPdtPop"  << endl;
      for(int timestep = 0; timestep < maxTimeStep; ++timestep){
        ctest << setw(16) << timePoints[timestep];
        for(int i(0); i<speciesProfileidx; ++i){
          ctest << setw(16) << speciesProfile[i][timestep];
        }
        ctest << setw(16) << totalIsomerPop[timestep] << setw(16) << totalPdtPop[timestep] << endl;
      }
      ctest << "}" << endl;
    }
    return true;
  }

  bool ReactionManager::BartisWidomPhenomenologicalRates(dMatrix& mesmerRates, MesmerFlags& mFlags, PersistPtr ppList)
  {
    const int smsize = int(m_eigenvectors->size());
    dMatrix eigenVec(smsize);  //copy ReactionOperator, the eigenvector Matrix (== V)

    for ( int i = 0 ; i < smsize ; ++i )
      for ( int j = 0 ; j < smsize ; ++j )
        eigenVec[i][j] = to_double((*m_eigenvectors)[i][j]) ;

    // constant variables
    const int nchem = static_cast<int>(m_isomers.size() + m_sources.size());  // number of isomers+pseudoisomers
    const int nchemIdx = smsize - nchem;       // idx for chemically significant eigenvalues & vectors

    dMatrix assymInvEigenVec(smsize);   // U^(-1)
    dMatrix assymEigenVec(smsize);      // U
    for(int i(0);i<smsize;++i){
      double tmp = m_eqVector[i];
      for(int j(0);j<smsize;++j){
        assymInvEigenVec[j][i] = eigenVec[i][j]/tmp;         //calculation of U^(-1) = (FV)^-1 = V^T * F^-1
        assymEigenVec[j][i] = m_eqVector[j] * eigenVec[j][i];  //calculation of U = FV
      }
    }

    //------------------------- TEST block ----------------------------------------
    dMatrix EigenVecIdentity(smsize);   // matrix for holding product of U^(-1) * U
    for(int i(0);i<smsize;++i){         // multiply U*U^(-1) for testing
      double test = 0.0;
      for(int j(0);j<smsize;++j){
        double sm = 0.0;
        for(int k(0);k<smsize;++k){
          sm += assymEigenVec[i][k] * assymInvEigenVec[k][j];
        }
        EigenVecIdentity[i][j] = sm;
        test += sm;
      }
      if((test/1.0) < 0.999 || (test/1.0) > 1.001)      // test that U*U^(-1) = 1
        ctest << "row " << i << " of the U*U^(-1) matrix does not equal unity. It sums to " << test << endl;
    }
    //------------------------- TEST block ----------------------------------------

    // EigenVecIdentity.showFinalBits(nchem);

    dMatrix Z_matrix(nchem);  // definitions of Y_matrix and Z_matrix taken from PCCP 2007(9), p.4085
    db2D Y_matrix;
    Reaction::molMapType::iterator ipos;  // set up an iterator through the isomer map
    Reaction::molMapType::iterator spos;  // set up an iterator through the source map
    sinkMap::iterator sinkpos;           // set up an iterator through the irreversible rxn map

    ctest << "\nBartis Widom eigenvalue/eigenvector analysis\n";
    ctest << endl << "Number of sinks in this system: " << m_sinkRxns.size() << endl;

    if(m_sinkRxns.size()>0){
      ctest << "\nThere should be " << nchem << " chemically significant eigenvalues (CSEs)" << endl;
    }
    else{
      ctest << "\nThere should be 1 zero eigenvalue (zero within numerical precision) and " << nchem-1
        << " chemically significant eigenvalues (CSEs)" << endl;
    }

    // check the separation between chemically significant eigenvalues (CSEs)
    // and internal energy relaxation eigenvalues (IEREs); if it's not good, print a warning

    const double last_CSE=(to_double(m_eigenvalues[nchemIdx]))* m_meanOmega;
    const double first_IERE=(to_double(m_eigenvalues[nchemIdx-1]))* m_meanOmega;
    const double CSE_IERE_separation=last_CSE/first_IERE;
    if(CSE_IERE_separation > 0.1){
      ctest << "\nWarning: CSEs not well separated from internal energy relaxation eigenvals (IEREs)" << endl;
      ctest << "\nThe last CSE = " << last_CSE << " and the first IERE = " << first_IERE << endl; 
      ctest << "(last CSE)/(first IERE) ratio = " << CSE_IERE_separation << ", which is less than an order of magnitude" << endl;
      ctest << "\nResults obtained from Bartis Widom eigenvalue-vector analysis may be unreliable" << endl;
    }

    for(int i(0); i<nchem; ++i){
      for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // calculate Z_matrix matrix elements for
        double sm = 0.0;                                                // all isomers in the system
        Molecule* isomer = ipos->first;
        const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
        const int colloptrsize = (!numberGroupedGrains)
          ? isomer->getColl().get_colloptrsize()
          : isomer->getColl().get_colloptrsize() - numberGroupedGrains + 1;
        // get colloptrsize for isomer
        int rxnMatrixLoc = ipos->second;                                // get location for isomer in the rxn matrix
        int seqMatrixLoc = m_SpeciesSequence[isomer];                   // get sequence position for isomer
        for(int j(0);j<colloptrsize;++j){
          sm += assymEigenVec[rxnMatrixLoc+j][nchemIdx+i];
        }
        Z_matrix[seqMatrixLoc][i] = sm;
      }
      for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // calculate Z_matrix matrix elements for
        double sm = 0.0;                                                // all sources in the system
        Molecule* pPseudoIsomer = spos->first ;
        int rxnMatrixLoc = spos->second;
        int seqMatrixLoc = m_SpeciesSequence[pPseudoIsomer];
        sm = assymEigenVec[rxnMatrixLoc][nchemIdx+i];
        Z_matrix[seqMatrixLoc][i] = sm;
      }
      if(m_sinkRxns.size()!=0) {
        for(sinkpos=m_sinkRxns.begin(); sinkpos!=m_sinkRxns.end(); ++sinkpos){ // calculate Y_matrix elements for sinks
          double sm = 0.0;
          vector<double> KofEs;                                         // vector to hold sink k(E)s
          vector<double> KofEsTemp;                                     // vector to hold sink k(E)s
          Reaction* sinkReaction = sinkpos->first;
          int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
          if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
            KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
          }
          else{                   // if the collision operator size is >1, there are k(E)s for the irreversible loss
            KofEs = sinkReaction->get_GrainKfmc();                      // assign sink k(E)s, the vector size == maxgrn
            Molecule* isomer = sinkReaction->get_reactant();
            const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();

            // DO NOT MOVE THIS SECTION --- INDEX SENSITIVE
            if (numberGroupedGrains != 0){
              for (int i(numberGroupedGrains - 1); i < colloptrsize; ++i)
                KofEsTemp.push_back(KofEs[i]);
            }
            else{
              for (int i(0); i < colloptrsize; ++i)
                KofEsTemp.push_back(KofEs[i]);
            }
            // DO NOT MOVE THIS SECTION --- INDEX SENSITIVE

            colloptrsize = (!numberGroupedGrains)
              ? colloptrsize
              : colloptrsize - numberGroupedGrains + 1;
          }
          int rxnMatrixLoc = sinkpos->second;                               // get sink location
          int seqMatrixLoc = m_SinkSequence[sinkReaction];                  // get sink sequence position
          for(int j(0);j<colloptrsize;++j){
            sm += assymEigenVec[rxnMatrixLoc+j][nchemIdx+i] * KofEsTemp[j];
          }
          Y_matrix[seqMatrixLoc][i] = sm;
          KofEs.clear();
        }
      }
    }

    //    Y_matrix.print((int)(m_sinkRxns.size()), (int)(m_SpeciesSequence.size())); // print out Y_matrix for testing

    dMatrix Zinv(Z_matrix), Zidentity(nchem), Kr(nchem);
    db2D Kp;

    if(Zinv.invertGaussianJordan()){
      cerr << "Inversion of Z_matrix failed.  Matrix before inversion is: ";
      Z_matrix.showFinalBits(nchem);
    }

    ctest << "\nZ_matrix: ";
    Z_matrix.showFinalBits(nchem, true);

    ctest << endl << "Z_matrix^(-1):" << endl;
    Zinv.showFinalBits(nchem, true);

    for(int i(0);i<nchem;++i){          // multiply Z_matrix*Z_matrix^(-1) for testing
      for(int j(0);j<nchem;++j){
        double sm = 0.0;
        for(int k(0);k<nchem;++k){
          sm += Z_matrix[i][k] * Zinv[k][j];
        }
        Zidentity[i][j] = sm;
      }
    }

    ctest << "\nZ_matrix * Z_matrix^(-1) [Identity matrix]:" << endl;
    Zidentity.showFinalBits(nchem, true);

    for(int i(0);i<nchem;++i){          // calculate Kr (definition taken from PCCP 2007(9), p.4085)
      for(int j(0);j<nchem;++j){
        double sm = 0.0;
        for(int k(0);k<nchem;++k){
          sm += Z_matrix[i][k] * to_double(m_eigenvalues[nchemIdx+k]) * Zinv[k][j];
        }
        Kr[i][j] = sm * m_meanOmega;
      }
    }
    ctest << "\nKr matrix:" << endl;
    Kr.showFinalBits(nchem, true);       // print out Kr_matrix

    if(m_sinkRxns.size()!=0){
      for(int i(0); i != int(m_sinkRxns.size()); ++i){    // calculate Kp (definition taken from PCCP 2007(9), p.4085)
        for(int j(0);j<nchem;++j){
          double sm = 0.0;
          for(int k(0);k<nchem;++k){
            sm += Y_matrix[i][k] * Zinv[k][j];
          }
          Kp[i][j] = sm;
        }
      }
      ctest << "\nKp matrix:" << endl;    // print out Kp_matrix
      Kp.print(int(m_sinkRxns.size()), int(m_SpeciesSequence.size()));
    }

    ctest << "\nFirst order & pseudo first order rate coefficients for loss rxns:\n{\n";
    Reaction::molMapType::iterator lossitr, rctitr, pdtitr;

    stringstream puSymbols;
    stringstream puNumbers;
    // print pseudo 1st order k loss for isomers
    for(lossitr=m_SpeciesSequence.begin(); lossitr!=m_SpeciesSequence.end(); ++lossitr){
      Molecule* iso = lossitr->first;
      int losspos = lossitr->second;
      ctest << iso->getName() << " loss = " << Kr[losspos][losspos] << endl;
      PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderLoss", Kr[losspos][losspos]);
      ppItem->XmlWriteAttribute("ref", iso->getName());

      if (mFlags.searchMethod == 3){
        puNumbers << Kr[losspos][losspos] << "\t";
        if (punchSymbolGathered == false){
          puSymbols << iso->getName() << " loss\t";
        }
      }
    }
    ctest << "}\n";

    if(m_SpeciesSequence.size()>m_sinkRxns.size()){
      ctest << "\nFirst order & pseudo first order rate coefficients for isomerization rxns:\n{\n";

      // print pseudo first order connecting ks
      for (rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
        Molecule* rct = rctitr->first;
        int rctpos = rctitr->second;
        for (pdtitr=m_SpeciesSequence.begin(); pdtitr!=m_SpeciesSequence.end(); ++pdtitr){
          Molecule* pdt = pdtitr->first;
          int pdtpos = pdtitr->second;
          if(rctpos != pdtpos){
            ctest << rct->getName() << " -> " << pdt->getName() << " = " << Kr[pdtpos][rctpos] << endl;

            PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", Kr[pdtpos][rctpos]);
            ppItem->XmlWriteAttribute("fromRef", rct->getName());
            ppItem->XmlWriteAttribute("toRef",   pdt->getName());
            ppItem->XmlWriteAttribute("reactionType", "isomerization");
          }

          if (mFlags.searchMethod == 3){
            puNumbers << Kr[pdtpos][rctpos] << "\t";
            if (punchSymbolGathered == false){
              puSymbols << rct->getName() << " -> " << pdt->getName() << "\t";
            }
          }
        }
      }
      ctest << "}\n";
    }

    if(m_sinkRxns.size()!=0){
      ctest << "\nFirst order & pseudo first order rate coefficients for irreversible rxns:\n{\n";
      sinkMap::iterator sinkitr;

      for(sinkitr=m_SinkSequence.begin(); sinkitr!=m_SinkSequence.end(); ++sinkitr){
        Reaction* sinkReaction = sinkitr->first;          // get Irreversible Rxn
        int colloptrsize = sinkReaction->getRctColloptrsize();
        int sinkpos = m_SinkSequence[sinkReaction];                   // get products & their position
        vector<Molecule*> pdts;
        sinkReaction->get_products(pdts);
        for(rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
          Molecule* rcts = rctitr->first;     // get reactants & their position
          int rctpos = rctitr->second;
          if(colloptrsize==1){
            ctest << rcts->getName() << " -> "  << pdts[0]->getName() << "(bim) = " << Kp[sinkpos][rctpos] << endl;
            if (mFlags.searchMethod == 3){
              puNumbers << Kp[sinkpos][rctpos] << "\t";
              if (punchSymbolGathered == false){
                puSymbols << rcts->getName() << " -> " << pdts[0]->getName() << "(bim)\t";
              }
            }
          }
          else{
            ctest << rcts->getName() << " -> "  << pdts[0]->getName() << " = " << Kp[sinkpos][rctpos] << endl;

            PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", Kp[sinkpos][rctpos]);
            ppItem->XmlWriteAttribute("fromRef", rcts->getName());
            ppItem->XmlWriteAttribute("toRef",   pdts[0]->getName());
            ppItem->XmlWriteAttribute("reactionType", "irreversible");
            if (mFlags.searchMethod == 3){
              puNumbers << Kp[sinkpos][rctpos] << "\t";
              if (punchSymbolGathered == false){
                puSymbols << rcts->getName() << " -> " << pdts[0]->getName() << "\t";
              }
            }
          }
        }
      }
      ctest << "}\n\n";
    }

    if (puSymbols.str().size()) {
      puSymbols << "\n";
      mFlags.punchSymbols = puSymbols.str();
      punchSymbolGathered = true;
    }

    if (puNumbers.str().size()) {
      puNumbers << "\n";
      mFlags.punchNumbers = puNumbers.str();
    }

    mesmerRates = Kr;
    return true;
  }

  // Set Initial population for individual species
  void ReactionManager::setInitialPopulation(PersistPtr anchor)
  {
    PersistPtr pp=anchor;
    double populationSum = 0.0;
    PersistPtr ppInitMol = pp->XmlMoveTo("molecule");
    while(ppInitMol){
      string sRef = ppInitMol->XmlReadValue("ref");
      if(sRef.size()){ // if got the name of the molecule
        Molecule* pMolecule = m_pMoleculeManager->find(sRef) ;
        double population = ppInitMol->XmlReadDouble("me:population", optional) ;
        if (population > 0.0){
          populationSum += population;
          pMolecule->getPop().setInitPopulation(population);
          ctest << "Initial population of " << pMolecule->getName() << " = " << population << endl;
        }
      }
      ppInitMol = ppInitMol->XmlMoveTo("molecule");
    }
    if (populationSum != 1.0){
      if (populationSum > 0.0){
        // Populations need to be Normalized.
      } else if (populationSum == 0.0){
        // Issue warning that there are no populations set, and calculate only rate coefficients.
        populationSum += 1.0;
        //      pPseudoIsomer->getPop().setInitPopulation(populationSum);
      } else {
        // Issue error and stop.
      }
    }
  }

  double ReactionManager::calcChiSquare(const dMatrix& mesmerRates, vector<conditionSet>& expRates){
    double chiSquare = 0.0;

    for (size_t i(0); i < expRates.size(); ++i){
      string ref1, ref2; double expRate(0.0), expErr(0.0); int seqMatrixLoc1(-1), seqMatrixLoc2(-1);
      expRates[i].get_conditionSet(ref1, ref2, expRate, expErr);

      // check and get the position of ref1 and ref2 inside m_SpeciesSequence vector
      Reaction::molMapType::iterator spcitr;
      for (spcitr = m_SpeciesSequence.begin(); spcitr != m_SpeciesSequence.end(); ++spcitr){
        if (ref1 == (spcitr->first)->getName()) {
          seqMatrixLoc1 = spcitr->second;
          break;
        }
      }
      if (seqMatrixLoc1 == -1){
        cerr << "No molecule named " << ref1 << " is available in the reaction species.";
        break;
      }
      for (spcitr = m_SpeciesSequence.begin(); spcitr != m_SpeciesSequence.end(); ++spcitr){
        if (ref2 == (spcitr->first)->getName()) {
          seqMatrixLoc2 = spcitr->second;
          break;
        }
      }
      if (seqMatrixLoc2 == -1){
        cerr << "No molecule named " << ref2 << " is available in the reaction species.";
        break;
      }

      double diff = (mesmerRates[seqMatrixLoc1][seqMatrixLoc2] - expRate);
      chiSquare += (diff * diff) / (expErr * expErr);
    }


    return chiSquare;
  }

  // This method constructs a transition matrix based on energy grains.
  //
  void ReactionManager::constructGrainMatrix(int msize){

    // Allocate space for the full system collision operator.
    if (m_reactionOperator) delete m_reactionOperator;
    m_reactionOperator = new qdMatrix(msize, 0.0) ;

    // Insert collision operators to reaction operator from individual wells.
    Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
    for (isomeritr = m_isomers.begin() ; isomeritr != m_isomers.end() ; ++isomeritr) {

      Molecule *isomer = isomeritr->first ;
      int colloptrsize = isomer->getColl().getNumberOfGroupedGrains() != 0
        ? isomer->getColl().get_colloptrsize() - isomer->getColl().getNumberOfGroupedGrains() + 1
        : isomer->getColl().get_colloptrsize();
      double omega = isomer->getColl().get_collisionFrequency();
      int idx = isomeritr->second ;

      isomer->getColl().copyCollisionOperator(m_reactionOperator, colloptrsize, idx, omega/m_meanOmega) ;

    }

    // Add connecting rate coefficients.
    for (size_t i(0) ; i < size() ; ++i) {
      m_reactions[i]->AddReactionTerms(m_reactionOperator,m_isomers,1.0/m_meanOmega) ;
    }

  }

  // This is a routine to construct the big basis matrix based on the alternative basis set method.
  // The full reaction operator is subject to a similarity transformation process by a set of eigenvectors.
  // If there are three wells and two sources in the system, and the eigenvectors of each well under the assumption
  // of the conservation of the wells are U_0, U_1 and U_2, respectively. The transformer matrix should look like
  //
  //        [  U_0   0    0   0   0 ]
  //        [   0   U_1   0   0   0 ]
  //    U = [   0    0   U_2  0   0 ]
  //        [   0    0    0   1   0 ]
  //        [   0    0    0   0   1 ]
  //
  // This transformer matrix operates on the reaction operator to give the basis matrix by doing
  //
  //     M'' = U^-1 M U
  //
  // One then needs to decide how many members of this basis matrix to include in the reduced basis matrix for
  // diagonalization.
  //
  void ReactionManager::constructBasisMatrix(void){

    // Determine the size and location of various blocks.

    // 1. Isomers.

    size_t msize(0) ;
    Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
    for (; isomeritr != m_isomers.end() ; ++isomeritr) {
      Molecule *isomer = isomeritr->first ;
      isomeritr->second = static_cast<int>(msize) ; //set location
      msize += isomer->getColl().get_nbasis() ;
    }

    // 2. Pseudoisomers.

    Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin() ;
    for (; pseudoIsomeritr != m_isomers.end() ; ++pseudoIsomeritr) {
      pseudoIsomeritr->second = static_cast<int>(msize) ; //set location
      msize++ ;
    }

    // Allocate space for the reaction operator.

    if (m_reactionOperator) delete m_reactionOperator;
    m_reactionOperator = new qdMatrix(msize, 0.0) ;

    // Insert collision operators: in the contracted basis these are the eignvalues
    // of the isomer collision operators.
    for (isomeritr = m_isomers.begin() ; isomeritr != m_isomers.end() ; ++isomeritr) {

      Molecule *isomer = isomeritr->first ;
      double omega = isomer->getColl().get_collisionFrequency() ;
      int idx = isomeritr->second ;

      isomer->getColl().copyCollisionOperatorEigenValues(m_reactionOperator, idx, omega) ;
    }

    // Add rate coefficients.
    for (size_t i(0) ; i < size() ; ++i) {
      m_reactions[i]->AddContractedBasisReactionTerms(m_reactionOperator,m_isomers) ;
    }

    // Print out system matrix.

    //ctest << endl << "System matrix:" << endl << endl ;
    //for (size_t i(0) ; i < msize ; ++i) {
    //  for (size_t j(0) ; j < msize ; ++j) {
    //    formatFloat(ctest, (*m_reactionOperator)[i][j],  6,  15) ; 
    //  }
    //  ctest << endl ;
    //}

  }

  //
  // The vector produced by this function contains the sqrt of
  // the normalized equilibrium distribution.
  //
  bool ReactionManager::produceReducedEquilibriumVector()
  {

    m_reducedEqVector.clear();
    m_reducedEqVector.resize(m_reducedEigenvectors->size());


    for (int i(0); i < int(m_locSizeMap.size()); ++i){  // iterate through source map to get
      Molecule* pMol = m_locSizeMap[i].mol;                            // eq Fractions
      int rxnMatrixLoc = m_locSizeMap[i].rml;
      double eqFrac = pMol->getPop().getEqFraction();
      m_reducedEqVector[rxnMatrixLoc] = sqrt(eqFrac);
    }

    return true;
  }

}//namespace


