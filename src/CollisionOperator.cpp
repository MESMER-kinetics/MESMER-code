//-------------------------------------------------------------------------------------------
//
// ReactionManager.cpp
//
// Author: Struan Robertson
// Date:   26/Feb/2011
//
// This file contains implementation of the master equation collision operator class.
//
//-------------------------------------------------------------------------------------------
#include "CollisionOperator.h"

#include "AssociationReaction.h"
#include "IrreversibleUnimolecularReaction.h"
#include "IsomerizationReaction.h"
#include "IrreversibleExchangeReaction.h"

namespace mesmer
{

  CollisionOperator::CollisionOperator() : m_pMoleculeManager(0), 
    m_pReactionManager(0), 
    m_isomers(),
    m_sources(),
    m_sinkRxns(),
    m_meanOmega(0.0),
    m_reactionOperator(0),
    m_eigenvectors(0),
    m_eigenvalues(),
    m_SpeciesSequence(),
    m_eqVector(),
    m_punchSymbolGathered(false) {}

  CollisionOperator::~CollisionOperator() {
    if (m_reactionOperator) delete m_reactionOperator;
  }

  // Initialize the collision operator object.
  bool CollisionOperator::initialize(MoleculeManager *pMoleculeManager, ReactionManager *pReactionManager) {

    if ( (m_pMoleculeManager = pMoleculeManager) && 
      (m_pReactionManager = pReactionManager)) return true ;

    return false ;
  };

  //
  // Main methods for constructing the Collision operator.
  //
  bool CollisionOperator::BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags, bool writeReport)
  {
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
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
      double TS_ZPE(INFIMUM);

      Reaction *pReaction = (*m_pReactionManager)[i] ;

      // Reset the the microcanonical re-calculation flags if required.
      if (!mFlags.useTheSameCellNumber) pReaction->resetCalcFlag();

      // Transition State
      // third check for the transition state in this reaction
      Molecule *pTransitionState = pReaction->get_TransitionState();
      if (pTransitionState){
        TS_ZPE = pTransitionState->getDOS().get_zpe();
        maxEnergy = max(maxEnergy, TS_ZPE) ;
      }

      // unimolecular species
      vector<Molecule *> unimolecules ;
      pReaction->get_unimolecularspecies(unimolecules) ;
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
      AssociationReaction *pAReaction = dynamic_cast<AssociationReaction*>(pReaction) ;
      if (pAReaction) {
        double pseudoIsomerZPE = pAReaction->get_pseudoIsomer()->getDOS().get_zpe();
        double excessReactantZPE = pAReaction->get_excessReactant()->getDOS().get_zpe();
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
      // For irreversible exchange reactions determine zero point energy location of the
      // associating pair.
      //
      IrreversibleExchangeReaction *pIEReaction = dynamic_cast<IrreversibleExchangeReaction*>(pReaction) ;
      if (pIEReaction) {
        double pseudoIsomerZPE = pIEReaction->get_pseudoIsomer()->getDOS().get_zpe();
        double excessReactantZPE = pIEReaction->get_excessReactant()->getDOS().get_zpe();
        double sourceTermZPE = pseudoIsomerZPE + excessReactantZPE;
        minEnergy = min(minEnergy, sourceTermZPE) ;
        maxEnergy = max(maxEnergy, sourceTermZPE) ;

        // There is no well for this reaction
      }

      //
      // For dissociation reactions determine zero point energy location of the barrier
      //
      IrreversibleUnimolecularReaction *pDissnRtn = dynamic_cast<IrreversibleUnimolecularReaction*>(pReaction) ;
      if (pDissnRtn) {
        const double rctZPE = pDissnRtn->get_reactant()->getDOS().get_zpe();
        double barrierZPE = rctZPE + pDissnRtn->get_ThresholdEnergy();
        minEnergy = min(minEnergy, barrierZPE) ;
        maxEnergy = max(maxEnergy, barrierZPE) ;

        // Calculate the lowest barrier associated with this well(species).
        if (barrierZPE < unimolecules[0]->getColl().getLowestBarrier()){
          unimolecules[0]->getColl().setLowestBarrier(barrierZPE);
        }
      }

    }

    // set grain parameters for the current Temperature/pressure condition
    if(!SetGrainParams(mEnv, mFlags, minEnergy, maxEnergy, writeReport))
      return false;

    // Calculate flux and k(E)s
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
      if(!(*m_pReactionManager)[i]->calcGrnAvrgMicroRateCoeffs())
        return false;
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
          if (isomer->getColl().isCemetery()) msize -= 1;
        }

        m_meanOmega += isomer->getColl().get_collisionFrequency() ;
      }
      m_meanOmega /= double(m_isomers.size());

      //
      // Find all source terms. Note: a source term contains the deficient reactant.
      // It is possible for there to be more than one source term.
      //
      m_sources.clear(); // Maps the location of source in the system matrix.
      for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
        (*m_pReactionManager)[i]->updateSourceMap(m_sources) ;
      }

      // Build reaction operator.
      //
      // One of two methods for building the reaction operator are available:
      // the conventional energy grained master equation method which is based
      // on energy grains and a contracted basis set method in which a basis
      // set is generated from the individual collision operators and a
      // representation of the reaction operator build upon this basis.

      if (!mEnv.useBasisSetMethod) {

        // Full energy grained reaction operator.

        constructGrainMatrix(msize);

      } else {

        // Contracted basis set reaction operator.

        constructBasisMatrix();

      }
    }

    return true;
  }

  // Sets grain parameters and determine system environment.
  bool CollisionOperator::SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne, bool writeReport)
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
    const double thermalEnergy = (mFlags.useTheSameCellNumber) ? MaximumTemperature * boltzmann_RCpK : 1.0/mEnv.beta ;
    mEnv.EMax += mEnv.EAboveHill * thermalEnergy;

    if (mEnv.GrainSize <= 0.0){
      mEnv.GrainSize = 100; //default 100cm-1
      cerr << "Grain size was invalid. Reset grain size to default: 100";
    }

    mEnv.MaxGrn = int((mEnv.EMax-mEnv.EMin)/mEnv.GrainSize + 0.5);
    mEnv.MaxCell = mEnv.GrainSize * mEnv.MaxGrn;

    if (writeReport) cinfo << "Cell number = " << mEnv.MaxCell << ", Grain number = " << mEnv.MaxGrn << endl;

    return true;
  }

  // This method constructs a transition matrix based on energy grains.
  //
  void CollisionOperator::constructGrainMatrix(int msize){

    // Determine the size and location of various blocks.

    // 1. Isomers.

    //size_t msize(0) ;
    //Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
    //for (; isomeritr != m_isomers.end() ; ++isomeritr) {
    //  Molecule *isomer = isomeritr->first ;
    //  isomeritr->second = static_cast<int>(msize) ; //set location
    //  msize += isomer->getColl().get_nbasis() ;
    //}

    // 2. Pseudoisomers.

    Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin() ;
    for (; pseudoIsomeritr != m_sources.end() ; ++pseudoIsomeritr) {
      pseudoIsomeritr->second = static_cast<int>(msize) ; //set location
      msize++ ;
    }

    // Allocate space for the full system collision operator.
    if (m_reactionOperator) delete m_reactionOperator;
    m_reactionOperator = new qdMatrix(msize, 0.0) ;

    // Insert collision operators to reaction operator from individual wells.
    Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
    for (isomeritr = m_isomers.begin() ; isomeritr != m_isomers.end() ; ++isomeritr) {

      Molecule *isomer = isomeritr->first ;
      int colloptrsize = isomer->getColl().getNumberOfGroupedGrains() != 0
        ? ( isomer->getColl().isCemetery()
        ? isomer->getColl().get_colloptrsize() - isomer->getColl().getNumberOfGroupedGrains()
        : isomer->getColl().get_colloptrsize() - isomer->getColl().getNumberOfGroupedGrains() + 1)
        : isomer->getColl().get_colloptrsize();
      double omega = isomer->getColl().get_collisionFrequency();
      int idx = isomeritr->second ;

      isomer->getColl().copyCollisionOperator(m_reactionOperator, colloptrsize, idx, omega/m_meanOmega) ;

    }

    // Add connecting rate coefficients.
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
      (*m_pReactionManager)[i]->AddReactionTerms(m_reactionOperator,m_isomers,1.0/m_meanOmega) ;
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
  void CollisionOperator::constructBasisMatrix(void){

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
    for (; pseudoIsomeritr != m_sources.end() ; ++pseudoIsomeritr) {
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
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
      (*m_pReactionManager)[i]->AddContractedBasisReactionTerms(m_reactionOperator,m_isomers) ;
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

  bool CollisionOperator::calculateEquilibriumFractions()
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
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {  //iterate through m_reactions

      Molecule* rct;
      Molecule* pdt;
      double Keq(0.0);

      //only need eq fracs for species in isom & assoc rxns
      if ((*m_pReactionManager)[i]->isEquilibratingReaction(Keq, &rct, &pdt)){

        int ploc(0), rloc(0) ;

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
      if (m_isomers.size()){
        Molecule* rct=(m_isomers.begin())->first;
        m_SpeciesSequence[rct] = counter;
      }
      else if (m_sources.size()){
        Molecule* rct=(m_sources.begin())->first;
        m_SpeciesSequence[rct] = counter;
      }
      else{
        return false;
      }
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
      string speciesName = key->getName();
      ctest << "Equilibrium Fraction for " << speciesName << " = " << key->getPop().getEqFraction() << endl;
    }

    // Calculate equilibrium vector.
    if (!produceEquilibriumVector()){
      cerr << "Calculation of equilibrium vector failed.";
      return false;
    }

    return true;
  }

  void CollisionOperator::diagReactionOperator(const MesmerFlags &mFlags, const MesmerEnv &mEnv,
    Precision precision, PersistPtr ppAnalysis)
  {
    // Allocate space for eigenvalues.
    const size_t smsize = m_reactionOperator->size() ;
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
    case DOUBLE: 
      {
        dMatrix dDiagM(smsize);
        for ( size_t i = 0 ; i < smsize ; ++i )
          for ( size_t j = 0 ; j < smsize ; ++j )
            dDiagM[i][j] = to_double((*m_reactionOperator)[i][j]) ;
        vector<double>  dEigenValue(smsize, 0.0);
        dDiagM.diagonalize(&dEigenValue[0]) ;
        for ( size_t i = 0 ; i < smsize ; ++i )
          m_eigenvalues[i] = dEigenValue[i];
        for ( size_t i = 0 ; i < smsize ; ++i )
          for ( size_t j = 0 ; j < smsize ; ++j )
            (*m_eigenvectors)[i][j] = dDiagM[i][j] ;
        break;
      }
    case DOUBLE_DOUBLE: 
      {
        ddMatrix ddDiagM(smsize);
        for ( size_t i = 0 ; i < smsize ; ++i )
          for ( size_t j = 0 ; j < smsize ; ++j )
            ddDiagM[i][j] = to_dd_real((*m_reactionOperator)[i][j]) ;
        vector<dd_real> ddEigenValue(smsize, 0.0);
        ddDiagM.diagonalize(&ddEigenValue[0]) ;
        for ( size_t i = 0 ; i < smsize ; ++i )
          m_eigenvalues[i] = ddEigenValue[i];
        for ( size_t i = 0 ; i < smsize ; ++i )
          for ( size_t j = 0 ; j < smsize ; ++j )
            (*m_eigenvectors)[i][j] = ddDiagM[i][j] ;
        break;
      }
    default: // diagonalize in quad double
      {
        (*m_eigenvectors) = (*m_reactionOperator) ;
        m_eigenvectors->diagonalize(&m_eigenvalues[0]) ;
      }

    }
    // diagonalize the whole matrix
    //-------------------------------------------------------------

    // This block prints Eigenvectors
    if (mFlags.printReactionOperatorNum){
      ctest << "Eigenvectors --- ";
      stringstream os;
      printEigenvectors(mFlags, os);
      ctest << os.str();
    }

    if(mFlags.printEigenValuesNum!=0)
    {
      size_t numberStarted = 0; //will apply when mFlags.printEigenValuesNum<0: print all
      if (mFlags.printEigenValuesNum > 0 && mFlags.printEigenValuesNum <= int(smsize))
        numberStarted = smsize - mFlags.printEigenValuesNum;

      PersistPtr ppEigenList = ppAnalysis->XmlWriteElement("me:eigenvalueList");
      ppEigenList->XmlWriteAttribute("number",toString(smsize));
      ppEigenList->XmlWriteAttribute("selection",toString(mFlags.printEigenValuesNum));//TODO improve this
      ctest << "\nTotal number of eigenvalues = " << smsize << endl;
      ctest << "Eigenvalues\n{\n";
      for (size_t i = numberStarted ; i < smsize; ++i) {
        qd_real tmp = (mEnv.useBasisSetMethod)? m_eigenvalues[i] : m_eigenvalues[i] * m_meanOmega ;
        formatFloat(ctest, tmp, 6, 15) ;
        ctest << endl ;
        ppEigenList->XmlWriteValueElement("me:eigenvalue", to_double(tmp), 6);
      }
      ctest << "}\n";
    }

  }

  void CollisionOperator::printReactionOperator(const MesmerFlags& mFlags)
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

  void CollisionOperator::printEigenvectors(const MesmerFlags& mFlags, std::ostream& os)
  {
    const size_t smsize = m_eigenvectors->size() ;

    switch (mFlags.printReactionOperatorNum)
    {
    case -1:
      os << "Printing all (" << smsize << ") columns/rows of the eigenvectors:\n";
      (*m_eigenvectors).showFinalBits(smsize, mFlags.print_TabbedMatrices);
      break;
    case -2:
      os << "Printing final 1/2 (" << smsize/2 << ") columns/rows of the eigenvectors:\n";
      (*m_eigenvectors).showFinalBits(smsize/2, mFlags.print_TabbedMatrices);
      break;
    case -3:
      os << "Printing final 1/3 (" << smsize/3 << ") columns/rows of the eigenvectors:\n";
      (*m_eigenvectors).showFinalBits(smsize/3, mFlags.print_TabbedMatrices);
      break;
    default: // the number is either smaller than -3 or positive
      if (abs(mFlags.printReactionOperatorNum) > int(smsize)){
        os << "Printing all (" << smsize << ") columns/rows of the eigenvectors:\n";
        (*m_eigenvectors).showFinalBits(smsize, mFlags.print_TabbedMatrices);
      }
      else{
        os << "Printing final " << abs(mFlags.printReactionOperatorNum) << " columns/rows of the eigenvectors:\n";
        (*m_eigenvectors).showFinalBits(abs(mFlags.printReactionOperatorNum), mFlags.print_TabbedMatrices);
      }
    }
  }  

  bool CollisionOperator::timeEvolution(MesmerFlags& mFlags, PersistPtr ppAnalysis, PersistPtr ppPopList)
  {
    ErrorContext c(__FUNCTION__);

    // Cut short if species profiles not needed.
    if(!mFlags.speciesProfileEnabled)
      return true;

    size_t smsize = m_eigenvectors->size();
    vector<double> r_0(smsize, 0.);
    if (!projectedInitialDistrbtn(r_0)) {
      cerr << "Projection of initial disttribution failed.";
      return false;
    }

    double shortestTime = 0.;
    // set the default maximum evolution time
    if (mFlags.shortestTimeOfInterest < 1.0e-20 || mFlags.shortestTimeOfInterest > 1.0)
      shortestTime = 1.0e-11;
    else
      shortestTime = mFlags.shortestTimeOfInterest;

    double maxEvoTime = 0.;
    // set the default maximum evolution time
    if (mFlags.maxEvolutionTime <= 0.001 || mFlags.maxEvolutionTime > 1.0e8)
      maxEvoTime = 1.2e5;
    else
      maxEvoTime = mFlags.maxEvolutionTime;

    // Calculates the time points
    vector<double> timePoints;
    for (int i = 0; i <= 300; ++i){
      double thetime = pow(10., static_cast<double>(i) / 10. - 20.);
      if (thetime < shortestTime) continue;
      if (thetime > maxEvoTime) break;
      timePoints.push_back(thetime);
    }

    //Initialises dt vector for calculating product yields
    vector<double> dt(timePoints.size()-1,0.0);
    dt[0] = timePoints[0];
    for (int i = 1; i < int(dt.size()); ++i){
      dt[i] = timePoints[i] - timePoints[i-1];
    }

    dMatrix totalEigenVecs(smsize); // copy full eigenvectors of the system
    for (size_t i = 0; i < smsize; ++i) {
      double tmp = to_double(m_eqVector[i]);
      for (size_t j = 0; j < smsize; ++j) {
        totalEigenVecs[i][j] = tmp*to_double((*m_eigenvectors)[i][j]);
      }
    }

    const size_t maxTimeStep = dt.size();
    db2D grnProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
    vector<double> p_t(smsize, 0.);

    for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
      double numColl = m_meanOmega * timePoints[timestep];
      for (size_t j(0); j < smsize; ++j) {
        p_t[j] = r_0[j] * exp(to_double(m_eigenvalues[j]) * numColl);
      } // now |p_t> = exp(Lambda*t)*V^(T)*|init> = exp(Lambda*t)*U^(-1)*|n_0>
      
      p_t *= totalEigenVecs ;

      for (size_t j(0); j < smsize; ++j) {
        grnProfile[j][timestep] = p_t[j];
      } // now |grnProfile(t)> = |grnProfile(i)> = F*V*exp(Lambda*t)*V^(T)*|init> = U*exp(Lambda*t)*U^(-1)*|n_0>

    }

    //------------------------------
    // print grained species profile
    if (mFlags.grainedProfileEnabled) {
      ctest << "\nGrained species profile (the first row is time points in unit of second):\n{\n";
      for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
        formatFloat(ctest, timePoints[timestep], 6,  15);
      }
      ctest << endl;
      for (size_t j(0); j < smsize; ++j) {
        for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
          formatFloat(ctest, grnProfile[j][timestep], 6,  15);
        }
        ctest << endl;
      }
      ctest << "}\n";

      PersistPtr ppGrainList = ppAnalysis->XmlWriteElement("me:grainPopulationList");
      size_t timestep = maxTimeStep/2; //temporary value
      { 
        PersistPtr ppGrainPop = ppGrainList->XmlWriteElement("me:grainPopulation");
        ppGrainPop->XmlWriteAttribute("time", toString(timePoints[timestep]));
        ppGrainPop->XmlWriteAttribute("logTime", toString(log10(timePoints[timestep])));
        for(size_t j = 0; j < smsize; j+=5)  
        {
          PersistPtr ppGrain = ppGrainPop->XmlWriteValueElement("me:grain", to_double(grnProfile[j][timestep]), 6);
          ppGrain->XmlWriteAttribute("index", toString(j));
        }
      }
    }

    ctest<<"mean collision frequency = " << m_meanOmega << "/s" << endl;

    vector<double> totalIsomerPop(maxTimeStep, 0.);
    vector<double> totalPdtPop(maxTimeStep, 0.);

    for(size_t timestep(0); timestep<maxTimeStep; ++timestep){
      for(size_t j(0);j<smsize;++j){
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

    if (mFlags.speciesProfileEnabled){
      ctest << endl << "Print time dependent species and product profiles" << endl << "{" << endl;
      int numberOfSpecies = static_cast<int>(m_isomers.size() + m_sources.size() + m_sinkRxns.size());

      //---------------------------------------------------------------------------------------------
      // Need to include the cemetery states too, so loop into isomers and see how many have cemetery.
      Reaction::molMapType::iterator iposC;
      for (iposC = m_isomers.begin(); iposC != m_isomers.end(); ++iposC){  // Iterate through the isomer map
        Molecule* isomer = iposC->first;
        if (isomer->getColl().isCemetery()) ++numberOfSpecies;
      }
      //---------------------------------------------------------------------------------------------

      db2D speciesProfile(numberOfSpecies, maxTimeStep);
      int speciesProfileidx(0);

      ctest << setw(16) << "Timestep/s";

      // Iterate through the source map, to determine the total source 
      // density as a function of time.
      vector<string> speciesNames;
      Reaction::molMapType::iterator spos;
      for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  
        Molecule* source = spos->first ;
        ctest << setw(16) << source->getName();
        speciesNames.push_back(source->getName());
        int rxnMatrixLoc = spos->second;
        for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
          speciesProfile[speciesProfileidx][timestep] = grnProfile[rxnMatrixLoc][timestep];
        }
        ++speciesProfileidx;
      }

      // Iterate through the isomer map, to calculate the total isomer 
      // density as a function of time.
      Reaction::molMapType::iterator ipos;
      for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  
        Molecule* isomer = ipos->first;                        
        string isomerName = isomer->getName();
        if (isomer->getColl().isCemetery()){
          isomerName += "(+)"; // active states
        }
        ctest << setw(16) << isomerName;
        speciesNames.push_back(isomerName);
        int rxnMatrixLoc = ipos->second;
        int colloptrsize = isomer->getColl().get_colloptrsize();
        const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
        if (numberGrouped != 0) { 
          const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
          colloptrsize -= (numberGrouped - nrg) ;
        }
        for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
          for(int i = 0; i < colloptrsize; ++i){
            speciesProfile[speciesProfileidx][timestep] += grnProfile[i+rxnMatrixLoc][timestep];
          }
        }
        ++speciesProfileidx;
      }

      int pdtProfileStartIdx = speciesProfileidx;

      // Taking account of the cemetery states in all wells.
      for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
        Molecule* isomer = ipos->first;
        if (isomer->getColl().isCemetery()){
          vector<double> grainKdmc = isomer->getColl().get_GrainKdmc();
          string cemName = isomer->getName() + "(-)";
          ctest << setw(16) << cemName;
          speciesNames.push_back(cemName);
          int rxnMatrixLoc = ipos->second;                       // get isomer location
          double TimeIntegratedCemeteryPop(isomer->getPop().getInitCemeteryPopulation());
          for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
            for(size_t i = 0; i < grainKdmc.size(); ++i){
              speciesProfile[speciesProfileidx][timestep] += m_meanOmega * grainKdmc[i]*grnProfile[i+rxnMatrixLoc][timestep]*dt[timestep];;
            }
            TimeIntegratedCemeteryPop += speciesProfile[speciesProfileidx][timestep];
            speciesProfile[speciesProfileidx][timestep]= TimeIntegratedCemeteryPop;
          }
          ++speciesProfileidx;
        }
      }

      sinkMap::iterator pos;      // iterate through sink map to get product profile vs t
      for (pos = m_sinkRxns.begin(); pos != m_sinkRxns.end(); ++pos){
        vector<double> KofEs;                             // vector to hold sink k(E)s
        Reaction* sinkReaction = pos->first;
        const int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
        vector<Molecule*> pdts;                               // in the sink reaction
        sinkReaction->get_products(pdts);

        size_t numberGrouped(0);
        string pdtName = pdts[0]->getName();
        if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
          KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
          pdtName += "(bim)";
          ctest << setw(16) << pdtName;
        }
        else{   // if the collision operator size is >1, there are k(E)s for the irreversible loss
          KofEs = sinkReaction->get_GrainKfmc();          // assign sink k(E)s, the vector size == maxgrn
          ctest << setw(16) << pdtName;
          numberGrouped = sinkReaction->get_reactant()->getColl().getNumberOfGroupedGrains();
        }
        speciesNames.push_back(pdtName);
        int rxnMatrixLoc = pos->second;                       // get sink location
        double TimeIntegratedProductPop(0.0);

        size_t idx(0) ;
        if (numberGrouped != 0){
          Molecule* rctMol = pos->first->get_reactant();
          const int nrg = rctMol->getColl().isCemetery() ? 0 : 1;
          idx = numberGrouped - nrg ;
        }
        for (size_t timestep(0); timestep < maxTimeStep; ++timestep){
          for(size_t i(0); i < colloptrsize - idx; ++i){
            speciesProfile[speciesProfileidx][timestep] += KofEs[i + idx]*grnProfile[i+rxnMatrixLoc][timestep]*dt[timestep];
          }
          TimeIntegratedProductPop += speciesProfile[speciesProfileidx][timestep];
          speciesProfile[speciesProfileidx][timestep] = TimeIntegratedProductPop;
        }
        ++speciesProfileidx ;
      }

      if (pdtProfileStartIdx < speciesProfileidx){
        for(size_t timestep(0); timestep < maxTimeStep; ++timestep){    // normalize product profile to account for small
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

      //Write to ctest and XML
      ctest << setw(16)<< "totalIsomerPop" << setw(16)<< "totalPdtPop"  << endl;
      for(size_t timestep(0); timestep < maxTimeStep; ++timestep){
        ctest << setw(16) << timePoints[timestep];
        PersistPtr ppPop =  ppPopList->XmlWriteElement("me:population");
        ppPop->XmlWriteAttribute("time", toString(timePoints[timestep]));
        ppPop->XmlWriteAttribute("logTime", toString(log10(timePoints[timestep])));
        for(int i(0); i<speciesProfileidx; ++i){
          ctest << setw(16) << speciesProfile[i][timestep];
          PersistPtr ppVal = ppPop->XmlWriteValueElement("me:pop", speciesProfile[i][timestep]);
          ppVal->XmlWriteAttribute("ref", speciesNames[i]);
        }
        ctest << setw(16) << totalIsomerPop[timestep] << setw(16) << totalPdtPop[timestep] << endl;
      }
      ctest << "}" << endl;
    }
    return true;
  }

  bool CollisionOperator::produceEquilibriumVector()
  {

    m_eqVector.clear();
    m_eqVector.resize(m_reactionOperator->size());

    Reaction::molMapType::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // Iterate through the source map to get
      Molecule* source = spos->first;                                 // the equilibrum fractions.
      int rxnMatrixLoc = spos->second;
      qd_real eqFrac = source->getPop().getEqFraction();
      m_eqVector[rxnMatrixLoc] = sqrt(eqFrac) ;
    }

    Reaction::molMapType::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // Iterate through the isomer map
      Molecule* isomer = ipos->first;                                 // to get the equilibrium fractions.
      int rxnMatrixLoc = ipos->second;
      qd_real eqFrac = isomer->getPop().getEqFraction();
      const int colloptrsize = isomer->getColl().get_colloptrsize();
      const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
      const int nrg = isomer->getColl().isCemetery() ? 0 : 1 ;
      vector<double> boltzFrac;
      isomer->getColl().normalizedGrnBoltzmannDistribution(boltzFrac, colloptrsize, numberGrouped);
      if (numberGrouped == 0) {
        for(int i(0);i<colloptrsize;++i){
          m_eqVector[rxnMatrixLoc + i] = sqrt(eqFrac * qd_real(boltzFrac[i]) ) ;
        }
      }
      else{
        for(int i(1-nrg), j(0);i<colloptrsize - numberGrouped + 1;++i, ++j){
          m_eqVector[rxnMatrixLoc + j] = sqrt(eqFrac * qd_real(boltzFrac[i]) ) ;
        }
      }
    }
    return true;
  }

  bool CollisionOperator::produceInitialPopulationVector(vector<double>& n_0) const {

    double populationSum = 0.0;

    Reaction::molMapType::const_iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
      Molecule* isomer = ipos->first;                        // to get isomer initial populations
      populationSum += isomer->getPop().getInitPopulation();
    }

    Reaction::molMapType::const_iterator spos;
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
        const int nrg = isomer->getColl().isCemetery() ? 0 : 1 ;
        if (numberGrouped == 0){
          for (int i = 0; i < colloptrsize; ++i){
            n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
          }
        }
        else{
          if (!nrg){
            isomer->getPop().setInitCemeteryPopulation(initFrac * boltzFrac[0]);
          }
          for (int i(1-nrg), j(0); i < colloptrsize - numberGrouped + 1; ++i, ++j){
            n_0[j + rxnMatrixLoc] = initFrac * boltzFrac[i];
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
      const int nrg = isomer->getColl().isCemetery() ? 0 : 1 ;
      if (numberGrouped == 0){
        for (int i = 0; i < colloptrsize; ++i){
          n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
        }
      }
      else{
        for (int i = 0; i < colloptrsize - numberGrouped + nrg; ++i){
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

  bool CollisionOperator::BartisWidomPhenomenologicalRates(qdMatrix& mesmerRates, MesmerFlags& mFlags, PersistPtr ppList)
  {
    // Constants.
    const size_t smsize   = m_eigenvectors->size() ;
    const size_t nchem    = m_isomers.size() + m_sources.size() ;  // number of isomers+pseudoisomers
    const size_t nchemIdx = smsize - nchem ;       // Location of chemically significant eigenvalues & vectors
    const size_t nsinks   = m_sinkRxns.size() ;    // Number of Sinks.

    ctest << "\nBartis Widom eigenvalue/eigenvector analysis\n" << endl ;
    ctest << "Number of sinks in this system: " << nsinks << endl;

    if(nsinks > 0){
      ctest << "\nThere should be " << nchem << " chemically significant eigenvalues (CSEs)" << endl;
    } else {
      ctest << "\nThere should be 1 zero eigenvalue (zero within numerical precision) and " << nchem-1
        << " chemically significant eigenvalues (CSEs)" << endl;
    }

    //
    // If there are no sinks, replace the equilibrium vector with the eigenvector whose
    // associated eigenvalue is zero, as this is a consistent estimate of the equilibrium 
    // with respect to the other eigenvalues. Also, as the system is conservative, set the 
    // smallest eigenvalue explicitly to zero.
    //
    if (nsinks < 1) {
      m_eigenvalues[smsize-1] = 0.0 ;
      for(size_t i(0) ; i<smsize ; ++i){
        m_eqVector[i] = (*m_eigenvectors)[i][smsize-1];
      }
    }

    //
    // Construct assymmetric eigenvectors required for the z matrix.
    //
    qdMatrix assymInvEigenVec(smsize);   // U^(-1)
    qdMatrix assymEigenVec(smsize);      // U
    for(size_t i(0) ; i<smsize ; ++i){
      qd_real tmp = m_eqVector[i];
      qd_real sm(0) ;
      for(size_t j(0) ; j<smsize ; ++j){
        assymInvEigenVec[j][i] = (*m_eigenvectors)[i][j]/tmp ;          //calculation of U^(-1) = (FV)^-1 = V^T * F^-1
        assymEigenVec[j][i] = m_eqVector[j] * (*m_eigenvectors)[j][i] ; //calculation of U = FV
        sm += assymEigenVec[j][i] ;
      }
    }

    //------------------------- TEST block ----------------------------------------
    for(size_t i(nchemIdx) ; i<smsize ; ++i){         // multiply U*U^(-1) for testing
      qd_real test = 0.0;
      for(size_t j(nchemIdx) ; j<smsize ; ++j){
        qd_real sm = 0.0;
        for(size_t k(0) ; k<smsize ; ++k){
          sm += assymEigenVec[i][k] * assymInvEigenVec[k][j];
        }
        test += sm;
      }
      if( test < 0.999 || test > 1.001)      // test that U*U^(-1) = 1
        ctest << "row " << i << " of the U*U^(-1) matrix does not equal unity. It sums to " << test << endl;
    }
    //------------------------- TEST block ----------------------------------------
    if (!mFlags.rateCoefficientsOnly){
      qdMatrix Z_matrix(nchem);  // definitions of Y_matrix and Z_matrix taken from PCCP 2007(9), p.4085
      qdb2D Y_matrix;
      Reaction::molMapType::iterator ipos;  // set up an iterator through the isomer map
      Reaction::molMapType::iterator spos;  // set up an iterator through the source map
      sinkMap::iterator sinkpos;           // set up an iterator through the irreversible rxn map

      // check the separation between chemically significant eigenvalues (CSEs)
      // and internal energy relaxation eigenvalues (IEREs); if it's not good, print a warning

      const double last_CSE   = (to_double(m_eigenvalues[nchemIdx]))* m_meanOmega;
      const double first_IERE = (to_double(m_eigenvalues[nchemIdx-1]))* m_meanOmega;
      const double CSE_IERE_separation = to_double(m_eigenvalues[nchemIdx]/m_eigenvalues[nchemIdx-1]);
      if(CSE_IERE_separation > 0.1){
        stringstream ss1 ;
        ss1 << "\nWarning: CSEs not well separated from internal energy relaxation eigenvals (IEREs)" << endl;
        ss1 << "\nThe last CSE = " << last_CSE << " and the first IERE = " << first_IERE << endl;
        ss1 << "(last CSE)/(first IERE) ratio = " << CSE_IERE_separation << ", which is less than an order of magnitude" << endl;
        ss1 << "\nResults obtained from Bartis Widom eigenvalue-vector analysis may be unreliable" << endl;
        ctest << ss1.str() ;
        ppList->XmlWriteValueElement("me:warning", ss1.str());
      }

      int numberOfCemeteries(0); // initialize the number of cemeteries to zero
      for(size_t i(0); i<nchem; ++i){

        numberOfCemeteries = 0; // re-initialize for every nchem calculation.

        // Calculate Z matrix elements for all the isomers in the system.

        for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
          qd_real sm(0.0) ; 
          Molecule* isomer = ipos->first;
          const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
          const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
          int colloptrsize = isomer->getColl().get_colloptrsize() ;       // get colloptrsize for isomer
          colloptrsize += (numberGroupedGrains) ? nrg - numberGroupedGrains : 0 ;
          int rxnMatrixLoc = ipos->second + colloptrsize - 1 ;            // get location for isomer in the rxn matrix
          int seqMatrixLoc = m_SpeciesSequence[isomer];                   // get sequence position for isomer
          for(int j(0) ; j<colloptrsize ; ++j){
            sm += assymEigenVec[rxnMatrixLoc-j][nchemIdx+i];
          }
          Z_matrix[seqMatrixLoc][i] = sm;
        }

        // Calculate Z_matrix matrix elements for all sources in the system.

        for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  
          Molecule* pPseudoIsomer = spos->first ;
          const int rxnMatrixLoc = spos->second;
          const int seqMatrixLoc = m_SpeciesSequence[pPseudoIsomer];
          Z_matrix[seqMatrixLoc][i] = assymEigenVec[rxnMatrixLoc][nchemIdx+i];
        }

        // Calculate Y_matrix elements for sinks.

        if(nsinks) {
          int seqMatrixLoc(0) ;
          for(sinkpos = m_sinkRxns.begin() ; sinkpos != m_sinkRxns.end() ; ++sinkpos, ++seqMatrixLoc) {
            qd_real sm(0.0);
            vector<double> KofEs;                                         // vector to hold sink k(E)s
            vector<double> KofEsTemp;                                     // vector to hold sink k(E)s
            Reaction* sinkReaction = sinkpos->first;
            int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
            if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
              KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
              KofEsTemp.push_back(KofEs[0]);
            } else {                   // if the collision operator size is >1, there are k(E)s for the irreversible loss
              KofEs = sinkReaction->get_GrainKfmc();                      // assign sink k(E)s, the vector size == maxgrn
              Molecule* isomer = sinkReaction->get_reactant();
              const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
              const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
              const int ll = (numberGroupedGrains != 0) ? numberGroupedGrains - nrg : 0 ;
              for (int i(ll) ; i < colloptrsize ; ++i)
                KofEsTemp.push_back(KofEs[i]);
              colloptrsize -= ll ;
            }
            int rxnMatrixLoc = sinkpos->second;                               // get sink location
            for(int j(0);j<colloptrsize;++j){
              sm += assymEigenVec[rxnMatrixLoc+j][nchemIdx+i] * KofEsTemp[j];
            }
            Y_matrix[seqMatrixLoc][i] = sm;
          }
        }

        // calculate Y_matrix matrix elements for cemetery states
        for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
          Molecule* isomer = ipos->first;
          if (isomer->getColl().isCemetery()){ // if it is in cemetery state
            qd_real sm = 0.0;
            vector<double> KofEs;                                         // vector to hold sink k(E)s
            KofEs = isomer->getColl().get_GrainKdmc();
            const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
            const int colloptrsize =  isomer->getColl().get_colloptrsize() - numberGroupedGrains;
            // get colloptrsize for isomer
            int rxnMatrixLoc = ipos->second;                                // get location for isomer in the rxn matrix
            int seqMatrixLoc = int(nsinks) + numberOfCemeteries;      // get sequence position for isomer
            for(int j(0);j<colloptrsize;++j){
              sm += assymEigenVec[rxnMatrixLoc+j][nchemIdx+i] * KofEs[j] * m_meanOmega;
            }
            Y_matrix[seqMatrixLoc][i] = sm;
            KofEs.clear();
            ++numberOfCemeteries;
          }
        }
      }

      // Print out Y_matrix for testing.
      if (nsinks + numberOfCemeteries){
        ctest << "Y_matrix:" << endl;
        Y_matrix.print((int)(nsinks) + numberOfCemeteries, (int)(m_SpeciesSequence.size())); 
      }

      qdMatrix Zinv(Z_matrix) ;
      if (nsinks + numberOfCemeteries) {

        // Apply standard inversion method.

        if(Zinv.invertGaussianJordan()){
          cerr << "Inversion of Z_matrix failed.  Matrix before inversion is: ";
          Z_matrix.showFinalBits(nchem);
        }

      } else {

        // Apply Gram-Schmit orthogonalization in order to invert the matrix.
        // This imposes detailed balance at the macroscopic level.
        //
        // SHR 25/Apr/2010 : It remains unclear that this is correct at the time
        // of writting, however for some systems it is difficult to realize mass
        // conservation without it.

        // Decompose the reduced eigenvector matrix.

        qdMatrix Fr(nchem), Fr_inv(nchem) ;
        for(size_t i(0) ; i<nchem ; ++i){
          Fr[i][i]     = sqrt(Z_matrix[i][nchem-1]) ;
          Fr_inv[i][i] = 1.0/Fr[i][i] ;
        }

        qdMatrix Er = Fr_inv * Z_matrix ;

        // Orthogonalize the reduced symmetric eigenvectro matrix.

        Er.GramSchimdt(nchem - 1) ;

        Z_matrix = Fr * Er ;

        // Transpose the orthonormal matrix and form inverse.

        Er.Transpose() ;

        Zinv = Er * Fr_inv ;

      }

      ctest << "\nZ_matrix: ";
      Z_matrix.showFinalBits(nchem, true);

      ctest << endl << "Z_matrix^(-1):" << endl;
      Zinv.showFinalBits(nchem, true);

      qdMatrix Zidentity = Z_matrix * Zinv ;

      ctest << "\nZ_matrix * Z_matrix^(-1) [Identity matrix]:" << endl;
      Zidentity.showFinalBits(nchem, true);

      // Construct phenomenological rate coefficient matrix.

      qdMatrix Egv(nchem) ;
      for (size_t i(0) ; i<nchem ; ++i){
        Egv[i][i] = m_eigenvalues[nchemIdx+i] * m_meanOmega ; 
      } 
      qdMatrix Kr = Z_matrix * Egv * Zinv ;

      ctest << "\nKr matrix:" << endl;
      Kr.showFinalBits(nchem, true);       // print out Kr_matrix

      // Construct loss matrix.

      qdb2D Kp;
      if (nsinks > 0) {
        for(size_t i(0); i != nsinks + numberOfCemeteries; ++i){    // calculate Kp (definition taken from PCCP 2007(9), p.4085)
          for(size_t j(0);j<nchem;++j){
            qd_real sm = 0.0;
            for(size_t k(0);k<nchem;++k){
              sm += Y_matrix[i][k] * Zinv[k][j];
            }
            Kp[i][j] = sm;
          }
        }
        ctest << "\nKp matrix:" << endl;    // print out Kp_matrix
        Kp.print(nsinks + numberOfCemeteries, m_SpeciesSequence.size());
      }

      // Write out phenomenological rate coefficients.
      PrintPhenomenologicalRates(Kr, Kp, numberOfCemeteries, mFlags, ppList) ;

      mesmerRates = Kr;
    }
    return true;    

  }

  // Write out phenomenological rate coefficients.
  bool CollisionOperator::PrintPhenomenologicalRates(qdMatrix& Kr, qdb2D& Kp, int numberOfCemeteries, MesmerFlags& mFlags, PersistPtr ppList) {

    Reaction::molMapType::iterator ipos;  // set up an iterator through the isomer map

    ctest << "\nFirst order & pseudo first order rate coefficients for loss rxns:\n{\n";
    Reaction::molMapType::iterator lossitr, rctitr, pdtitr;

    stringstream puSymbols;
    stringstream puNumbers;
    // print pseudo 1st order k loss for isomers
    for(lossitr=m_SpeciesSequence.begin(); lossitr!=m_SpeciesSequence.end(); ++lossitr){
      Molecule* iso = lossitr->first;
      int losspos = lossitr->second;
      string isomerName = iso->isCemetery() ? iso->getName() + "(+)" : iso->getName();
      ctest << isomerName << " loss = " << Kr[losspos][losspos] << endl;
      PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderLoss", to_double(Kr[losspos][losspos]));
      ppItem->XmlWriteAttribute("ref", isomerName);

      puNumbers << Kr[losspos][losspos] << "\t";
      if (m_punchSymbolGathered == false){
        puSymbols << isomerName << " loss\t";
      }
    }
    ctest << "}\n";

    if(m_SpeciesSequence.size()>1){
      ctest << "\nFirst order & pseudo first order rate coefficients for isomerization rxns:\n{\n";

      // print pseudo first order connecting ks
      for (rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
        Molecule* rct = rctitr->first;
        string rctName = rct->isCemetery() ? rct->getName() + "(+)" : rct->getName();
        int rctpos = rctitr->second;
        for (pdtitr=m_SpeciesSequence.begin(); pdtitr!=m_SpeciesSequence.end(); ++pdtitr){
          Molecule* pdt = pdtitr->first;
          string pdtName = pdt->isCemetery() ? pdt->getName() + "(+)" : pdt->getName();
          int pdtpos = pdtitr->second;
          if(rctpos != pdtpos){
            ctest << rctName << " -> " << pdtName << " = " << Kr[pdtpos][rctpos] << endl;

            PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", to_double(Kr[pdtpos][rctpos]));
            ppItem->XmlWriteAttribute("fromRef", rctName);
            ppItem->XmlWriteAttribute("toRef",   pdtName);
            ppItem->XmlWriteAttribute("reactionType", "isomerization");
          }

          puNumbers << Kr[pdtpos][rctpos] << "\t";
          if (m_punchSymbolGathered == false){
            puSymbols << rctName << " -> " << pdtName << "\t";
          }
        }
      }
      ctest << "}\n";
    }

    if(m_sinkRxns.size()!=0){
      ctest << "\nFirst order & pseudo first order rate coefficients for irreversible rxns:\n{\n";
      sinkMap::iterator sinkitr = m_sinkRxns.begin();

      for (int sinkpos(0) ; sinkitr!=m_sinkRxns.end() ; ++sinkitr, ++sinkpos) {
        Reaction* sinkReaction = sinkitr->first;          // get Irreversible Rxn
        vector<Molecule*> pdts;
        sinkReaction->get_products(pdts);
        string pdtsName = pdts[0]->getName();
        if (pdts.size() == 2) {pdtsName += + "+"; pdtsName += pdts[1]->getName();}
        for(rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
          Molecule* rcts = rctitr->first;     // get reactants & their position
          int rctpos = rctitr->second;
          if(sinkReaction->getRctColloptrsize()==1){
            ctest << rcts->getName() << " -> "  << pdtsName << "(bim) = " << Kp[sinkpos][rctpos] << endl;
            puNumbers << Kp[sinkpos][rctpos] << "\t";
            if (!m_punchSymbolGathered) {
              puSymbols << rcts->getName() << " -> " << pdtsName << "(bim)\t";
            }
          } else {
            string rctName = rcts->isCemetery() ? rcts->getName() + "(+)" : rcts->getName();

            ctest << rctName << " -> "  << pdtsName << " = " << Kp[sinkpos][rctpos] << endl;

            PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", to_double(Kp[sinkpos][rctpos]));
            ppItem->XmlWriteAttribute("fromRef", rctName);
            ppItem->XmlWriteAttribute("toRef",   pdtsName);
            ppItem->XmlWriteAttribute("reactionType", "irreversible");
            puNumbers << Kp[sinkpos][rctpos] << "\t";
            if (m_punchSymbolGathered == false){
              puSymbols << rctName << " -> " << pdtsName << "\t";
            }
          }
        }
      }
      ctest << "}\n\n";
    }

    Reaction::molMapType::iterator speciesitr;
    int tempNumCemetery(0);
    if(numberOfCemeteries){
      ctest << "\nFirst order & pseudo first order rate coefficients for cemetery states:\n{\n";
      for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
        Molecule* isomer = ipos->first;
        if (isomer->isCemetery()){ // if it is in cemetery state
          string cemName = isomer->getName() + "(-)";
          for (speciesitr=m_SpeciesSequence.begin(); speciesitr!=m_SpeciesSequence.end(); ++speciesitr){
            Molecule* rct = speciesitr->first;
            string rctName = rct->isCemetery() ? rct->getName() + "(+)" : rct->getName();
            int rctpos = speciesitr->second;
            ctest << rctName << " -> " << cemName << " = " << Kp[m_sinkRxns.size()+tempNumCemetery][rctpos] << endl;

            PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", to_double(Kp[m_sinkRxns.size()+tempNumCemetery][rctpos]));
            ppItem->XmlWriteAttribute("fromRef", rctName);
            ppItem->XmlWriteAttribute("toRef",   cemName);
            ppItem->XmlWriteAttribute("reactionType", "deactivation");

            puNumbers << Kp[m_sinkRxns.size()+tempNumCemetery][rctpos] << "\t";
            if (m_punchSymbolGathered == false){
              puSymbols << rctName << " -> " << cemName << "\t";
            }
          }
          ++tempNumCemetery;
        }
      }
      ctest << "}\n\n";
    }

    if (puSymbols.str().size()) {
      puSymbols << "\n";
      mFlags.punchSymbols = puSymbols.str();
      m_punchSymbolGathered = true;
    }

    if (puNumbers.str().size()) {
      puNumbers << "\n";
      mFlags.punchNumbers = puNumbers.str();
    }

    return true;
  }

  //
  // Calculates the Bartis-Widom macroscopic rate coefficients, using the contracted basis set eigenvectors.
  //
  bool CollisionOperator::BartisWidomBasisSetRates(qdMatrix& mesmerRates, MesmerFlags& mFlags) {

    // Constants.
    const size_t smsize   = m_eigenvectors->size() ;
    const size_t nchem    = m_isomers.size() + m_sources.size() ;  // number of isomers+pseudoisomers
    // const size_t nchemIdx = smsize - nchem ;                       // idx for chemically significant eigenvalues & vectors

    // Print out eigenvector matrix.

    //ctest << endl << "Eigenvector matrix:" << endl << endl ;
    //for (size_t i(0) ; i < smsize ; ++i) {
    //  for (size_t j(0) ; j < smsize ; ++j) {
    //    formatFloat(ctest, (*m_eigenvectors)[i][j],  6,  15) ;
    //  }
    //  ctest << endl ;
    //}

    qdMatrix Z(nchem), Zinv(nchem), Kr(nchem);

    if (m_sinkRxns.size()==0){

      //
      // Conservative system.
      //

      // 1. Isomers.

      size_t location(0) ;
      Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
      for (size_t i(0); isomeritr != m_isomers.end() ; ++i, ++isomeritr) {
        location = isomeritr->second ;
        for (size_t j(1); j<=nchem; ++j){
          Z[i][nchem - j] = (*m_eigenvectors)[location][smsize - j] ;
        }
      }

      // Invert Z matrix. 

      ctest << endl << "BW coefficient matrix:" << endl << endl ;
      for (size_t i(0) ; i < nchem ; ++i) {
        for (size_t j(0) ; j < nchem ; ++j) {
          formatFloat(ctest, Z[i][j],  6,  15) ;
          Zinv[j][i] = Z[i][j] ;
        }
        ctest << endl ;
      }

      // Calculate symmetric rate matrix.

      m_eigenvalues[smsize - 1] = 0.0 ;

      for (size_t i(0) ; i < nchem ; ++i) {
        for (size_t j(0) ; j < nchem ; ++j) {
          qd_real sm = 0.0;
          for (size_t k(0) ; k < nchem ; ++k) {
            // sm += Z[i][k] * to_double(m_eigenvalues[nchemIdx+k]) * Zinv[k][j];
            sm += Zinv[i][k]*Z[k][j] ;
          }
          Kr[i][j] = sm ; // * m_meanOmega;
        }
      }

      // Apply similarity transform. 

      //for (size_t i(0) ; i < nchem ; ++i) {
      //  for (size_t j(0) ; j < nchem ; ++j) {
      //    Kr[i][j] *= Z[i][nchem]/Z[j][nchem];
      //  }
      //}

      string rcm(string("Rate coefficient matrix:"));
      Kr.print(rcm, ctest) ;

    } else {

      //
      // Non-conservative system.
      //

    }

    mesmerRates = Kr;

    return true;

  }

  int CollisionOperator::getSpeciesSequenceIndex(const std::string ref)
  {
    Reaction::molMapType::iterator spcitr;
    for (spcitr = m_SpeciesSequence.begin(); spcitr != m_SpeciesSequence.end(); ++spcitr)
    {
      if (ref == (spcitr->first)->getName())
        return spcitr->second;
    }
    cerr << "No molecule named " << ref << " is available in the reaction species.";
    return -1;
  }

  void CollisionOperator::locateSinks()
  {
    m_sinkRxns.clear();
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {

      Reaction* pReaction = (*m_pReactionManager)[i];
      ReactionType reactionType = pReaction->getReactionType() ;

      bool Irreversible = (reactionType == IRREVERSIBLE_ISOMERIZATION || reactionType == IRREVERSIBLE_EXCHANGE || reactionType == DISSOCIATION );
      if (Irreversible && m_sinkRxns.find(pReaction) == m_sinkRxns.end()) {   
        // Add an irreversible rxn to the map.
        Molecule* rctnt = pReaction->get_reactant();
        if(reactionType == IRREVERSIBLE_ISOMERIZATION || reactionType == DISSOCIATION ){
          m_sinkRxns[pReaction] = m_isomers[rctnt];
        } else { // Irreversible exchange reaction.
          m_sinkRxns[pReaction] = m_sources[rctnt];
        }
      }
    }

  }

  // Accessor to get specified eigenvalue.
  double CollisionOperator::getEigenvalue(size_t idEigenvalue) const {

    // Check id is sensible.
    if (idEigenvalue > m_eigenvalues.size()) {
      throw std::runtime_error("Eigenvalue ID greater than collision operator size.");
    }

    return -m_meanOmega*to_double(m_eigenvalues[m_eigenvalues.size()-idEigenvalue]) ;
  }

  // Calculate Yields.
  void CollisionOperator::calculateYields(YieldMap &yieldMap, double &time) const {

    //
    // Yields are calculated by integrating the term Sum_i ki pi.
    // This effectively involves integration of pi between 0 and
    // infinity. This integral leads to the expression Sum_i ki M^(-1)pi_0.
    // The inversion is effected by inverting the eigenvalue expression, 
    // i.e. M^(-1) = FV(eigenvalues)^(-1)V^TF^(-1).
    //
    if(m_sinkRxns.size() == 0){
      // No Sinks so throw an error.
      throw std::runtime_error("No sinks defined, therefore no yields can be calculated.");
    }

    // Get initial distribution.
    size_t smsize = m_eigenvalues.size() ;
    vector<double> p_0(smsize, 0.0) ; 
    if (!produceInitialPopulationVector(p_0)){
      throw std::runtime_error("Calculation of initial conditions vector failed.");
    }

    vector<qd_real> wrk(smsize, 0.0) ;
    for (size_t j(0); j < smsize; ++j) {
      wrk[j] =  p_0[j]/m_eqVector[j] ;
    }

    (*m_eigenvectors).Transpose() ;
    wrk *= (*m_eigenvectors) ;

    if (time > 0.0) {

      // Experimental time.

      for (size_t j(0); j < smsize; ++j) {
        wrk[j] *= (exp(m_meanOmega*m_eigenvalues[j]*time) - 1.0)/(m_meanOmega*m_eigenvalues[j]) ;
      }
    } else {

      // Infinite time limit.

      for (size_t j(0); j < smsize; ++j) {
        wrk[j] /= fabs(m_meanOmega*m_eigenvalues[j]) ;
      }
    }

    (*m_eigenvectors).Transpose() ;
    wrk *= (*m_eigenvectors) ;

    double sum(0.0) ;
    for (size_t j(0); j < smsize; ++j) {
      wrk[j] *= m_eqVector[j] ;
      sum    += to_double(wrk[j]) ;
    }

    sinkMap::const_iterator sinkitr = m_sinkRxns.begin();
    for (; sinkitr != m_sinkRxns.end() ; ++sinkitr) {

      // Locate the sink reaction.

      Reaction* sinkReaction = sinkitr->first ;
      size_t rxnMatrixLoc = sinkitr->second ;

      // Calculate the total flux through this channel.

      // First, determine the mirco rate coefficients for this channel:
      vector<double> ktemp ;
      size_t colloptrsize = sinkReaction->getRctColloptrsize();  
      size_t idx(0) ;
      if(colloptrsize == 1){  
        // If the collision operator size is 1, there is one canonical loss rate coefficient.
        ktemp.push_back(sinkReaction->get_fwdGrnCanonicalRate());
      } else {
        // If the collision operator size is >1, there are k(E)s for the irreversible loss.
        ktemp = sinkReaction->get_GrainKfmc();
        Molecule* isomer = sinkReaction->get_reactant();
        const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
        const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
        idx = (numberGroupedGrains != 0) ? numberGroupedGrains - nrg : 0 ;
        colloptrsize -= idx ;
      }

      // Now form the yield fraction. Note more than one channel may produce the same product.
      double yield(0.0) ;
      for (size_t i(0); i < colloptrsize; ++i) {
        yield += to_double(ktemp[idx + i] * wrk[rxnMatrixLoc + i]) ;
      }
      yieldMap[sinkReaction] = yield ;
    }

  }

  bool CollisionOperator::parseDataForGrainProfileAtTime(PersistPtr ppData)
  {
    //This is called from System::parse()
    do
    {
      const char* pRef = ppData->XmlReadValue("ref", optional);
      if(!pRef)
      {
        cerr << "Need to specify the species with a \"ref\" attribute on me:printGrainProfileAtTime" << endl;
        return false;
      }

      Molecule* pMol = m_pMoleculeManager->find(pRef);
      if(!pMol)
        return false; //error message is in find()
      double tim;
      vector<double> times;
      PersistPtr pp = ppData;
      while( pp = pp->XmlMoveTo("me:time"))
      {
        const char* ptimtxt =pp->XmlRead();
        stringstream ss(ptimtxt);
        ss >>tim;
        times.push_back(tim);
      }
      if(times.empty())
      {
        cerr << "Need to specify at least one time in a \"time\" element in me:printGrainProfileAtTime";
        return false;
      }
      m_GrainProfileAtTimeData.push_back(make_pair(pMol, times));
      //go for next species
      ppData = ppData->XmlMoveTo("me:printGrainProfileAtTime");
    } while(ppData);
    return true;
  }

  bool CollisionOperator::printGrainProfileAtTime(PersistPtr ppAnalysis) {

    // Check there is something to do.
    if (!m_GrainProfileAtTimeData.size())
      return true ;

    // Use GrainProfileAtTimeData to calculate population
    // at each grain energy of each pMol at each time (Struan)

    size_t smsize = m_eigenvectors->size();
    vector<double> r_0(smsize, 0.); // initial distribution
    if (!projectedInitialDistrbtn(r_0)) {
      cerr << "Projection of initial disttribution failed.";
      return false;
    }

    // Copy full eigenvectors of the system.
    dMatrix totalEigenVecs(smsize); 
    for (size_t i(0) ; i < smsize; ++i) {
      double tmp = to_double(m_eqVector[i]);
      for (size_t j(0) ; j < smsize; ++j) {
        totalEigenVecs[i][j] = tmp*to_double((*m_eigenvectors)[i][j]);
      }
    }

    PersistPtr ppGrainList = ppAnalysis->XmlWriteElement("me:grainPopulationList");
    // Iterate over species requested for output
    for (size_t iMol(0); iMol < m_GrainProfileAtTimeData.size(); ++iMol) { 

      // Find the location of the species in the density vector.
      Molecule*  pMol = m_GrainProfileAtTimeData[iMol].first ;
      int iLoc(-1);
      size_t slsize(0);
      if (m_isomers.find(pMol) != m_isomers.end()) {
        iLoc   = m_isomers[pMol] ;
        slsize = pMol->getColl().get_colloptrsize(); ; 
      } else if (m_sources.find(pMol) != m_sources.end()) {
        iLoc = m_sources[pMol] ; 
        slsize = 1 ; 
      } else {
        cerr << "Could not calculate species profile for " << pMol->getName() << "." << endl;
        continue;
      }

      const vector<double> Times(m_GrainProfileAtTimeData[iMol].second) ;

      for (size_t iTime(0); iTime < Times.size(); ++iTime){
        double numColl = m_meanOmega * Times[iTime];
        vector<double> p_t(smsize,0.0) ;

        // |p_t> = exp(Lambda*t)*V^(T)*|init> = exp(Lambda*t)*U^(-1)*|n_0>
        for (size_t j(0) ; j < smsize; ++j) {
          p_t[j] = r_0[j] * exp(to_double(m_eigenvalues[j]) * numColl);
        } 

        // |p_t> =  F*V*exp(Lambda*t)*V^(T)*|init> = U*exp(Lambda*t)*U^(-1)*|n_0> 

        p_t *= totalEigenVecs ;

        // Copy densities for output.

        vector<double> density(p_t.begin() + iLoc, p_t.begin() + (iLoc + slsize - 1)) ;

        // Output density to XML (Chris)
        PersistPtr ppGrainPop = ppGrainList->XmlWriteElement("me:grainPopulation");
        { 
          ppGrainPop->XmlWriteAttribute("ref", pMol->getName());
          ppGrainPop->XmlWriteAttribute("time", toString(Times[iTime]));
          ppGrainPop->XmlWriteAttribute("logTime", toString(log10(Times[iTime])));
          ppGrainPop->XmlWriteAttribute("units", "cm-1");
          for(size_t j(0); j < slsize-1; ++j)  
          {
            PersistPtr ppGrain = ppGrainPop->XmlWriteValueElement("me:grain", density[j], 6);
            ppGrain->XmlWriteAttribute("energy", toString((j+0.5) * pMol->getEnv().GrainSize)); //cm-1
          }
        }
      }
    }

    return true;
  }

  bool CollisionOperator::projectedInitialDistrbtn(vector<double>& r_0) const {

    // This method calculates the projection of the initial distribution on to the
    // eigenspace of the collision matrix.

    vector<double> n_0 = r_0 ; 
    if (!produceInitialPopulationVector(n_0)){
      cerr << "Calculation of initial conditions vector failed.";
      return false;
    }

    // Convert the initial population vector into Boltzmann weighted population vector.
    // All transitions in the reaction matrix are Boltzmann weighted for symmetry.
    // |n_0> = F^(-1)*|n_0>
    for (size_t j(0) ; j < n_0.size() ; ++j) {
      n_0[j] /= to_double(m_eqVector[j]) ;
    }

    // Multiply the initial population with the inverse of the eigenvector
    // which converts the populations into the "decay modes" domain.
    // |r_0> = V^(T)*F^(-1)*|n_0> = U^(-1)*|n_0>
    for (size_t i(0) ; i < r_0.size() ; ++i) {
      double sum = 0.;
      for (size_t j(0) ; j < r_0.size() ; ++j) {
        sum += n_0[j] * to_double((*m_eigenvectors)[j][i]);
      }
      r_0[i] = sum;  
    }

    return true;
  }

}  //namespace
