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
    m_pSystemCollisionOperator(0),
    m_eigenvalues(),
    m_isomers(),
    m_sources(),
    m_initialPopulations(),
    m_meanOmega(0.0)
  {};

  //
  // Add a new reaction to the map.
  //
  bool ReactionManager::addreactions(PersistPtr ppReacList, const MesmerEnv& Env)
  {
    PersistPtr ppReac = ppReacList->XmlMoveTo("reaction");
    while(ppReac)
    {
      //Read reaction ID

      const char* id = ppReac->XmlReadValue("id");
      if(!id)
        cinfo << "Reaction ID not found.\n";

      // Read reactant and product types.

      string rct1Name, rct1Type, rct2Name, rct2Type ;
      string pdt1Name, pdt1Type, pdt2Name, pdt2Type ;
      bool readStatus(true), bRct2(false), bPdt1(false), bPdt2(false) ;

      PersistPtr ppReactant1  = ppReac->XmlMoveTo("reactant");
      readStatus = GetMoleculeInfo(ppReactant1, rct1Name, rct1Type) ;

      PersistPtr ppReactant2  = ppReactant1->XmlMoveTo("reactant");
      if(ppReactant2) {
        readStatus = (readStatus && GetMoleculeInfo(ppReactant2, rct2Name, rct2Type)) ;
        bRct2 = true;
      }

      PersistPtr ppProduct1 = ppReac->XmlMoveTo("product");
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
      // Create a new Reaction.
      //
      Reaction *preaction ;
      if     (!bRct2 && bPdt1 && pdt1Type == "modelled" && !bPdt2)
        preaction = new IsomerizationReaction(m_pMoleculeManager, Env, id) ;
      else if( bRct2 && bPdt1 && !bPdt2)
        preaction = new AssociationReaction(m_pMoleculeManager, Env, id) ;
      else if(!bRct2 && bPdt1 && (pdt1Type == "sink" || pdt2Type == "sink"))
        preaction = new IrreversibleReaction(m_pMoleculeManager, Env, id) ;
      else if( bRct2 && bPdt1 &&  bPdt2)
        preaction = new ExchangeReaction(m_pMoleculeManager, Env, id) ;
      else {
        cinfo << "Unknown reaction type.\n";
        return false ;
      }
      /*The information of the products of a dissociation reaction is necessary, as in the xml output, Mesmer needs to know
      the products to draw the potential energy surface. In addition, for dissociation reaction with QM tunneling,
      Mesmer also needs to know the barrier height on the products side. */

      //
      // Initialize Reaction from input stream.
      //
      if(!preaction->InitializeReaction(ppReac)){
        delete preaction;
        return false;
      }

      //
      // Add reaction to map.
      //

      //need to check if there is duplicate reaction name/species: CHL

      m_reactions.push_back(preaction) ;

      ppReac = ppReac->XmlMoveTo("reaction");
    }

    return true;
  }

  void ReactionManager::resetCalcFlags(){
    for (size_t i(0) ; i < size() ; ++i) {
      m_reactions[i]->resetCalcFlag();
    }
  }

  bool ReactionManager::SetGrainParams(MesmerEnv &Env, const double minEne, const double maxEne)
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

    Env.EMin = minEne;
    Env.EMax = maxEne;

    /*For testing purposes, set the maxGrn based on the highest temperature we use in all calculations.*/
    double MaximumTemperature = Env.MaximumTemperature;

    /*EAboveHill: Max energy above the highest hill. The temperature refers to the current condition.*/
    if (Env.useTheSameCellNumber){
      Env.EMax += Env.EAboveHill * MaximumTemperature * boltzmann_RCpK;
    }
    else{
      Env.EMax += Env.EAboveHill / Env.beta;
    }

    if(Env.GrainSize <= 0.0){
      Env.GrainSize = 100; //default 100cm-1
      cerr << "Grain size was invalid. Reset grain size to default: 100";
    }

    Env.MaxGrn = (int)((Env.EMax-Env.EMin)/Env.GrainSize + 0.5);
    Env.MaxCell = Env.GrainSize * Env.MaxGrn;

    cerr << "Cell number = " << Env.MaxCell << ", Grain number = " << Env.MaxGrn << endl;

    return true;
  }

  bool ReactionManager::BuildSystemCollisionOperator(MesmerEnv &Env)
  {
    // reset the DOS calculation flags before building the system collision operator
    resetCalcFlags();

    //
    // Find all the unique wells and lowest zero point energy.
    //
    m_isomers.clear();

    double minEnergy = 0.0 ; //this is the minimum of ZPE amongst all wells
    double maxEnergy = 0.0 ; //this is the maximum of ZPE amongst all hills
    double populationSum = 0.0;
    BathGasMolecule *pBathGasMolecule = dynamic_cast<BathGasMolecule*>(m_pMoleculeManager->get_BathGasMolecule());

    // populate isomerMap with unimolecular species and determine minimum/maximum energy on the PES
    for (size_t i(0) ; i < size() ; ++i) {
      vector<ModelledMolecule *> unimolecules ;
      m_reactions[i]->get_unimolecularspecies(unimolecules) ;

      // populate isomerMap with unimolecular species
      for (size_t j(0) ; j < unimolecules.size() ; ++j) {
        // wells
        CollidingMolecule *pCollidingMolecule = dynamic_cast<CollidingMolecule*>(unimolecules[j]) ;
        if(m_isomers.find(pCollidingMolecule) == m_isomers.end()){ // New isomer
          double population = m_initialPopulations[pCollidingMolecule->getName()];
          if (population){
            populationSum += population;
            pCollidingMolecule->setInitPopulation(population);
            ctest << "Initial population of " << pCollidingMolecule->getName() << " = " << population << endl;
          }
          m_isomers[pCollidingMolecule] = 0 ; //initialize to a trivial location
          minEnergy = min(minEnergy, pCollidingMolecule->get_zpe()) ;
          maxEnergy = max(maxEnergy, pCollidingMolecule->get_zpe()) ;
        }
      }

      // SuperMolecules
      // first check for any SuperMolecule in this reaction
      SuperMolecule* pSuper = m_reactions[i]->get_bi_molecularspecies();
      if (pSuper){
        minEnergy = min(minEnergy, pSuper->get_zpe()) ;
        maxEnergy = max(maxEnergy, pSuper->get_zpe()) ;
      }

      // Transition State
      // third check for the transition state in this reaction
      TransitionState *pTransitionState = m_reactions[i]->get_TransitionState();
      if (pTransitionState){
        minEnergy = min(minEnergy, pTransitionState->get_zpe()) ;
        maxEnergy = max(maxEnergy, pTransitionState->get_zpe()) ;
      }
    }

    // set grain parameters for the current Temperature/pressure condition
    if(!SetGrainParams(Env, minEnergy, maxEnergy))
      return false;

    // Set grain ZPE for all species involved in the reactions according to the minimum energy of the system.
    for (size_t i(0) ; i < size() ; ++i) {
      // first check for any SuperMolecule in this reaction
      SuperMolecule* pSuper = m_reactions[i]->get_bi_molecularspecies();
      // the grain ZPE of SuperMolecule has to be calculated from zpeReactant1 + zpeReactant2 - minEnergy
      if (pSuper){
        double zpe = pSuper->get_relative_ZPE();
        pSuper->set_grainValues(zpe);
      }

      // second check for unimolecular species in this reaction
      std::vector<ModelledMolecule *> unimolecules ;
      m_reactions[i]->get_unimolecularspecies(unimolecules) ;
      for (size_t j(0) ; j < unimolecules.size() ; ++j) {
        CollidingMolecule *pCollidingMolecule = dynamic_cast<CollidingMolecule*>(unimolecules[j]) ;
        double zpe = pCollidingMolecule->get_relative_ZPE() ;
        pCollidingMolecule->set_grainValues(zpe);
      }
    }

    // Caluclate TSFlux and k(E)s
    for (size_t i(0) ; i < size() ; ++i) {
      m_reactions[i]->calcGrnAvrgMicroRateCoeffs() ;
    }

    if (!Env.rateCoefficientsOnly){
      //
      // Shift all wells to the same origin, calculate the size of the system collision operator,
      // calculate the mean collision frequency and initialize all collision operators.
      //
      int msize(0) ; // size of the collision matrix
      m_meanOmega = 0.0;

      Reaction::isomerMap::iterator isomeritr = m_isomers.begin() ;
      for (; isomeritr != m_isomers.end() ; ++isomeritr) {

        CollidingMolecule *isomer = isomeritr->first ;
        isomeritr->second = msize ; //set location

        int grnZpe = isomer->get_grnZpe() ; //set grain ZPE (with respect to the minimum of all wells)

        int colloptrsize = Env.MaxGrn - grnZpe ;
        isomer->set_colloptrsize(colloptrsize) ;
        msize += colloptrsize ;

        isomer->initCollisionOperator(Env.beta, pBathGasMolecule) ;
        m_meanOmega += isomer->get_collisionFrequency() ;
      }
      m_meanOmega /= m_isomers.size();

      //
      // Find all source terms.
      // Note: 1. A source term contains the deficient reactant.  It is possible for
      //          there to be more than one source term.
      //       2. In the current construction of Mesmer, the source is a SuperMolecule
      //          representing both reactants.
      m_sources.clear(); // Maps the location of source in the system matrix.
      for (size_t i(0) ; i < size() ; ++i) {
        AssociationReaction *pReaction = dynamic_cast<AssociationReaction*>(m_reactions[i]) ;
        if (pReaction) {
          SuperMolecule *pSuperMolecule = pReaction->get_bi_molecularspecies();
          if (pSuperMolecule && m_sources.find(pSuperMolecule) == m_sources.end()){ // New source
            double population = m_initialPopulations[(pSuperMolecule->getMember1())->getName()];
            if (population){
              populationSum += population;
              (pSuperMolecule->getMember1())->setInitPopulation(population);
            }
            if (populationSum == 0.0){
              populationSum += 1.0;
              (pSuperMolecule->getMember1())->setInitPopulation(populationSum);
            }
            m_sources[pSuperMolecule] = msize ;
            pReaction->putSourceMap(&m_sources) ;
            ++msize ;
          }
        }
      }

      // Allocate space for the full system collision operator.
      m_pSystemCollisionOperator = new dMatrix(msize) ;

      // Insert collision operators for individual wells.
      for (isomeritr = m_isomers.begin() ; isomeritr != m_isomers.end() ; ++isomeritr) {

        CollidingMolecule *isomer = isomeritr->first ;
        int colloptrsize = isomer->get_colloptrsize() ;
        double omega = isomer->get_collisionFrequency() ;
        int idx = isomeritr->second ;

        isomer->copyCollisionOperator(m_pSystemCollisionOperator, colloptrsize, idx, omega/m_meanOmega) ;

      }

      // Add connecting rate coefficients.
      for (size_t i(0) ; i < size() ; ++i) {
        m_reactions[i]->AddReactionTerms(m_pSystemCollisionOperator,m_isomers,1.0/m_meanOmega) ;
      }
    }

    return true;
  }

  void ReactionManager::diagCollisionOperator(const MesmerEnv &Env)
  {
    // Allocate space for eigenvalues.
    const int smsize = int(m_pSystemCollisionOperator->size()) ;

    m_eigenvalues.clear();
    m_eigenvalues.resize(smsize, 0.0);
    m_pSystemCollisionOperator->diagonalize(&m_eigenvalues[0]) ;

    int numberStarted = 0;
    int numberPrinted = smsize; // Default prints all of the eigenvalues
    if (Env.printEigenValuesNum > 0 && Env.printEigenValuesNum <= smsize){ //at least prints 1 eigenvalue
      numberPrinted = Env.printEigenValuesNum;
      numberStarted = smsize - Env.printEigenValuesNum;
    }

    ctest << "\nTotal number of eigenvalues = " << smsize << endl;
    ctest << "Eigenvalues\n{\n";
    for (int i = numberStarted ; i < smsize; ++i) {
      formatFloat(ctest, m_meanOmega * m_eigenvalues[i], 6, 15) ;
      ctest << endl ;
    }
    ctest << "}\n";
  }

  //
  // Extract molecule inforamtion from XML stream.
  //
  bool ReactionManager::GetMoleculeInfo(PersistPtr pp, string& MolName, string& MolType)
  {
    PersistPtr ppmol = pp->XmlMoveTo("molecule");
    if(!ppmol) {
      cerr << "Ill formed molecule tag." << endl;
      return false;
    }
    MolName = ppmol->XmlReadValue("ref");
    if(MolName.size()){
      MolType = ppmol->XmlReadValue("me:type");
    } else {
      cerr << "Cannot find molecule name." << endl;
    }

    return true ;
  }

  bool ReactionManager::calculateEquilibriumFractions(const double beta)
  { /* Consider a three well system: e.g., A <-> B <-> C where A <-> B has Keq = K1 & B <-> C has Keq = K2.
       This routine uses the fact that the normalized equilibrated system may be described
       by a 3x3 matrix and a vector which satisfy the following:
                                            |-K1  1   0| |A|   |0|
                                            | 0  -K2  1| |B| = |0|
                                            | 1   1   1| |C|   |1|
       The equilibrium fraction of each isomer (or psuedo isomer, in the case of a source term) may be
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

      // check if a reaction is association or isomerization
      AssociationReaction* arc = dynamic_cast<AssociationReaction*>(m_reactions[i]);
      IsomerizationReaction* irc = dynamic_cast<IsomerizationReaction*>(m_reactions[i]);
      ModelledMolecule* rct;
      ModelledMolecule* pdt;
      long double Keq;

      if(arc){     //if it's an association reaction
        rct = arc->get_pseudoIsomer();        // get the reactants
        vector<ModelledMolecule*> unispecies;
        arc->get_unimolecularspecies(unispecies);
        pdt = unispecies[0];                  // get the products
        Keq = arc->calcEquilibriumConstant();
      }
      else if(irc){     //if it's an isomerization reaction
        vector<ModelledMolecule*> isomers;
        irc->get_unimolecularspecies(isomers);
        rct = isomers[0];                     // get the reactants
        pdt = isomers[1];                     // get the products
        Keq = irc->calcEquilibriumConstant();
      }

      if (arc || irc){                  //only need eq fracs for species in isom & assoc rxns

        int ploc, rloc, rval, pval;

        map<ModelledMolecule*, int>::iterator rctitr = m_SpeciesSequence.find(rct);   //check if the reactant is in the map
        if(rctitr==m_SpeciesSequence.end())        //if the reactant isnt in the map
          rval = 0;
        else{
          rloc = rctitr->second;        //if the reactant is in the map, get the location
          rval = 1;
        }

        map<ModelledMolecule*, int>::iterator pdtitr = m_SpeciesSequence.find(pdt);   //check if the product is in the map
        if(pdtitr==m_SpeciesSequence.end())        //if the product isnt in the map
          pval = 0;
        else{
          ploc = pdtitr->second;        //if the product is in the map, get the location
          pval = 1;
        }

        if(rval==0 && pval==0){             // if neither reactant nor product are in the m_SpeciesSequence map
          m_SpeciesSequence[rct] = counter;            // update the eqMatrix elements
          counter += 1;
          m_SpeciesSequence[pdt] = counter;
          eqMatrix[counter-1][counter-1] =+ -1.0 * Keq;
          eqMatrix[counter-1][counter] =+ 1.0;
          counter += 1;
        }
        else if(rval==0 && pval==1){        // if reactant isnt in m_SpeciesSequence map & product is
          m_SpeciesSequence[rct] = counter;            // update the eqMatrix matrix elements
          eqMatrix[counter-1][ploc] =+ 1.0;
          eqMatrix[counter-1][counter] =+ -1.0 * Keq;
          counter += 1;
        }
        else if(rval==1 && pval==0){        // if reactant is in m_SpeciesSequence map & product isnt
          m_SpeciesSequence[pdt] = counter;            // update the eqMatrix matrix elements
          eqMatrix[counter-1][rloc] =+ -1.0 * Keq;
          eqMatrix[counter-1][counter] =+ 1.0;
          counter += 1;
        }
        else if(rval==1 && pval==1){        // if both reactant & product are in m_SpeciesSequence map

          double pdtRowSum(0.0), rctRowSum(0.0);

          for(int j(0);j<counter;++j){           // calculate pdt & rct rowSums of EqMatrix to see if the rxn is redundant
            pdtRowSum += eqMatrix[ploc][j];
            rctRowSum += eqMatrix[rloc][j];
          }

          if(pdtRowSum!=0.0 && rctRowSum!=0.0){ // connection is redundant
            eqMatrix[counter-1][ploc] =+ 1.0;
            eqMatrix[counter-1][rloc] =+ -1.0 * Keq;
          }
          else if(rctRowSum==0.0){              // connection is not redundant, pdts lack specification
            eqMatrix[rloc][ploc] =+ 1.0;
            eqMatrix[rloc][rloc] =+ -1.0 * Keq;
          }
          else if(pdtRowSum==0.0){
            eqMatrix[ploc][ploc] =+ 1.0;        // connection is not redundant, rcts lack specification
            eqMatrix[ploc][rloc] =+ -1.0 * Keq;
          }
        }
      }
    }
    for(int i=0; i < counter; ++i){         // add ones to the final row of the matrix
      eqMatrix[counter-1][i]= 1.0;
    }

    //    ctest << "matrix elements for calculating isomer equilibrium fractions:" << endl;
    //    eqMatrix.showFinalBits(counter);

    dMatrix backup(eqMatrix);  //backup EqMatrix for error reporting

    if(eqMatrix.invert()){
      cerr << "Inversion of matrix for calculating Eq fractions failed.  Matrix before inversion is: ";
      backup.showFinalBits(counter);
    }

    //    ctest << "inverse:" << endl;
    //    eqMatrix.showFinalBits(counter);

    map<ModelledMolecule*, int>::iterator itr1;

    for(itr1= m_SpeciesSequence.begin(); itr1!=m_SpeciesSequence.end(); ++itr1){  //assign Eq fraction to appropriate ModelledMolecule
      int position = itr1->second;                          //in the Eq frac map
      ModelledMolecule* key = itr1->first;
      key->setEqFraction(eqMatrix[position][counter-1]);    //set Eq fraction to last column in eqMatrix
      ctest << "Equilibrium Fraction for " << key->getName() << " = " << key->getEqFraction() << endl;
    }
    return true;
  }

  bool ReactionManager::produceInitialPopulationVector(vector<long double>& initDist){

    long double populationSum = 0.0;

    Reaction::isomerMap::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
      CollidingMolecule* isomer = ipos->first;                        // to get isomer initial populations
      populationSum += isomer->getInitPopulation();
    }

    Reaction::sourceMap::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map to get
      SuperMolecule* source = spos->first;                            // source initial populations
      populationSum += (source->getMember1())->getInitPopulation();
    }

    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      CollidingMolecule* isomer = ipos->first;                        // get initial population of each isomer
      long double initFrac = isomer->getInitPopulation();
      if (initFrac != 0.0){                                           // if isomer initial populations are nonzero
        initFrac /= populationSum;                                    // normalize initial pop fraction
        int location = ipos->second;
        const int colloptrsize = isomer->get_colloptrsize();
        vector<long double> boltzFrac;
        isomer->normalizedGrainDistribution(boltzFrac, colloptrsize);
        for (int i = 0; i < colloptrsize; ++i){
          initDist[i + location] = initFrac * boltzFrac[i];
        }
      }
    }

    // if there is no source term and the populationSum is still zero, set population = 1.0 for the first isomer
    int sizeSource = static_cast<int>(m_sources.size());
    if (populationSum == 0. && sizeSource == 0){
      ipos = m_isomers.begin();
      CollidingMolecule* isomer = ipos->first;
      isomer->setInitPopulation(1.0); // set initial population for the first isomer
      long double initFrac = isomer->getInitPopulation();
      cinfo << "No population was assigned, and there is no source term."  << endl
            << "Initialize a Boltzmann distribution in the first isomer." << endl;
      int location = ipos->second;
      const int colloptrsize = isomer->get_colloptrsize();
      vector<long double> boltzFrac;
      isomer->normalizedGrainDistribution(boltzFrac, colloptrsize);
      for (int i = 0; i < colloptrsize; ++i){
        initDist[i + location] = initFrac * boltzFrac[i];
      }
    }

    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      SuperMolecule* source = spos->first;
      long double initFrac = (source->getMember1())->getInitPopulation() / populationSum;
      int location = spos->second;
      initDist[location] = initFrac;
      if (populationSum == 0. && spos == m_sources.begin()){
        cinfo << "No population was assigned. Initialize the first source term to 1.0." << endl;
        initDist[location] = 1.0;
      }
    }

    return true;
  }

  bool ReactionManager::produceEquilibriumVector()
  { //the vector produced by this function contains the sqrt of the normalized equilibrium distribution

    eqVector.clear();
    eqVector.resize(m_pSystemCollisionOperator->size());

    Reaction::sourceMap::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map to get
      SuperMolecule* source = spos->first;                            // eq Fractions
      int location = spos->second;
      long double eqFrac = (source->getMember1())->getEqFraction();
      eqVector[location] = sqrt(eqFrac);
    }

    Reaction::isomerMap::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
      CollidingMolecule* isomer = ipos->first;                        // to get eq Fractions
      int location = ipos->second;
      long double eqFrac = isomer->getEqFraction();
      const int colloptrsize = isomer->get_colloptrsize();
      vector<long double> boltzFrac;
      isomer->normalizedGrainDistribution(boltzFrac, colloptrsize);
      for(int i(0);i<colloptrsize;++i){
        eqVector[location + i]= sqrt(eqFrac * boltzFrac[i]);
      }
    }
    return true;
  }

  bool ReactionManager::timeEvolution(int maxTimeStep, const MesmerEnv mEnv)
  {
    int smsize = int(m_pSystemCollisionOperator->size());

    long double maxEvoTime = 0.;
   // set the default maximum evolution time
    if (mEnv.maxEvolutionTime <= 0. || mEnv.maxEvolutionTime > 1.0e8)
      maxEvoTime = 1.0e8;
    else
      maxEvoTime = mEnv.maxEvolutionTime;

    /* calculate the time points */
    vector<double> timePoints;
    for (int i = 0; i < maxTimeStep; ++i){
      long double time = pow(10., static_cast<double>(i) / 10. - 11.);
      if (time > maxEvoTime)
        break;
      timePoints.push_back(time);
    }

    vector<long double> initDist(smsize, 0.); // initial distribution
    if (!produceInitialPopulationVector(initDist)){
      cerr << "Calculation of initial conditions vector failed.";
      return false;
    }

    if (!produceEquilibriumVector()){
      cerr << "Calculation of equilibrium vector failed.";
      return false;
    }

    // |initDist> = F^(-1)*|p(0)>
    for (int j = 0; j < smsize; ++j) {
      initDist[j] = initDist[j]/eqVector[j];
    }

    dMatrix sysCollOptr(*m_pSystemCollisionOperator); // copy the system collision operator, which holds the eigenvectors
    vector<long double> work1(smsize, 0.);

    for (int i = 0; i < smsize; ++i) {
      long double sum = 0.;
      for (int j = 0; j < smsize; ++j) {
        sum += initDist[j] * sysCollOptr[j][i];
      }
      work1[i] = sum;  // now |work1> = V^(T)*|init> = U^(-1)*|p(0)>
    }

    for (int i = 0; i < smsize; ++i) {
      long double tmp = eqVector[i];
      for (int j = 0; j < smsize; ++j) {
        sysCollOptr[i][j] *= tmp;
      }
    }

    ldb2D grainedProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
    vector<long double> work2(smsize, 0.);

    for (int timestep = 0; timestep < maxTimeStep; ++timestep){
      long double numColl = m_meanOmega * timePoints[timestep];
      for (int j = 0; j < smsize; ++j) {
        work2[j] = work1[j] * exp(m_eigenvalues[j] * numColl);
      } // now |wk2> = exp(Dt)*V^(T)*|init> = exp(Dt)*U^(-1)*|p(0)>
      for (int j = 0; j < smsize; ++j) {
        long double sum = 0.;
        for (int l = 0; l < smsize; ++l) {
          sum += work2[l] * sysCollOptr[j][l];
        }
        grainedProfile[j][timestep] = sum;
      } // now |grainedProfile(t)> = |grainedProfile(i)> = F*V*exp(Dt)*V^(T)*|init> = U*exp(Dt)*U^(-1)*|p(0)>
    }

    vector<long double> totalIsomerPop(maxTimeStep, 0.);
    vector<long double> totalProductPop(maxTimeStep, 0.);

    for(int timestep(0); timestep<maxTimeStep; ++timestep){
      for(int j(0);j<smsize;++j){
        totalIsomerPop[timestep] += grainedProfile[j][timestep];
      }
      totalProductPop[timestep] = 1.0 - totalIsomerPop[timestep];
    }

    //------------------------------
    // print grained species profile
    if (mEnv.grainedProfileEnabled) {
      ctest << "\nGrained species profile (the first row is time points in unit of second):\n{\n";
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        formatFloat(ctest, timePoints[timestep], 6,  15);
      }
      ctest << endl;
      for (int j = 0; j < smsize; ++j) {
        for (int timestep = 0; timestep < maxTimeStep; ++timestep){
          formatFloat(ctest, grainedProfile[j][timestep], 6,  15);
        }
        ctest << endl;
      }
      ctest << "}\n";
    }
    //------------------------------

    int numberOfSpecies = static_cast<int>(m_isomers.size() + m_sources.size());
    db2D speciesProfile(numberOfSpecies, maxTimeStep);

    //----------------------
    // print species profile
    ctest << "\nSpecies profile (the first row is time points in unit of second):\n{\nTimesteps   ";
    for (int timestep = 0; timestep < maxTimeStep; ++timestep){
      formatFloat(ctest, timePoints[timestep], 6,  15);
    }
    ctest << endl;
    //-----------------------------------------------
    // speciesProfile will contain the sum of all grains corresponding to an individual species at each time step
    // and is sorted so that it has the same ordering as the system collision operator
    Reaction::isomerMap::iterator ipos;
    std::map<int, CollidingMolecule*> numMap1;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      numMap1[ipos->second] = ipos->first;
    }
    std::map<int, CollidingMolecule*>::iterator jpos;
    int j = 0;
    for (jpos = numMap1.begin(); jpos != numMap1.end(); ++jpos, ++j){
      int idx = jpos->first;
      int colloptrsize = (jpos->second)->get_colloptrsize();
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        for (int i = 0; i < colloptrsize; ++i){
          speciesProfile[j][timestep] += grainedProfile[i + idx][timestep];
        }
      }
    }

    // print isomer terms
    for (jpos = numMap1.begin(), j = 0; jpos != numMap1.end(); ++jpos, ++j){
      ctest << (jpos->second)->getName() << " ";
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        formatFloat(ctest, speciesProfile[j][timestep], 6,  15);
      }
      ctest << endl;
    }

    // printing source terms
    Reaction::sourceMap::iterator kpos;
    std::map<int, SuperMolecule*> numMap2;
    for (kpos = m_sources.begin(); kpos != m_sources.end(); ++kpos){
      numMap2[kpos->second] = kpos->first;
    }
    std::map<int, SuperMolecule*>::iterator lpos;
    for (lpos = numMap2.begin(); lpos != numMap2.end(); ++lpos, ++j){
      int idx = lpos->first;
      ctest << (lpos->second)->getName() << " ";
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        speciesProfile[j][timestep] = grainedProfile[idx][timestep];
        formatFloat(ctest, speciesProfile[j][timestep], 6,  15);
      }
      ctest << endl;
    }
    ctest << "}\n";

    //-----------------------------------------------
    return true;
  }

  // Set Initial population for individual species
  void ReactionManager::setInitialPopulation(PersistPtr anchor)
  {
    PersistPtr pp=anchor;
    const char* txt;
    PersistPtr ppInitMol = pp->XmlMoveTo("molecule");
    while(ppInitMol){
      double population = 0.0;
      string sRef = ppInitMol->XmlReadValue("ref");
      if(sRef.size()){ // if got the name of the molecule
        txt = ppInitMol->XmlReadValue("me:population");
        population = atof(txt);
        m_initialPopulations[sRef] = population;
      }
      ppInitMol = ppInitMol->XmlMoveTo("molecule");
    }
  }

  bool ReactionManager::BartisWidomPhenomenologicalRates()
  {
    dMatrix eigenVec(*m_pSystemCollisionOperator);  //copy SystemCollisionOperator, the eigenvector Matrix (== V)
    int smsize = int(m_pSystemCollisionOperator->size());
    dMatrix assymInvEigenVec(smsize);   // U^(-1)
    dMatrix assymEigenVec(smsize);      // U
    dMatrix EigenVecIdentity(smsize);   // matrix for holding product of U^(-1) * U
    double tmp;
    double sm, test;
    map<ModelledMolecule*, int> SpeciesSequence;  //initialize a map to keep track of species sequence in matrices

    int nchem = static_cast<int>(m_isomers.size() + m_sources.size());  // number of isomers+psuedoisomers
    int nchemIdx = smsize - nchem;       // idx for chemically significant eigenvalues & vectors

    for(int i(0);i<smsize;++i){
      tmp = eqVector[i];
      for(int j(0);j<smsize;++j){
        assymInvEigenVec[j][i] = eigenVec[i][j]/tmp;         //calculation of U^(-1) = (FV)^-1 = V^T * F^-1
        assymEigenVec[j][i] = eqVector[j] * eigenVec[j][i];  //calculation of U = FV
      }
    }

    for(int i(0);i<smsize;++i){          // multiply U*U^(-1) for testing
      test = 0.0;
      for(int j(0);j<smsize;++j){
        sm = 0.0;
        for(int k(0);k<smsize;++k){
          sm += assymEigenVec[i][k] * assymInvEigenVec[k][j];
        }
        EigenVecIdentity[i][j] = sm;
        test += sm;
      }
      if((test/1.0) < 0.999 || (test/1.0) > 1.001)      // test that U*U^(-1) = 1
        ctest << "row " << i << " of the U*U^(-1) matrix does not equal unity. It sums to " << test << endl;
    }

    // EigenVecIdentity.showFinalBits(nchem);

    dMatrix Z(nchem), Y(nchem);          // definitions of Z & Y taken from PCCP 2007(9), p.4085
    Reaction::isomerMap::iterator ipos;  // set up an iterator through the isomer map
    Reaction::sourceMap::iterator spos;  // set up an iterator through the source map

    for(int i(0); i<nchem; ++i){
      for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // calculate Z matrix elements for
        sm = 0.0;                                                       // all isomers in the system
        CollidingMolecule* isomer = ipos->first;
        const int colloptrsize = isomer->get_colloptrsize();
        int location = ipos->second;
        int position = m_SpeciesSequence[isomer];
        for(int j(0);j<colloptrsize;++j){
          sm += assymEigenVec[location+j][nchemIdx+i];
        }
        Z[position][i] = sm;
      }
      for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // calculate Z matrix elements for
        sm = 0.0;                                                       // all sources in the system
        SuperMolecule* source = spos->first;
        ModelledMolecule* psuedo_isomer = source->getMember1();
        int location = spos->second;
        int position = m_SpeciesSequence[psuedo_isomer];
        sm = assymEigenVec[location][nchemIdx+i];
        Z[position][i] = sm;
      }
    }

//    ctest << "Z matrix:" << endl;
//    Z.showFinalBits(nchem);

    dMatrix Zinv(Z), Zidentity(nchem), Kr(nchem);

    if(Zinv.invert()){
      cerr << "Inversion of Z matrix failed.  Matrix before inversion is: ";
      Z.showFinalBits(nchem);
    }

//    ctest << "inverse of Z:" << endl;
//    Zinv.showFinalBits(nchem);

    for(int i(0);i<nchem;++i){          // multiply Z*Z^(-1) for testing
      for(int j(0);j<nchem;++j){
        sm = 0.0;
        for(int k(0);k<nchem;++k){
          sm += Z[i][k] * Zinv[k][j];
        }
        Zidentity[i][j] = sm;
      }
    }

    ctest << "\nZ*Z^(-1):" << endl;
    Zidentity.showFinalBits(nchem);

    for(int i(0);i<nchem;++i){          // calculate Kr (definition taken from PCCP 2007(9), p.4085)
      for(int j(0);j<nchem;++j){
        sm = 0.0;
        for(int k(0);k<nchem;++k){
          sm += Z[i][k] * m_eigenvalues[nchemIdx+k] * Zinv[k][j];
        }
        Kr[i][j] = sm * m_meanOmega;
      }
    }

    ctest << "\nKr:" << endl;
    Kr.showFinalBits(nchem);

    ctest << "\nFirst order and psuedo first order rate coefficients:\n{\n";
    map<ModelledMolecule*, int>::iterator lossitr;
    map<ModelledMolecule*, int>::iterator rctitr;
    map<ModelledMolecule*, int>::iterator pdtitr;
    // print k loss for isomers
    for(lossitr=m_SpeciesSequence.begin(); lossitr!=m_SpeciesSequence.end(); ++lossitr){
      ModelledMolecule* iso = lossitr->first;
      int losspos = lossitr->second;
      ctest << iso->getName() << " loss = " << Kr[losspos][losspos] << endl;
    }
    // print k for connecting rates
    for (rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){  // print connecting rates
      ModelledMolecule* rct = rctitr->first;
      int rctpos = rctitr->second;
      for (pdtitr=m_SpeciesSequence.begin(); pdtitr!=m_SpeciesSequence.end(); ++pdtitr){
        ModelledMolecule* pdt = pdtitr->first;
        int pdtpos = pdtitr->second;
        if(rctpos != pdtpos)
          ctest << rct->getName() << " -> " << pdt->getName() << " =  " << Kr[pdtpos][rctpos] << endl;
      }
    }
    ctest << "}\n";

    return true;
}

}//namespace



