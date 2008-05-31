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
    m_populations(),
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
      if (!bRct2 && bPdt1 && !bPdt2)
        preaction = new IsomerizationReaction(m_pMoleculeManager, Env, id) ;
      else if(bRct2 && bPdt1 && !bPdt2)
        preaction = new AssociationReaction(m_pMoleculeManager, Env, id) ;
      else if(!bRct2 && bPdt1 && bPdt2)
        preaction = new DissociationReaction(m_pMoleculeManager, Env, id) ;
      else if(bRct2 && bPdt1 && bPdt2)
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
          double population = m_populations[pCollidingMolecule->getName()];
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
            double population = m_populations[(pSuperMolecule->getMember1())->getName()];
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

  bool ReactionManager::calculateEquilibriumFractions(const double beta){
    // All species we care about in here is the source term(s) and isomers.
    // A sink term's equilibrium population is not calculated, because there is no equilibrium for a reaction with a sink term.

    // Loop through all isomers in the system giving their partition function values as their equilibrium fraction.
    Reaction::isomerMap::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      CollidingMolecule* pI1 = ipos->first;
      long double eqFrac(1.0);
      const long double prtFn1 = pI1->rovibronicGrnCanPrtnFn();
      Reaction::isomerMap::iterator jpos;
      for (jpos = m_isomers.begin(); jpos != m_isomers.end(); ++jpos){
        if (jpos->first != pI1){
          CollidingMolecule* pI2 = jpos->first;
          const long double HeatDiff = pI2->get_zpe() - pI1->get_zpe();
          const long double prtFn2 = pI2->rovibronicGrnCanPrtnFn();
          const long double prtFn21 = prtFn2 / prtFn1 * exp(-beta * HeatDiff);
          eqFrac += prtFn21;
        }
      }
      Reaction::sourceMap::iterator kpos;
      for (kpos = m_sources.begin(); kpos != m_sources.end(); ++kpos){
        SuperMolecule* pS2 = kpos->first;
        const long double HeatDiff = pS2->get_zpe() - pI1->get_zpe();
        const long double prtFn2 = pS2->rovibronicGrnCanPrtnFn();
        const long double trCon2 = translationalContribution((pS2->getMember1())->getMass(), (pS2->getMember2())->getMass(), beta);
        const long double excess2 = pS2->getExcessReactantConc();
        const long double prtFn21 = prtFn2 * trCon2 / (prtFn1 * excess2) * exp(-beta * HeatDiff);
        eqFrac += prtFn21;
      }
      eqFrac = 1.0 / eqFrac;
      pI1->setEqFraction(eqFrac);
      ctest << "Equilibrium Fraction for " << pI1->getName() << " = " << eqFrac << endl;
    }

    Reaction::sourceMap::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      SuperMolecule* pSx = spos->first;
      long double eqFrac(1.0);
      const long double prtFn1 = pSx->rovibronicGrnCanPrtnFn();
      const long double trCon1 = translationalContribution((pSx->getMember1())->getMass(), (pSx->getMember2())->getMass(), beta);
      const long double excess1 = pSx->getExcessReactantConc();
      Reaction::isomerMap::iterator jpos;
      for (jpos = m_isomers.begin(); jpos != m_isomers.end(); ++jpos){
        CollidingMolecule* pIy = jpos->first;
        const long double HeatDiff = pIy->get_zpe() - pSx->get_zpe();
        const long double prtFn2 = pIy->rovibronicGrnCanPrtnFn();
        const long double prtFn21 = prtFn2 * excess1/ (prtFn1 * trCon1) * exp(-beta * HeatDiff);
        eqFrac += prtFn21;
      }
      Reaction::sourceMap::iterator kpos;
      for (kpos = m_sources.begin(); kpos != m_sources.end(); ++kpos){
        if (kpos->first != pSx){
          SuperMolecule* pSy = kpos->first;
          const long double HeatDiff = pSy->get_zpe() - pSx->get_zpe();
          const long double prtFn2 = pSy->rovibronicGrnCanPrtnFn();
          const long double trCon2 = translationalContribution((pSy->getMember1())->getMass(), (pSy->getMember2())->getMass(), beta);
          const long double excess2 = pSy->getExcessReactantConc();
          eqFrac += prtFn2 * trCon2 * excess1 / (prtFn1 * trCon1 * excess2) * exp(-beta * HeatDiff);
        }
      }
      eqFrac = 1.0 / eqFrac;
      (pSx->getMember1())->setEqFraction(eqFrac);
      ctest << "Equilibrium Fraction for " << pSx->getName() << " = " << eqFrac << endl;
    }
    // description of the calculation: _2008_05_30__12_48_35_ on the end of the file.
    return true;
  }

  bool ReactionManager::produceInitialPopulationVector(vector<double>& eqFracCoeff, vector<double>& initDist){
    double populationSum = 0.0;
    Reaction::isomerMap::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      CollidingMolecule* isomer = ipos->first;
      populationSum += isomer->getInitPopulation();
    }
    Reaction::sourceMap::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      SuperMolecule* source = spos->first;
      populationSum += (source->getMember1())->getInitPopulation();
    }

    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      CollidingMolecule* isomer = ipos->first;
      double initFrac = isomer->getInitPopulation() / populationSum;
      int location = ipos->second;
      const double eqFrac = isomer->getEqFraction();
      const int colloptrsize = isomer->get_colloptrsize();
      vector<double> boltzFrac;
      isomer->normalizedGrainDistribution(boltzFrac, colloptrsize);
      for (int i = 0; i < colloptrsize; ++i){
        eqFracCoeff[i + location] = sqrt(boltzFrac[i] * eqFrac);
      }
      if (initFrac != 0.0){
        for (int i = 0; i < colloptrsize; ++i){
          initDist[i + location] = sqrt(initFrac * boltzFrac[i] / eqFrac);
        }
      }
    }

    // if there is no source term and the populationSum is still zero, set population = 1.0 for the first isomer
    int sizeSource = static_cast<int>(m_sources.size());
    if (populationSum == 0. && sizeSource == 0){
      ipos = m_isomers.begin();
      CollidingMolecule* isomer = ipos->first;
      isomer->setInitPopulation(1.0); // set initial population for the first isomer
      double initFrac = isomer->getInitPopulation();
      cinfo << "No population was assigned with no source term. Initialize the first isomer term to 1.0." << endl;
      int location = ipos->second;
      const double eqFrac = isomer->getEqFraction();
      const int colloptrsize = isomer->get_colloptrsize();
      vector<double> boltzFrac;
      isomer->normalizedGrainDistribution(boltzFrac, colloptrsize);
      for (int i = 0; i < colloptrsize; ++i){
        eqFracCoeff[i + location] = sqrt(boltzFrac[i] * eqFrac);
      }
      for (int i = 0; i < colloptrsize; ++i){
        initDist[i + location] = sqrt(initFrac * boltzFrac[i] / eqFrac);
      }
    }

    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      SuperMolecule* source = spos->first;
      double initFrac = (source->getMember1())->getInitPopulation() / populationSum;
      int location = spos->second;
      const double eqFrac = (source->getMember1())->getEqFraction();
      eqFracCoeff[location] = sqrt(eqFrac);
      initDist[location] = sqrt(initFrac) / eqFracCoeff[location];
      if (populationSum == 0. && spos == m_sources.begin()){
        cinfo << "No population was assigned. Initialize the first source term to 1.0." << endl;
        initDist[location] = 1.0 / eqFracCoeff[location];
      }
    }

    return true;
  }

  bool ReactionManager::timeEvolution(int maxTimeStep, const MesmerEnv mEnv)
  {
    int smsize = int(m_pSystemCollisionOperator->size());

    double maxEvoTime = 0.;
    if (mEnv.maxEvolutionTime <= 0. || mEnv.maxEvolutionTime > 1.0e8)
      maxEvoTime = 1.0e8;
    else
      maxEvoTime = mEnv.maxEvolutionTime;

    /* calculate the time points */
    vector<double> timePoints;
    for (int i = 0; i < maxTimeStep; ++i){
      double time = pow(10., static_cast<double>(i) / 10. - 11.);
      if (time > maxEvoTime)
        break;
      timePoints.push_back(time);
    }

    if (!calculateEquilibriumFractions(mEnv.beta)){
      cerr << "Failed calculating equilibrium fractions.";
      return false;
    }

    vector<double> initDist(smsize, 0.); // initial distribution
    vector<double> eqFracCoeff(smsize, 0.);   // equilibrium fraction coefficients
    if (!produceInitialPopulationVector(eqFracCoeff, initDist)){
      cerr << "Failed producing initial population vector.";
      return false;
    }

    // Coefficients due to the initial distribution
    dMatrix sysCollOptr(*m_pSystemCollisionOperator); // copy the system collision operator, which holds the eigenvectors
    vector<double> work1(smsize, 0.);
    for (int i = 0; i < smsize; ++i) {
      double sum = 0.;
      for (int j = 0; j < smsize; ++j) {
        sum += initDist[j] * sysCollOptr[j][i];
      }
      work1[i] = sum;
    }

    // Multiply equilibrium matrix with eigenvector matrix. 
    for (int i = 0; i < smsize; ++i) {
      double tmp = eqFracCoeff[i];
      for (int j = 0; j < smsize; ++j) {
        sysCollOptr[i][j] *= tmp;
      }
    }

    // populations calculated here
    db2D grainedProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
    vector<double> work2(smsize, 0.);
    for (int timestep = 0; timestep < maxTimeStep; ++timestep){
      double numColl = m_meanOmega * timePoints[timestep];
      for (int j = 0; j < smsize; ++j) {
        work2[j] = work1[j] * exp(m_eigenvalues[j] * numColl);
      }
      for (int j = 0; j < smsize; ++j) {
        double sum = 0.;
        for (int l = 0; l < smsize; ++l) {
          sum += work2[l] * sysCollOptr[j][l];
        }
        grainedProfile[j][timestep] = sum;
      }
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
    // Sum up individual species
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
      ctest << (jpos->second)->getName() << " ";
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        for (int i = 0; i < colloptrsize; ++i){
          speciesProfile[j][timestep] += grainedProfile[i + idx][timestep];
        }
        formatFloat(ctest, speciesProfile[j][timestep], 6,  15);
      }
      ctest << endl;
    }

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
        m_populations[sRef] = population;
      }
      ppInitMol = ppInitMol->XmlMoveTo("molecule");
    }
  }

  // _2008_05_30__12_48_35_
  //
  // For a system with one source term and two isomer terms as the following
  // excess + A ---> B
  // B <--> C
  // The equilibrium population fraction of A will be expressed as
  //                         1.0
  // PA = ------------------------------------------
  //            Q_rve_B [excess]    Q_rve_C [excess]
  //      1.0 + ----------------  + ----------------
  //            Q_rvet_Reactants    Q_rvet_Reactants
  //
  // where in the expression, Q_rve is rovibronic contribution of the partition function,
  // and t is the translational contribution of the partition function.
  // [excess] is the number density of the excess reactant.
  //
  // Also, the equilibrium population fraction of B will be
  //                      1.0
  // PB = ----------------------------------
  //            Q_rvet_Reactants     Q_rve_C
  //      1.0 + ----------------  +  -------
  //            Q_rve_B [excess]     Q_rve_B

}//namespace
