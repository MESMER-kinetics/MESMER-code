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
      // Note: 1. A source term is probably the only deficient reactant that initiates
      //          the whole process of reactions in the master equation. In this case
      //          we think there may be more than one source terms.
      //       2. In the current construction of Mesmer, the source is a SuperMolecule
      //          representing both reactants.
      m_sources.clear(); // Maps the location of source in the system matrix.
      for (size_t i(0) ; i < size() ; ++i) {
        AssociationReaction *pReaction = dynamic_cast<AssociationReaction*>(m_reactions[i]) ;
        if (pReaction) {
          SuperMolecule *pSuperMolecule = pReaction->get_bi_molecularspecies();
          if (pSuperMolecule && m_sources.find(pSuperMolecule) == m_sources.end()){ // New source
            m_sources[pSuperMolecule] = msize ;
            pReaction->putSourceMap(&m_sources) ;
            ++msize ;
          }
        }
      }

      // Allocate space for system collision operator.
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
      CollidingMolecule* pI0 = ipos->first;
      long double eqFrac(1.0);
      const long double prtFn1 = pI0->rovibronicGrnCanPrtnFn();
      Reaction::isomerMap::iterator jpos;
      for (jpos = m_isomers.begin(); jpos != m_isomers.end(); ++jpos){
        if (jpos->first != pI0){
          CollidingMolecule* pI1 = jpos->first;
          const long double HeatDiff = pI1->get_zpe() - pI0->get_zpe();
          const long double prtFn2 = pI1->rovibronicGrnCanPrtnFn();
          eqFrac += prtFn2 / prtFn1 * exp(-beta * HeatDiff);
        }
      }
      Reaction::sourceMap::iterator kpos;
      for (kpos = m_sources.begin(); kpos != m_sources.end(); ++kpos){
        SuperMolecule* pS1 = kpos->first;
        const long double HeatDiff = pS1->get_zpe() - pI0->get_zpe();
        long double prtFn2 = pS1->rovibronicGrnCanPrtnFn();
        prtFn2 *= translationalContribution((pS1->getMember1())->getMass(), (pS1->getMember2())->getMass(), beta);
        eqFrac += prtFn2 / prtFn1 * exp(-beta * HeatDiff);
      }
      pI0->setEqFraction(1.0 / eqFrac);
    }

    Reaction::sourceMap::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      SuperMolecule* pSx = spos->first;
      long double eqFrac(1.0);
      long double prtFn1 = pSx->rovibronicGrnCanPrtnFn();
      prtFn1 *= translationalContribution((pSx->getMember1())->getMass(), (pSx->getMember2())->getMass(), beta);
      Reaction::isomerMap::iterator jpos;
      for (jpos = m_isomers.begin(); jpos != m_isomers.end(); ++jpos){
        CollidingMolecule* pSy = jpos->first;
        const long double HeatDiff = pSy->get_zpe() - pSx->get_zpe();
        const long double prtFn2 = pSy->rovibronicGrnCanPrtnFn();
        eqFrac += prtFn2 / prtFn1 * exp(-beta * HeatDiff);
      }
      Reaction::sourceMap::iterator kpos;
      for (kpos = m_sources.begin(); kpos != m_sources.end(); ++kpos){
        if (kpos->first != pSx){
          SuperMolecule* pSy = kpos->first;
          const long double HeatDiff = pSy->get_zpe() - pSx->get_zpe();
          long double prtFn2 = pSy->rovibronicGrnCanPrtnFn();
          prtFn2 *= translationalContribution((pSy->getMember1())->getMass(), (pSy->getMember2())->getMass(), beta);
          eqFrac += prtFn2 / prtFn1 * exp(-beta * HeatDiff);
        }
      }
      (pSx->getMember1())->setEqFraction(1.0 / eqFrac);
    }
    return true;
  }

  bool ReactionManager::produceInitialPopulationVector(vector<double>& eqFracCoeff, vector<double>& initDist){
    Reaction::isomerMap::iterator ipos;
    for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
      CollidingMolecule* isomer = ipos->first;
      int location = ipos->second;
      const double eqFrac = isomer->getEqFraction();
      const int colloptrsize = isomer->get_colloptrsize();
      vector<double> boltzFrac;
      isomer->normalizedGrainDistribution(boltzFrac, colloptrsize);
      for (int i = 0; i < colloptrsize; ++i){
        eqFracCoeff[i + location] = sqrt(boltzFrac[i] * eqFrac);
      }
    }
    Reaction::sourceMap::iterator spos;
    for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
      SuperMolecule* source = spos->first;
      int location = spos->second;
      const double eqFrac = (source->getMember1())->getEqFraction();
      eqFracCoeff[location] = sqrt(eqFrac);
      initDist[location] = 1.0 / eqFracCoeff[location];
    }

    return true;
  }

  bool ReactionManager::timeEvolution(int maxTimeStep, const double beta)
  {
    int smsize = int(m_pSystemCollisionOperator->size());

    /* calculate the time points */
    vector<double> timePoints;
    for (int i = 0; i < maxTimeStep; ++i){
      timePoints.push_back(pow(10., static_cast<double>(i) / 10. - 11.));
    }

    if (!calculateEquilibriumFractions(beta)){
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
    dMatrix sysCollOptr(*m_pSystemCollisionOperator); // copy the system collision operator
    vector<double> work1(smsize, 0.);
    for (int i = 0; i < smsize; ++i) {
      double sum = 0.;
      for (int j = 0; j < smsize; ++j) {
        sum += initDist[j] * sysCollOptr[j][i];
      }
      work1[i] = sum;
    }

    // Multiply root matrix with eigenvector matrix. Note VT contains transpose of V
    for (int i = 0; i < smsize; ++i) {
      double tmp = eqFracCoeff[i];
      for (int j = 0; j < smsize; ++j) {
        sysCollOptr[i][j] *= tmp;
      }
    }

    // populations calculated here
    db2D speciesProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
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
        speciesProfile[j][timestep] = sum;
      }
    }

    // print species profile
    ctest << "\nSpecies profile (the first row is time points in unit of second):\n{\n";
    for (int timestep = 0; timestep < maxTimeStep; ++timestep){
      formatFloat(ctest, timePoints[timestep], 6,  15);
    }
    ctest << endl;
    for (int j = 0; j < smsize; ++j) {
      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        formatFloat(ctest, speciesProfile[j][timestep], 6,  15);
      }
      ctest << endl;
    }
    ctest << "}\n";
    return true;
  }
}//namespace
