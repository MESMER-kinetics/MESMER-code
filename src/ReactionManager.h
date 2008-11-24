#ifndef GUARD_ReactionManager_h
#define GUARD_ReactionManager_h

//-------------------------------------------------------------------------------------------
//
// ReactionManager.h
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This header file contains the declaration of the ReactionManager class.
// This class will contain the reactions that go to make up a system.
//
//-------------------------------------------------------------------------------------------

#include "AssociationReaction.h"
#include "IrreversibleUnimolecularReaction.h"
#include "IsomerizationReaction.h"
#include "IrreversibleExchangeReaction.h"

namespace mesmer
{
  class ReactionManager
  {
  public:

    // Type defs
    typedef  size_t  size_type ;
    typedef std::map<std::string , double> populationMap ;
    typedef std::map<Reaction* , int, Reaction::ReactionPtrLess> sinkMap ;

    ReactionManager(MoleculeManager *pMoleculeManager);

    // Destructor.
    ~ReactionManager(){} ;

    // Add a new reaction to the map.
    bool addreactions(PersistPtr ReacList, const MesmerEnv& mEnv, MesmerFlags& mFlags) ;

    // Remove a reaction from the map.
    void remove(){} ;

    void resetCalcFlags();

    // Total number of reaction in map.
    size_type size() const {return m_reactions.size() ; } ;

    // Find a particular reaction.
    Reaction*       operator[](const size_type i)       { return m_reactions[i] ; } ;
    const Reaction* operator[](const size_type i) const { return m_reactions[i] ; } ;

    // Find a reaction from its id
    Reaction* find(const std::string& id) const ;

    // Build collision operator for system.
    bool BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags) ;

    // Diagonalize the collision operator.
    void diagCollisionOperator(const MesmerFlags &mFlags, const int precision) ;

    // Calculate the time evolution of the system
    bool timeEvolution(const MesmerFlags mFlags, const int precision, dMatrix& rates);

    // Set Initial population for individual species
    void setInitialPopulation(PersistPtr);

    bool calculateEquilibriumFractions(const double beta);

    double calcChiSquare(const dMatrix& mesmerRates, vector<conditionSet>& expRates);

  private:

    std::vector<Reaction *> m_reactions ;

    MoleculeManager        *m_pMoleculeManager ;

    qdMatrix               *m_pReactionOperator ;

    std::vector<qd_real>    m_eigenvalues;

    std::vector<qd_real>     m_eqVector;

    // Maps the location of individual reactant collision operator and source terms in the system matrix.
    Reaction::isomerMap    m_isomers;
    Reaction::sourceMap    m_sources;
    sinkMap                m_sinkRxns;
    populationMap          m_initialPopulations;

    // map modelled molecules (isomers + sources) with their sequence in the EqMatrix and Rate Coefficient matrix
    Reaction::sourceMap    m_SpeciesSequence;

    sinkMap m_SinkSequence;

    double m_meanOmega;

    // Extract molecule information from XML stream.
    bool GetMoleculeInfo(PersistPtr pp, std::string& MolName, std::string& MolType) ;

    // sets grain parameters and determines system environment
    bool SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne);

    bool produceInitialPopulationVector(vector<qd_real>& initDist);

    bool produceEquilibriumVector();

    void printCollisionOperator(const MesmerFlags &mFlags);

    template <class T>
    void constructSpeciesProfile(TMatrix<T> sysCollOptr, vector<T> eigenValues, const vector<T>& n_0, const vector<T> eqvector, const MesmerFlags mFlags){

      T maxEvoTime = 0.;
      // set the default maximum evolution time
      if (mFlags.maxEvolutionTime <= 0. || mFlags.maxEvolutionTime > 1.0e8)
        maxEvoTime = 1.0e8;
      else
        maxEvoTime = mFlags.maxEvolutionTime;

      // calculate the time points
      vector<T> timePoints;
      for (int i = 0; i <= 160; ++i){
        T time = pow(10., static_cast<T>(i) / 10. - 11.);
        if (time > maxEvoTime)
          break;
        timePoints.push_back(time);
      }

      //initialize dt vector for calculating product yields
      vector<T> dt(timePoints.size()-1,0.0);
      dt[0] = timePoints[0];
      for (int i = 1; i < int(dt.size()); ++i){
        dt[i] = timePoints[i] - timePoints[i-1];
      }

      const int smsize = int(sysCollOptr.size());
      vector<T> r_0(smsize, 0.);            // which holds the eigenvectors, |r_0> = U^-1 |n_0>

      for (int i = 0; i < smsize; ++i) {
        T sum = 0.;
        for (int j = 0; j < smsize; ++j) {
          sum += n_0[j] * sysCollOptr[j][i];
        }
        r_0[i] = sum;  // now |r_0> = V^(T)*|init> = U^(-1)*|n_0>
      }

      for (int i = 0; i < smsize; ++i) {
        T tmp = eqvector[i];
        for (int j = 0; j < smsize; ++j) {
          sysCollOptr[i][j] *= tmp;
        }
      }

      const int maxTimeStep = int(dt.size());
      a2d_t<T> grnProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
      vector<T> work2(smsize, 0.);

      for (int timestep = 0; timestep < maxTimeStep; ++timestep){
        T meanOmegaTemp = m_meanOmega;
        T numColl = meanOmegaTemp * timePoints[timestep];
        for (int j = 0; j < smsize; ++j) {
          work2[j] = r_0[j] * exp(eigenValues[j] * numColl);
        } // now |wk2> = exp(Dt)*V^(T)*|init> = exp(Dt)*U^(-1)*|n_0>
        vector<T> sumVector;
        for (int j = 0; j < smsize; ++j) {
          T sum(0.0);
          for (int l = 0; l < smsize; ++l) {
            sum += work2[l] * sysCollOptr[j][l];
          }
          sumVector.push_back(sum);
        } // now |grnProfile(t)> = |grnProfile(i)> = F*V*exp(Dt)*V^(T)*|init> = U*exp(Dt)*U^(-1)*|n_0>
        for (int j = 0; j < smsize; ++j) {
          grnProfile[j][timestep] = sumVector[j];
        }
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

      vector<T> totalIsomerPop(maxTimeStep, 0.);
      vector<T> totalPdtPop(maxTimeStep, 0.);

      for(int timestep(0); timestep<maxTimeStep; ++timestep){
        for(int j(0);j<smsize;++j){
          T someGrain = grnProfile[j][timestep];
          totalIsomerPop[timestep] += someGrain;
        }
        T popTime = totalIsomerPop[timestep];
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
        int sinkpos(0);
        m_sinkRxns.clear();                      // Populate map: m_sinkRxns[Rxns] = location of rct
        for (size_t i(0) ; i < size() ; ++i) {
          bool Irreversible = (m_reactions[i])->isIrreversible() ;
          if (Irreversible && m_sinkRxns.find(m_reactions[i]) == m_sinkRxns.end()) {   // add an irreversible rxn to the map
            Reaction* reaction = m_reactions[i];
            ModelledMolecule* source = reaction->get_reactant();
            CollidingMolecule* isomer = dynamic_cast<CollidingMolecule*>(source);
            if(isomer){
              int location = m_isomers[isomer];
              m_sinkRxns[reaction] = location;
              m_SinkSequence[reaction] = sinkpos;               // populate SinkSequence map with Irreversible Rxns
              ++sinkpos;
            }
            else if(source){
              int location = m_sources[source];
              m_sinkRxns[reaction] = location;
              m_SinkSequence[reaction] = sinkpos;
              ++sinkpos;
            }
          }
        }

        int numberOfSpecies = static_cast<int>(m_isomers.size() + m_sources.size() + m_sinkRxns.size());
        a2d_t<T> speciesProfile(numberOfSpecies, maxTimeStep);
        int speciesProfileidx(0);

        ctest << setw(16) << "Timestep/s ";

        Reaction::sourceMap::iterator spos;
        for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map
          ModelledMolecule* source = spos->first ;                        // to get source profile vs t
          ctest << setw(16) << source->getName();
          int location = spos->second;
          for (int timestep = 0; timestep < maxTimeStep; ++timestep){
            T gPf = grnProfile[location][timestep];
            speciesProfile[speciesProfileidx][timestep] = gPf;
          }
          ++speciesProfileidx;
        }

        Reaction::isomerMap::iterator ipos;
        for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
          CollidingMolecule* isomer = ipos->first;                        // to get isomer profile vs t
          ctest << setw(16) << isomer->getName();
          int location = ipos->second;
          const int colloptrsize = isomer->get_colloptrsize();
          for (int timestep = 0; timestep < maxTimeStep; ++timestep){
            T sumGP(0.0);
            for(int i = 0; i < colloptrsize; ++i){
              T gpop = grnProfile[i+location][timestep];
              sumGP += gpop;
            }
            speciesProfile[speciesProfileidx][timestep] = sumGP;
          }
          ++speciesProfileidx;
        }

        sinkMap::iterator pos;      // iterate through sink map to get product profile vs t
        int pdtProfileStartIdx = speciesProfileidx;
        for (pos = m_sinkRxns.begin(); pos != m_sinkRxns.end(); ++pos){
          vector<double> KofEs;                             // vector to hold sink k(E)s
          Reaction* sinkReaction = pos->first;
          const int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
          vector<ModelledMolecule*> pdts;                               // in the sink reaction
          sinkReaction->get_products(pdts);
          if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
            KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
            ctest << setw(11) << pdts[0]->getName()<< setw(5) << "(bim)";
          }
          else{   // if the collision operator size is >1, there are k(E)s for the irreversible loss
            KofEs = sinkReaction->get_GrainKfmc();          // assign sink k(E)s, the vector size == maxgrn
            ctest << setw(16) << pdts[0]->getName();
          }
          int location = pos->second;                       // get sink location
          T TimeIntegratedProductPop(0.0);
          for (int timestep = 0; timestep < maxTimeStep; ++timestep){
            for(int i = 0; i < colloptrsize; ++i){
              T tempK = KofEs[i];
              speciesProfile[speciesProfileidx][timestep] += tempK * grnProfile[i+location][timestep] * dt[timestep];
            }
            TimeIntegratedProductPop += speciesProfile[speciesProfileidx][timestep];
            speciesProfile[speciesProfileidx][timestep]= TimeIntegratedProductPop;
          }
          ++speciesProfileidx;
          KofEs.clear();
        }

        if (pdtProfileStartIdx < speciesProfileidx){
          for(int timestep = 0; timestep < maxTimeStep; ++timestep){    // normalize product profile to account for small
            T normConst(0.0);                          // numerical errors in TimeIntegratedProductPop
            T pdtYield(0.0);
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
    }

    template <class T>
    bool BartisWidomPhenomenologicalRates(TMatrix<T>& mesmerRates, TMatrix<T>& TMev, vector<T>& TVev, vector<T> eqvec)
    {
      const int smsize = int(TMev.size());
      TMatrix<T> eigenVec(smsize);  //copy ReactionOperator, the eigenvector Matrix (== V)

      for ( int i = 0 ; i < smsize ; ++i )
        for ( int j = 0 ; j < smsize ; ++j )
          eigenVec[i][j] = TMev[i][j] ;

      // constant variables
      const int nchem = static_cast<int>(m_isomers.size() + m_sources.size());  // number of isomers+pseudoisomers
      const int nchemIdx = smsize - nchem;       // idx for chemically significant eigenvalues & vectors

      TMatrix<T> assymInvEigenVec(smsize);   // U^(-1)
      TMatrix<T> assymEigenVec(smsize);      // U
      for(int i(0);i<smsize;++i){
        T tmp = eqvec[i];
        for(int j(0);j<smsize;++j){
          assymInvEigenVec[j][i] = eigenVec[i][j]/tmp;         //calculation of U^(-1) = (FV)^-1 = V^T * F^-1
          assymEigenVec[j][i] = eqvec[j] * eigenVec[j][i];  //calculation of U = FV
        }
      }

      //------------------------- TEST block ----------------------------------------
      TMatrix<T> EigenVecIdentity(smsize);   // matrix for holding product of U^(-1) * U
      for(int i(0);i<smsize;++i){         // multiply U*U^(-1) for testing
        T test = 0.0;
        for(int j(0);j<smsize;++j){
          T sm = 0.0;
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

      TMatrix<T> Z_matrix(nchem);  // definitions of Y_matrix and Z_matrix taken from PCCP 2007(9), p.4085
      a2d_t<T> Y_matrix;
      Reaction::isomerMap::iterator ipos;  // set up an iterator through the isomer map
      Reaction::sourceMap::iterator spos;  // set up an iterator through the source map
      sinkMap::iterator sinkpos;           // set up an iterator through the irreversible rxn map

      ctest << "\nBartis Widom eigenvalue/eigenvector analysis\n";
      ctest << endl << "Number of sinks in this system: " << m_sinkRxns.size() << endl;

      for(int i(0); i<nchem; ++i){
        for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // calculate Z_matrix matrix elements for
          T sm = 0.0;                                                // all isomers in the system
          CollidingMolecule* isomer = ipos->first;
          const int colloptrsize = isomer->get_colloptrsize();            // get colloptrsize for isomer
          int location = ipos->second;                                    // get location for isomer
          int position = m_SpeciesSequence[isomer];                       // get sequence position for isomer
          for(int j(0);j<colloptrsize;++j){
            sm += assymEigenVec[location+j][nchemIdx+i];
          }
          Z_matrix[position][i] = sm;
        }
        for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // calculate Z_matrix matrix elements for
          T sm = 0.0;                                                // all sources in the system
          ModelledMolecule* pPseudoIsomer = spos->first ;
          int location = spos->second;
          int position = m_SpeciesSequence[pPseudoIsomer];
          sm = assymEigenVec[location][nchemIdx+i];
          Z_matrix[position][i] = sm;
        }
        if(m_sinkRxns.size()!=0) {
          for(sinkpos=m_sinkRxns.begin(); sinkpos!=m_sinkRxns.end(); ++sinkpos){ // calculate Y_matrix elements for sinks
            T sm = 0.0;
            vector<double> KofEs;                                         // vector to hold sink k(E)s
            Reaction* sinkReaction = sinkpos->first;
            const int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
            if(colloptrsize == 1)  // if the collision operator size is 1, there is one canonical loss rate coefficient
              KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
            else                   // if the collision operator size is >1, there are k(E)s for the irreversible loss
              KofEs = sinkReaction->get_GrainKfmc();                      // assign sink k(E)s, the vector size == maxgrn
            int location = sinkpos->second;                               // get sink location
            int position = m_SinkSequence[sinkReaction];                  // get sink sequence position
            for(int j(0);j<colloptrsize;++j){
              sm += assymEigenVec[location+j][nchemIdx+i] * KofEs[j];
            }
            Y_matrix[position][i] = sm;
            KofEs.clear();
          }
        }
      }

      //    Y_matrix.print((int)(m_sinkRxns.size()), (int)(m_SpeciesSequence.size())); // print out Y_matrix for testing

      TMatrix<T> Zinv(Z_matrix), Zidentity(nchem), Kr(nchem);
      a2d_t<T> Kp;

      if(Zinv.invertAdjointCofactors()){
        cerr << "Inversion of Z_matrix failed.";
      }
      ctest << "\nZ_matrix: ";
      Z_matrix.showFinalBits(nchem, true);

      ctest << endl << "Z_matrix^(-1):" << endl;
      Zinv.showFinalBits(nchem, true);

      for(int i(0);i<nchem;++i){          // multiply Z_matrix*Z_matrix^(-1) for testing
        for(int j(0);j<nchem;++j){
          T sm = 0.0;
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
          T sm = 0.0;
          for(int k(0);k<nchem;++k){
            sm += Z_matrix[i][k] * TVev[nchemIdx+k] * Zinv[k][j];
          }
          T meanOmegaTemp = m_meanOmega;
          Kr[i][j] = sm * meanOmegaTemp;
        }
      }
      ctest << "\nKr matrix:" << endl;
      Kr.showFinalBits(nchem, true);       // print out Kr_matrix

      if(m_sinkRxns.size()!=0){
        for(int i(0); i != int(m_sinkRxns.size()); ++i){    // calculate Kp (definition taken from PCCP 2007(9), p.4085)
          for(int j(0);j<nchem;++j){
            T sm = 0.0;
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
      Reaction::sourceMap::iterator lossitr, rctitr, pdtitr;

      // print pseudo 1st order k loss for isomers
      for(lossitr=m_SpeciesSequence.begin(); lossitr!=m_SpeciesSequence.end(); ++lossitr){
        ModelledMolecule* iso = lossitr->first;
        int losspos = lossitr->second;
        ctest << iso->getName() << " loss = " << Kr[losspos][losspos] << endl;
      }

      ctest << "}\n";
      ctest << "\nFirst order & pseudo first order rate coefficients for isomerization rxns:\n{\n";

      // print pseudo first order connecting ks
      for (rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
        ModelledMolecule* rct = rctitr->first;
        int rctpos = rctitr->second;
        for (pdtitr=m_SpeciesSequence.begin(); pdtitr!=m_SpeciesSequence.end(); ++pdtitr){
          ModelledMolecule* pdt = pdtitr->first;
          int pdtpos = pdtitr->second;
          if(rctpos != pdtpos)
            ctest << rct->getName() << " -> " << pdt->getName() << " = " << Kr[pdtpos][rctpos] << endl;
        }
      }
      ctest << "}\n";

      if(m_sinkRxns.size()!=0){
        ctest << "\nFirst order & pseudo first order rate coefficients for irreversible rxns:\n{\n";
        sinkMap::iterator sinkitr;

        for(sinkitr=m_SinkSequence.begin(); sinkitr!=m_SinkSequence.end(); ++sinkitr){
          Reaction* sinkReaction = sinkitr->first;          // get Irreversible Rxn
          int colloptrsize = sinkReaction->getRctColloptrsize();
          int sinkpos = m_SinkSequence[sinkReaction];                   // get products & their position
          vector<ModelledMolecule*> pdts;
          sinkReaction->get_products(pdts);
          for(rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
            ModelledMolecule* rcts = rctitr->first;     // get reactants & their position
            int rctpos = rctitr->second;
            if(colloptrsize==1)
              ctest << rcts->getName() << " -> "  << pdts[0]->getName() << "(bim) = " << Kp[sinkpos][rctpos] << endl;
            else
              ctest << rcts->getName() << " -> "  << pdts[0]->getName() << " = " << Kp[sinkpos][rctpos] << endl;
          }
        }
        ctest << "}\n\n";
      }

      mesmerRates = Kr;
      return true;
    }

  } ;

}//namespace

#endif // GUARD_ReactionManager_h
