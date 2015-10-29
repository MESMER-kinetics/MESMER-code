//-------------------------------------------------------------------------------------------
// gWellProperties.h
//
// Author: Chi-Hsiu Liang
//
// This file contains property groups of class Molecule. These groups give molecules variables
// and functions to perform tasks; from the definitions of these groups a molecule can play
// roles when it is required to do so. Classes in this file do not depend on each other and
// thus they can be seperated. Any of them can be added into a molecule (with a new() to construct
// an object and then pass the pointer to the molecule) when the role of the molecule requires
// the information in that group.
//-------------------------------------------------------------------------------------------

#ifndef GUARD_gWellProperties_h
#define GUARD_gWellProperties_h

#include <vector> 
#include "MolecularComponents.h"
#include "gDensityOfStates.h"

namespace mesmer
{

  class gWellProperties:public MolecularComponent
  {
    //-------------------------------------------------------------------------------------------------
    // Collisional transfer related properties
    //-------------------------------------------------------------------------------------------------

  private:

    double m_collisionFrequency ; // Current value of collision frequency.
    size_t m_ncolloptrsize ;      // Size of the collision operator matrix.
    double m_lowestBarrier;       // lowest barrier associatied with this species
    size_t m_numGroupedGrains;    // Number of grains grouped into a reservoir grain.

    DistributionCalculator* m_pDistributionCalculator;

    std::map<std::string, EnergyTransferModel*> m_EnergyTransferModels; //with different bath gases

    std::vector<double>  m_grainDist ; // Grain distribution (not normalized)
    qdMatrix            *m_egvec ;     // Eigenvectors used to diagonalize (P - I) matrix.
    std::vector<qd_real> m_egval;

    // Calculate collision frequency.
    double collisionFrequency(MesmerEnv env, Molecule *pBathGasMolecule) ;

    // Calculate raw transition matrix.
    template<class T> 
    bool rawTransitionMatrix(MesmerEnv& env, vector<double> &gEne,  vector<double> &gDOS, TMatrix<T>* egme) ;

    // Write out collision operator diaganostics.
    template<class T> 
    void writeCollOpProps(vector<double>& ene, TMatrix<T>* egme) const;

    // Construct reservoir state.
    template<class T> 
    void constructReservoir(MesmerEnv& env, vector<double> &gEne, vector<double> &gDOS, TMatrix<T>* egme) ;

  public:

    //
    // Constructor, destructor and initialization
    //
    gWellProperties(Molecule* pMol);
    virtual ~gWellProperties();
    bool initialization();

    // Returns an existing model associated with the named bath gas or makes a new one
    EnergyTransferModel* addBathGas(const char* pbathGasName, EnergyTransferModel* pModel);

    // Initialize the Collision Operator.
    bool initCollisionOperator(MesmerEnv& env, Molecule *pBathGasMolecule) ;

    // Calculate a reaction matrix element.
    qd_real matrixElement(int eigveci, int eigvecj, std::vector<double> &k) const;

    // Accessor a collision operator eigenvector.
    void eigenVector(int eigveci, std::vector<double> &evec) const ;

    // Calculate collision operator.
    template<class T> 
    bool collisionOperator (MesmerEnv& env, TMatrix<T> **egme) ;

    // Diagonalize collision operator before adding reaction terms to get eigenvectors and eigenvalues.
    void diagonalizeCollisionOperator(qdMatrix *egme);

    template<class T> 
    void copyCollisionOperator(qdMatrix *CollOptr, TMatrix<T> *egme, const size_t locate, const double RducdOmega) const ;

    void copyCollisionOperatorEigenValues(qdMatrix *CollOptr, const size_t locate, const double RducdOmega) const ;

    void normalizedInitialDistribution(vector<double> &grainFrac) ;
    void normalizedGrnBoltzmannDistribution(vector<double> &grainFrac);
    void normalizedCellBoltzmannDistribution(vector<double> &grainFrac);

    // Accessors.

    double get_collisionFrequency() const {return m_collisionFrequency ; } ;

    void set_colloptrsize(int ncolloptrsize) { m_ncolloptrsize = ncolloptrsize; };

    size_t get_colloptrsize() const {return m_ncolloptrsize ; } ;

    size_t get_nbasis() const ;

    const int get_grnZPE();

    const double getLowestBarrier() { return m_lowestBarrier;}

    void setLowestBarrier(double value){ m_lowestBarrier = value;}

    const size_t reservoirShift() {return m_numGroupedGrains == 0 ? 0 : m_numGroupedGrains - 1; }

  };

  //
  // Calculate collision operator
  //
  template<class T> 
  bool gWellProperties::collisionOperator(MesmerEnv& env, TMatrix<T> **CollOp)
  {
    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    const size_t reducedCollOptrSize = m_ncolloptrsize - reservoirShift();

    // Allocate memory.
    TMatrix<T>* egme = new TMatrix<T>(m_ncolloptrsize);

    // Calculate raw transition matrix.
    if (!rawTransitionMatrix(env, gEne, gDOS, egme)) return false;

    if (m_host->getFlags().showCollisionOperator != 0){
      ctest << "\nCollision operator of " << m_host->getName() << " before normalization:\n";
      egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    //Normalisation
    egme->normalizeProbabilityMatrix();

    if (m_host->getFlags().showCollisionOperator >= 1){
      ctest << "\nCollision operator of " << m_host->getName() << " after normalization:\n";
      egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    // If requested write out column sums etc. to check normalization results.
    if (m_host->getFlags().reactionOCSEnabled){
      writeCollOpProps(gEne, egme) ;
    }

    // Construct reservoir state if specifed.
    if (m_numGroupedGrains > 1) {
      constructReservoir(env, gEne, gDOS, egme) ;
    }

    vector<T> popDist; // Grained population distribution.
    popDist.push_back(0.0);
    T prtnFn(0.0), beta(env.beta) ;
    for (size_t idx(0); idx < m_ncolloptrsize; ++idx) {
      const T tmp(T(gDOS[idx])*exp(-beta*T(gEne[idx])));
      prtnFn += tmp ;
      if (idx < std::max(m_numGroupedGrains,size_t(1))){
        popDist[0] += tmp;
      }
      else {
        popDist.push_back(tmp);
      }
    }

    // Symmetrization of the collision matrix.
    for (size_t i(1); i < reducedCollOptrSize; ++i) {
      for (size_t j(0); j < i; ++j){
        (*egme)[j][i] *= sqrt(popDist[i] / popDist[j]) ;
        (*egme)[i][j] = (*egme)[j][i];
      }
    }

    // Account for collisional loss by subrtacting unity from the leading diagonal.
    // SHR: note the slightly complex lower limit below improves accuracy at lower 
    // temperatures where reservoir states are used.
    for (size_t i((m_numGroupedGrains > 1) ? 1 : 0); i < reducedCollOptrSize; ++i) {
      (*egme)[i][i] -= 1.0;
    }

    if (m_host->getFlags().showCollisionOperator >= 2){
      ctest << "Collision operator of " << m_host->getName() << " after :\n";
      egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    m_ncolloptrsize = reducedCollOptrSize;

    if (*CollOp) 
      delete *CollOp;  // Delete any existing matrix.
    (*CollOp) = new TMatrix<T>(reducedCollOptrSize);
    (**CollOp) = (*egme) ;

    delete egme;

    return true;
  }

  //
  // Calculate raw transition matrix.
  //
  template<class T> 
  bool gWellProperties::rawTransitionMatrix(MesmerEnv& env, vector<double> &gEne, vector<double> &gDOS, TMatrix<T>* egme)
  {
    EnergyTransferModel* pEnergyTransferModel = m_EnergyTransferModels[env.bathGasName];
    if (!pEnergyTransferModel)
    {
      cerr << "No energyTransferModel for " << m_host->getName() << " with " << env.bathGasName << endl;
      return false;
    }

    T beta = T(env.beta) ;

    // Use number of states to weight the downward transition
    if (m_host->getFlags().useDOSweightedDT){
      // The collision operator.
      for (size_t i = 0; i < m_ncolloptrsize; ++i) {
        T ei = T(gEne[i]);
        T ni = T(gDOS[i]);
        for (size_t j = i; j < m_ncolloptrsize; ++j) {
          T ej = T(gEne[j]);
          T nj = T(gDOS[j]);
          // Transfer to lower Energy -
          // double transferDown = exp(-alpha*(ej - ei)) * (ni/nj);
          // (*m_egme)[i][j] = transferDown;
          T transferDown = T(pEnergyTransferModel->calculateTransitionProbability(to_double(ej), to_double(ei)));
          (*egme)[i][j] = transferDown * (ni / nj);

          // Transfer to higher Energy (via detailed balance) -
          // double transferUp = exp(-(alpha + beta)*(ej - ei));
          // (*m_egme)[j][i] = transferUp;
          (*egme)[j][i] = transferDown * exp(-beta*(ej - ei));
        }
      }
    }
    else {
      // The collision operator.
      for (size_t i = 0; i < m_ncolloptrsize; ++i) {
        T ei = T(gEne[i]);
        T ni = T(gDOS[i]);
        for (size_t j = i; j < m_ncolloptrsize; ++j) {
          T ej = T(gEne[j]);
          T nj = T(gDOS[j]);
          // Transfer to lower Energy -
          T transferDown = T(pEnergyTransferModel->calculateTransitionProbability(to_double(ej), to_double(ei))) ;
          (*egme)[i][j] = transferDown;

          // Transfer to higher Energy (via detailed balance) -
          (*egme)[j][i] = transferDown * (nj / ni) * exp(-beta*(ej - ei));
        }
      }
    }

    return true;
  }

  //
  // Write out collision operator diaganostics.
  //
  template<class T> 
  void gWellProperties::writeCollOpProps(vector<double>& ene, TMatrix<T>* egme) const {
    ctest << endl << "Collision operator column sums and energy transfer parameters" << endl << "{" << endl;
    ctest << " Column Sums           E   <Delta E>  <Delta E>d  <Delta E>u" << endl;
    for (size_t i(0); i < m_ncolloptrsize; ++i) {
      T columnSum(0.0);
      T meanEnergyTransfer(0.0);
      T meanEnergyTransferDown(0.0);
      T meanEnergyTransferUp(0.0);
      for (size_t j(0); j < m_ncolloptrsize; ++j){
        T trnsPrb = (*egme)[j][i] ;
        T eneMom  = (ene[j] - ene[i])*trnsPrb ;
        columnSum += trnsPrb ;
        meanEnergyTransfer += eneMom ;
        if (ene[j] < ene[i]) {
          meanEnergyTransferDown += eneMom ;
        }
        else {
          meanEnergyTransferUp += eneMom ;
        }
      }
      ctest << formatFloat(columnSum, 3, 12)
        << formatFloat(ene[i], 3, 12)
        << formatFloat(meanEnergyTransfer, 3, 12)
        << formatFloat(meanEnergyTransferDown, 3, 12)
        << formatFloat(meanEnergyTransferUp, 3, 12)
        << endl;
    }
    ctest << "}" << endl;
  }


  //
  // Construct reservoir state. This method calculates the total transition probability 
  // into the reservoir and then uses detailed balance to construct the probability of
  // excitiation from the the reservoir: 
  // k_a * x_r = k_d(E) * f(E) / Q_a * x_a
  // where Q_a is equal to x_a and cancelled out.
  // So, k_a = k_d(E) * f(E) / x_r;
  // Communication between regular grains is effected using the related ME block from 
  // the full grain solution, which is simply copied. 
  // Note upward transitions are determined as part of symmetrization.
  //
  template<class T> 
  void gWellProperties::constructReservoir(MesmerEnv& env, vector<double> &gEne, vector<double> &gDOS, TMatrix<T>* egme) {

    // Sum up the downward transition probabilities into the reservoir grain.
    T sumOfDeactivation(0.0), ptfReservoir(0.0);
    const T beta(env.beta) ;
    for (size_t j(0); j < m_ncolloptrsize; ++j) {
      if (j < m_numGroupedGrains){
        // Summing up the partition function of reservoir state.
        ptfReservoir += T(gDOS[j])*exp(-beta*T(gEne[j])) ;
      }
      else {
        T downwardSum(0.0);
        for (size_t i(0); i < m_numGroupedGrains; ++i) {
          downwardSum += (*egme)[i][j]; // Sum of the normalized downward prob.
        }
        T ptfj = T(gDOS[j])*exp(-beta*T(gEne[j])) ;
        sumOfDeactivation += downwardSum * ptfj;
        (*egme)[0][j - m_numGroupedGrains + 1] = downwardSum;
      }
    }
    sumOfDeactivation /= ptfReservoir; 
    (*egme)[0][0] = -sumOfDeactivation;

    // Shift active state block.
    for (size_t i(m_numGroupedGrains), ii(1); i < m_ncolloptrsize; ++i, ++ii) {
      for (size_t j(m_numGroupedGrains), jj(1); j < m_ncolloptrsize; ++j, ++jj) {
        (*egme)[ii][jj] = (*egme)[i][j];
      }
    }

  }

  //
  // Copy collision operator to diagonal block of system matrix.
  //
  template<class T> 
  void gWellProperties::copyCollisionOperator(qdMatrix *CollOptr, TMatrix<T> *egme, const size_t locate, const double RducdOmega) const
  {
    // Find size of system matrix.

    const size_t smsize = CollOptr->size();

    // Check there is enough space in system matrix.

    if (locate + m_ncolloptrsize > smsize)
      throw (std::runtime_error("Error in the size of the system matrix."));

    // Copy collision operator to the diagonal block indicated by "locate"
    // and multiply by the reduced collision frequencey.

    for (size_t i(0), ii(locate) ; i < m_ncolloptrsize; ++i, ++ii) {
      for (size_t j(0), jj(locate) ; j < m_ncolloptrsize; ++j, ++jj) {
        (*CollOptr)[ii][jj] = RducdOmega * (*egme)[i][j];
      }
    }
  }


}//namespace

#endif // GUARD_gWellProperties_h
