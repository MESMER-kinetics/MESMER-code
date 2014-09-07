//-------------------------------------------------------------------------------------------
// MolecularComponents.h
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

#ifndef GUARD_MolecularComponents_h
#define GUARD_MolecularComponents_h

#include <memory>
#include "MicroRate.h"
#include "DensityOfStates.h"
#include "Distribution.h"
#include "MesmerEnv.h"
#include "MesmerFlags.h"
#include "Rdouble.h"
#include "EnergyTransferModel.h"
#include "vector3.h"
#include "dMatrix.h"
#include "Persistence.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{

  enum RotationalTop {
    LINEAR,
    NONLINEAR,
    SPHERICAL,
    OBLATE,
    PROLATE,
    ASYMMETRIC,
    UNDEFINED_TOP
  } ;

  // Forward class declarations.
  class Molecule;

  class MolecularComponent{
  public:
    Molecule* getHost() { return m_host; }
    const Molecule* getHost() const { return m_host; }

  protected:
    Molecule* m_host;
    MolecularComponent():m_host(NULL){}
  };

  class gBathProperties:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Bath gas related properties
    //-------------------------------------------------------------------------------------------------

  private:
    double         m_Sigma ;            // Lennard-Jones sigma.
    double         m_Epsilon ;          // Lennard-Jones epsilon.

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_Sigma_chk;
    int m_Epsilon_chk;
    //================================================

  public:

    //
    // Constructor, destructor and initialization
    //
    gBathProperties(Molecule* pMol);
    virtual ~gBathProperties();

    double getSigma() ;
    double getEpsilon() ;
    void   setSigma(double value);
    void   setEpsilon(double value);
  };

  class gDensityOfStates: public MolecularComponent
  {
    friend class Molecule ;
    //-------------------------------------------------------------------------------------------------
    // Cell density of states related properties
    //-------------------------------------------------------------------------------------------------

  private:
    std::vector<DensityOfStatesCalculator*> m_DOSCalculators;

    double m_RotCstA ;          // Moment of inertia A.
    double m_RotCstB ;          // Moment of inertia B.
    double m_RotCstC ;          // Moment of inertia C.
    double m_Sym ;              // Rotational symmetry number.

    Rdouble m_ZPE ;             // Zero Point Energy. (kJ/mol)

    double m_scaleFactor ;      // scale factor for input real/imaginary vibrational frequencies
    int    m_SpinMultiplicity ; // spin multiplicity

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_RC_chk;
    int m_Sym_chk;
    int m_ZPE_chk;
    int m_scaleFactor_chk;
    int m_SpinMultiplicity_chk;
    std::string m_EnergyConvention;
    //================================================

    std::vector<double> m_eleExc  ;      // Electronic excitation (E.g. OH, NO, NS).
    std::vector<double> m_VibFreq ;      // Values of vibrational frequencies.
    dMatrix            *m_Hessian ;      // Hessian matrix (If supplied, used to calculate vibrational frequncies).
    dMatrix            *m_Modes ;        // Vectors representing modes that are to be projected from Hessian.
	size_t              m_nModes ;       // Number of projected modes.
	std::string         m_HessianUnits ; // Hessian matrix units.
	size_t              m_MaximumCell ;  // Number of cells used in DOS calculations.

	//------------------------
    // Cell density of states
    //------------------------
    std::vector<double> m_cellDOS ;   // Cell density of states array.

    //------------------------
    // Grain density of states
    //------------------------
    std::vector<double> m_grainEne ;  // Grain average energy array.
    std::vector<double> m_grainDOS ;  // Grain density of states array.

    //
    // Constructor, destructor and initialization
    //
  public:
    gDensityOfStates(Molecule* pMol);
    ~gDensityOfStates();

    // Get the number of degrees of freedom for this species.
    unsigned int getNoOfDegOfFreeedom() ;

    // Get cell density of states. No recalculation if bcalc==false.
    bool getCellDensityOfStates(std::vector<double> &cellDOS, int startingCell = 0, bool bcalc=true) ;

    // Set cell  density of states.
    void setCellDensityOfStates(std::vector<double> &cellDOS) { m_cellDOS = cellDOS ; } ;

    // Get Electronic excitations
    void getEleExcitation(vector<double> &elecExci) { elecExci = m_eleExc ; } ;

    // Calculate Density of states
    bool calcDensityOfStates();

    // Calculate classical energy
    double getClassicalEnergy();

    // Accessors.
    double get_zpe();
    void set_zpe(const double value){ m_ZPE = value; m_ZPE_chk = 0;};
    void set_zpe(const double valueL, const double valueU, const double stepsize){
      m_ZPE.set_range(valueL, valueU, stepsize, "ZPE");
      m_ZPE_chk = 0;
    }

    std::string getEnergyConvention()const {
      return m_EnergyConvention.empty() ? "arbitary" : m_EnergyConvention;
    }

    double get_Sym(void);
    RotationalTop test_rotConsts(void);
    RotationalTop get_rotConsts(std::vector<double> &mmtsInt);
    void get_VibFreq(std::vector<double>& vibFreq);
    bool removeVibFreq(double freq); 

    int getSpinMultiplicity();

    int get_cellOffset(void);

    //----------------------------------
    // Grain density of states functions
    //----------------------------------

    // Get grain density of states.
    void getGrainDensityOfStates(std::vector<double> &grainDOS, const int startGrnIdx = 0, const int ignoreCellNumber = 0) ;

    // Get grain energies.
    void getGrainEnergies(std::vector<double> &grainEne) ;

    // Get Grain canonical partition function.
    double rovibronicGrnCanPrtnFn() ;

    // Calculate standard thermodynamic quantities as a function of temperature.
    bool thermodynamicsFunctions(double temp, double unitFctr, double& enthalpy, double& entropy, double& gibssFreeEnergy) ;

    bool RemoveDOSCalculator(const string& id);
    bool AddDOSCalculator(const string& id);
    DensityOfStatesCalculator* GetDOSCalculator(const string& id);

    // Get scale factor for vibrational frequencies
    double get_scaleFactor();

    // Methods for projecting out modes from the Hessian

    bool hasHessian() const { return m_Hessian ; } ;

    // This method is used to project a mode from the stored Hessian and
    // re-calculate the remaining frequencies.

    bool projectMode(std::vector<double> &mode) ;

  private:

    bool initialization() ;

    bool ReadDOSMethods();

    bool ReadZeroPointEnergy(PersistPtr &ppPropList) ;

    // This function checks if any of the DPoint values is different then a DOS recalculation will take place
    bool needReCalculateDOS(void){ return !m_ZPE.isUnchanged() ; }

    // This function explicitly tell all DPoint values in this Molecule that a DOS recalculation is completed.
    void recalculateDOScompleted(void){ m_ZPE.setUnchanged() ; }

    // Test the rovibrational density of states.
    void testDensityOfStates() ;

    // Calculate vibrational frequencies from molecular Hessian.
    bool FrqsFromHessian() ;

	// Helper function to create projector.
    void UpdateProjector(vector<double> &eigenvector) ;

    // Helper function to shift translation projection vector.
    void ShiftTransVector(vector<double> &eigenvector) ;

    // Function to calculate the rotational mode vectors.
    void RotationVector(vector<double> &aa, size_t loca, double sgna, vector<double> &bb, size_t locb, double sgnb, vector<double> &massWeights, vector<double> &mode) ;

    // Function to calculate the vibrational frequencies from a projected Hessian matrix.
    bool calculateFreqs(vector<double> &freqs, bool projectTransStateMode = false) ;

    // This method is used to orthogonalize a mode against existing projected modes
	// and then add it to the projected set.
    bool orthogonalizeMode(vector<double> &mode) ;

  };

  class gTransitionState:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Transition state related properties
    //-------------------------------------------------------------------------------------------------

  private:
    Rdouble m_ImFreq;            // Imaginary frequency of this barrier (For tunneling in QM calculations)
    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_ImFreq_chk;
    //================================================

  public:
    //
    // Constructor, destructor and initialization
    //
    gTransitionState(Molecule* pMol);
    virtual ~gTransitionState();

    double get_ImFreq();
    void set_imFreq(const double value){ m_ImFreq = value; m_ImFreq_chk = 0;}

  };

  class gPopulation:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Population and equilibrium fraction
    //-------------------------------------------------------------------------------------------------

  private:

    double m_initPopulation ;   // initial population of the molecule.
    double m_eqFraction ;       // equilibrium fraction of the species
    map <int,double> grainPopulations; // a map which holds any initial grain populations that have been specified

  public:

    //
    // Constructor, destructor and initialization
    //
    gPopulation(Molecule* pMol);

    double getInitPopulation() const { return m_initPopulation;};

    void setInitPopulation(double value) {
      if(grainPopulations.size()==0){
        m_initPopulation = value;
      } // only let the population be specified if there's no grain pop specified
      else{cerr << "initial grain population in this isomer has already been defined... ignoring population specifications " << endl;}
    };

    void getInitGrainPopulation(map<int,double>& inputMap) { 
      map<int,double>::iterator it;
      for(it=grainPopulations.begin(); it!=grainPopulations.end(); ++it){
        inputMap[it->first]=it->second;
      }
    };

    //  note: any given molecule should have EITHER a total population OR a grain population, BUT NOT BOTH
    void setInitGrainPopulation(int grain, double value) { 
      map<int,double>::iterator it;
      it = grainPopulations.find(grain);
      if(it==grainPopulations.end() && m_initPopulation==0){
        grainPopulations[grain] = value;
        m_initPopulation += value;
      }
      else if(it!=grainPopulations.end()){  // ignore redefinitions of grain populations
        cerr << "initial population of grain " << grain << " has been defined twice... ignoring redefinition " << endl;
      }
      else if(m_initPopulation!=0){ // only let the grain population be specified if there's no total population specified
        cerr << "initial population in this isomer has already been defined... ignoring grain population specifications " << endl;
      }
    };
    double getEqFraction() const { return m_eqFraction;};
    void setEqFraction(double value){ m_eqFraction = value;};

  };

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
	void normalizedCellBoltzmannDistribution(vector<double> &grainFrac, const int totalCellNumber);

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

  class gStructure:public MolecularComponent
  {
    //-------------------------------------------------------------------------------------------------
    // Chemical Structure related properties
    //-------------------------------------------------------------------------------------------------

  private:
    double m_MolecularWeight ;
    vector<double> m_PrincipalMI ; // amuAng2
	dMatrix *m_AxisAlignment ;

    struct atom
    {
      std::string id;
      std::string element;
      OpenBabel::vector3 coords;
      std::vector<std::string> connects; // Other atom ids.
    };
    std::map<std::string, atom> Atoms;

	std::map<std::string, std::pair<std::string, std::string> > Bonds;
    std::vector<std::string> m_atomicOrder ;
    bool m_HasCoords;

	bool m_verbose ; // Controls debug output.

	// Rotatable bonds
	vector<string> m_RotBondIDs ;

    enum AxisLabel {X = 0, Y = 1, Z = 2} ;

    // No default construction.
    gStructure();

    // Returns an ordered array of coordinates.
    void getAtomicCoords(vector<double> &coords, AxisLabel cartLabel) const ;

	// Method to shift coordinates to the centre of mass/principal axis frame. 
    bool AlignCoords() ;

    //Calculates moment of inertia of a set of atoms about an axis define by at1 and at2.
    double CalcMomentAboutAxis(std::vector<std::string> atomset, OpenBabel::vector3 at1, OpenBabel::vector3 at2);

    // Calculates internal rotation eigenvector about an axis define by at1 and at2.
    bool CalcInternalRotVec(std::vector<string> atomset, OpenBabel::vector3 at1, OpenBabel::vector3 at2, vector<double> &mode, bool ApplyMWeight) ;

    // Returns in atomset the IDs of all the atoms attached to atomID via bonds, but
    // does not include any atoms already in atomset or atoms beyond them.
    void GetAttachedAtoms(std::vector<std::string>& atomset, const std::string& atomID);

    // For a bond between atom at1 and atom at2 find all the atoms connected to at1
    // excluding those connect via at2.
    void findRotorConnectedAtoms(vector<string> &atomset, const string at1, const string at2) ;

    // Apply inertia weighting to the raw internal rotation velocity vector.
    void ApplyInertiaWeighting(vector<string> &atomset, vector<double> &velocity, double fctr) const ;

	// Calculates the GRIT for the current set of coordinates.
	double getGRIT(std::string bondID) ;

  public:

    gStructure(Molecule* pMol);

	~gStructure() { if(m_AxisAlignment) delete m_AxisAlignment ; } ;

	void set_Verbose(bool verbose) { m_verbose = verbose ; } ;

	//Returns true if atoms have coordinates
    bool ReadStructure();

    int NumAtoms() { return Atoms.size(); }

    bool IsAtom() {
      if(Atoms.empty())
        ReadStructure();
      return Atoms.size()==1;
    }

    double CalcMW();

    std::pair<std::string,std::string> GetAtomsOfBond(const std::string& bondID) {
      return Bonds[bondID];
    }

    OpenBabel::vector3 GetAtomCoords(const std::string atomID)
    {
      return Atoms[atomID].coords;
    }

    // Calculate the reduce moment of inertia about axis defined by specifed atoms.
    double reducedMomentInertia(pair<string,string>& bondats) ;

    // Calculate the angular dependent reduce moment of inertia about axis defined by specifed atoms.
   void reducedMomentInertiaAngular(string bondID, double phase, vector<double>& angles,
                                    vector<double>& redInvMOI, PersistPtr ppConfigData=NULL) ;

	// Calculate the internal rotation eigenvector. Based on the internal rotation 
	// mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010).
	// Typically this vector is used to project out an internal rotational mode
	// from a Hessian - in this situation mass weighting is required and this is
    // governed by the "ApplyWeight" flag. This vector is more recently been used 
    // in the calculation of internal rotation moments of inertia and in this 
    // case "ApplyWeight" should be false.
	void internalRotationVector(string bondID, vector<double>& mode, bool ApplyMWeight = true) ;

    // Read librarymols.xml to obtain the ab initio energy and enthalpy at zero K
    // for an isolated atom of each atom type in the molecule.
    // Return the sums of (E - Hf0) over the molecule in kJ/mol.
    // First parameter is true when atom-based thermochemistry is used, see DOI: 10.1002/chem.200903252
    // If useHf298 is true, calculates sum over (E - Hf298 + H0-H298)
    double CalcSumEMinusHf0(bool UsingAtomBasedThermo, bool useHf298);

    //Calculate moment of inertia matrix
    vector<double> CalcRotConsts(); 

    double getMass() const { return m_MolecularWeight;};

    void setMass(double value) { m_MolecularWeight = value;};

    int getAtomicOrder(std::string AtomID) const { 
      size_t i(0) ; 
      for  (; AtomID != m_atomicOrder[i] && i < m_atomicOrder.size() ; i++ ) ;		
      return (i < m_atomicOrder.size()) ? int(i) : -1 ;
    } ;

	// Returns an ordered array of masses.
    void getAtomicMasses(vector<double> &AtomicMasses) const ;

    // Returns an ordered array of X coordinates.
    void getXCoords(vector<double> &coords) const ;

    // Returns an ordered array of Y coordinates.
    void getYCoords(vector<double> &coords) const ;

    // Returns an ordered array of Z coordinates.
    void getZCoords(vector<double> &coords) const ;

	// Returns the alignment matix if it exists
	void getAlignmentMatrix(dMatrix &rAlignmentMatrix) const { 
	  if (m_AxisAlignment) 
		rAlignmentMatrix = *m_AxisAlignment ;
	} ;

	// Add rotatable bond ID (needed to calculate GRIT).
	void addRotBondID(std::string id) { m_RotBondIDs.push_back(id) ; } ;

    // Export to xmol and CML format.
    void exportToXYZ(const char* txt=NULL, bool last=false, PersistPtr ppConfigData=NULL) ;
    void exportToCML(const char* txt=NULL, bool last=false, PersistPtr ppConfigData=NULL) ;
  };


  //-------------------------------------------------------------------------------------------------
  // Other related functions
  //-------------------------------------------------------------------------------------------------

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS(gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, std::vector<double>& rctsCellDOS);

}//namespace

#endif // GUARD_MolecularComponents_h
