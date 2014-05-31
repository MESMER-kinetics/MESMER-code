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

    std::vector<double> m_eleExc  ;      // Electronic excitation (E.g. OH, NO, NS, otherwise no member).
    std::vector<double> m_VibFreq ;      // Values of vibrational frequencies.
    dMatrix            *m_Hessian ;      // Hessian matrix (If supplied, used to calculate vibrational frequncies).
    dMatrix            *m_Modes ;        // Vectors representing modes that are to be projected from Hessian.
	size_t              m_nModes ;       // Number of projected modes.
	std::string         m_HessianUnits ; // Hessian matrix units.
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
    void getEleExcitation(vector<double> &elecExci);

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
    
    //EnergyTransferModel* m_pEnergyTransferModel ; //loaded for each PT condition
    std::map<std::string, EnergyTransferModel*> m_EnergyTransferModels; //with different bath gases

    std::vector<double> m_grainDist ;          // Grain distribution (not normalized)
    dMatrix             *m_egme ;              // Matrix containing the energy grained collision operator.
    dMatrix             *m_egvec ;             // Eigenvectors used to diagonalize (P - I) matrix.
    std::vector<double>  m_egval;

    // Calculate collision frequency.
    double collisionFrequency(MesmerEnv env, Molecule *pBathGasMolecule) ;

    // Calculate collision operator.
    bool collisionOperator (MesmerEnv& env) ;

    // Calculate raw transition matrix.
    bool rawTransitionMatrix(MesmerEnv& env, vector<double> &gEne,  vector<double> &gDOS, dMatrix *egme) ;

    // Calculate collision operator with reservoir state.
    bool collisionOperatorWithReservoirState(MesmerEnv& env) ;

    double getBoltzmannWeightedEnergy(int numberOfGrains, const vector<double>& gEne, const vector<double>& gDos, double beta, double& totalDOS);

    // Diagonalize collision operator before adding reaction terms to get eigenvectors and eigenvalues.
    void diagonalizeCollisionOperator();

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

    void copyCollisionOperator(qdMatrix *CollOptr, const int locate, const double RducdOmega) const ;

    void copyCollisionOperatorEigenValues(qdMatrix *CollOptr, const int locate, const double RducdOmega) const ;

    void normalizedInitialDistribution(vector<double> &grainFrac) ;
    void normalizedGrnBoltzmannDistribution(vector<double> &grainFrac);
	void normalizedCellBoltzmannDistribution(vector<double> &grainFrac, const int totalCellNumber);

    // Accessors.

    double get_collisionFrequency() const {return m_collisionFrequency ; } ;

    void set_colloptrsize(int ncolloptrsize) ;
    size_t get_colloptrsize() const {return m_ncolloptrsize ; } ;

    size_t get_nbasis() const ;

    const int get_grnZPE();

    const dMatrix* getEigenVectors(){ return m_egvec; }

    const double getLowestBarrier() { return m_lowestBarrier;}
    void setLowestBarrier(double value){ m_lowestBarrier = value;}

    const int reservoirShift() {return m_numGroupedGrains == 0 ? 0 : m_numGroupedGrains - 1; }

  };


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
   void reducedMomentInertiaAngular(string bondID, double phase, vector<double>& angles, vector<double>& redInvMOI) ;

	// Calculate the internal rotation eigenvector. Based on the internal rotation 
	// mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010).
	// Typically this vector is used to project out an internal rotational mode
	// from a Hessian.
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

	// Export to xmol format.
	void exportToXYZ() const ;

  };


  //-------------------------------------------------------------------------------------------------
  // Other related functions
  //-------------------------------------------------------------------------------------------------

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS(gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, std::vector<double>& rctsCellDOS);

}//namespace

#endif // GUARD_MolecularComponents_h
