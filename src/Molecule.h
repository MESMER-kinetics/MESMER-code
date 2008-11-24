#ifndef GUARD_Molecule_h
#define GUARD_Molecule_h

//-------------------------------------------------------------------------------------------
//
// Molecule.h
//
// Author: Struan Robertson
// Date:   5/Jan/2003
//
// This header file contains the declaration of the Molecule class.
//
//-------------------------------------------------------------------------------------------

#include <memory>
#include "MicroRate.h"
#include "DensityOfStates.h"
#include "Distribution.h"
#include "MesmerEnv.h"
#include "MesmerFlags.h"
#include "Pointer.h"

namespace mesmer
{

  //**************************************************
  /// Basic molecule: has name and some collision parameters.
  /// Used for bath gases and unmodelled product molecules.
  class Molecule
  {
  private:
    const MesmerEnv&     m_Env;
    MesmerFlags&         m_Flags;

    int            m_errorflag;         // count the errors during initialization
    PersistPtr     m_ppPersist;         // Conduit for I/O
    std::string    m_Name ;             // Molecule name.
    std::string    m_Description;       // Longer description for the structure

    //Molecule(const Molecule&) ;
    //Molecule& operator=(const Molecule&) ;

  public:
    Molecule(const MesmerEnv& Env, MesmerFlags& m_Flags) ;
    virtual ~Molecule();

    // Initialize Molecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    PersistPtr  getPersistentPointer();
    void setPersistentPointer(PersistPtr value);
    std::string getName() const;
    std::string getDescription() const;
    const MesmerEnv& getEnv() const;
    MesmerFlags& getFlags();
    int    getErrorFlag() ;
    void   setName(string name) ;
    void   setFlag(bool value) ;
  };

  //**************************************************
  class BathGasMolecule : public Molecule
  {
  public:
    BathGasMolecule(const MesmerEnv& Env, MesmerFlags& Flags);
    virtual ~BathGasMolecule();
    // Initialize BathGasMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    double getMass() ;
    void   setMass(double value);
    double getSigma() ;
    void   setSigma(double value);
    double getEpsilon() ;
    void   setEpsilon(double value);


  private:
    double         m_Mass ;             // Mass.
    double         m_Sigma ;            // Lennard-Jones sigma.
    double         m_Epsilon ;          // Lennard-Jones epsilon.

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_Mass_chk;
    int m_Sigma_chk;
    int m_Epsilon_chk;
    //================================================

  };

  //**************************************************
  ///Molecule with multiple energy states. Probably not instantiated itself
  class ModelledMolecule : public Molecule
  {

  private:
    //
    // Molecular properties.
    //
    double m_Mass ;             // Mass.
    double m_RotCstA ;          // Moment of inertia A.
    double m_RotCstB ;          // Moment of inertia B.
    double m_RotCstC ;          // Moment of inertia C.
    double m_Sym ;              // Rotational symmetry number.
    DPoint m_ZPE ;              // Zero Point Energy. (kJ/mol)
    double m_scaleFactor ;      // scale factor for input real/imaginary vibrational frequencies
    int    m_SpinMultiplicity ; // spin multiplicity

    double m_initPopulation ;   // initial population of the molecule.
    double m_eqFraction ;       // equilibrium fraction of the species

    DensityOfStatesCalculator *m_pDensityOfStatesCalculator ;

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_Mass_chk;
    int m_RC_chk;
    int m_Sym_chk;
    int m_ZPE_chk;
    int m_scaleFactor_chk;
    int m_SpinMultiplicity_chk;
    int m_VibFreq_chk;
    std::string m_EnergyConvention;

    //================================================

    std::vector<double> m_eleExc;   // Electronic excitation(for OH, NO, NS otherwise no member).
    std::vector<double> m_VibFreq ; // Values of vibrational frequencies.

  protected:

    //
    // Cell and grain averages.
    //
    std::vector<double> m_cellDOS ;   // Cell density of states array.
    std::vector<double> m_grainEne ;  // Grain average energy array.
    std::vector<double> m_grainDOS ;  // Grain density of states array.

  public:

    //----------------
    ModelledMolecule(const MesmerEnv& Env, MesmerFlags& Flags);
    virtual ~ModelledMolecule();

    // Initialize ModelledMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    // Get cell density of states.
    void getCellDensityOfStates(std::vector<double> &cellDOS, int startingCell = 0) ;

    // Set cell  density of states.
    void setCellDensityOfStates(std::vector<double> &cellDOS) { m_cellDOS = cellDOS ; } ;

    // Get grain density of states.
    void getGrainDensityOfStates(std::vector<double> &grainDOS, const int startGrnIdx = 0, const int ignoreCellNumber = 0) ;

    // Get grain energies.
    void getGrainEnergies(std::vector<double> &grainEne) ;

    // Get Electronic excitations
    void getEleExcitation(vector<double> &elecExci);

    // Get Grain canonical partition function.
    double rovibronicGrnCanPrtnFn() ;

    // Calculate Density of states
    bool calcDensityOfStates();

    // Calculate classical energy
    double getClassicalEnergy();

    // This function checks if any of the DPoint values is different then a DOS recalculation will take place
    bool needReCalculateDOS(void){
      if (!m_ZPE.isConstant()) return true;
      return false;
    }

    // This function explicitly tell all DPoint values in this ModelledMolecule that a DOS recalculation is completed.
    void recalculateDOScompleted(void){
      m_ZPE.updateValue();
    }

    // Accessors.
    double getMass() ;
    void   setMass(double value);
    virtual double get_zpe();
    double get_scaleFactor();
    void set_zpe(const double value){ m_ZPE = value; m_ZPE_chk = 0;};
    void set_zpe(const double valueL, const double valueU, const double stepsize){
      m_ZPE.set_range(valueL, valueU, stepsize);
      m_ZPE_chk = 0;
    }
    const std::string& getEnergyConvention()const { return m_EnergyConvention; }
    void set_scaleFactor(double value);
    double get_Sym(void);
    int test_rotConsts(void);
    int  get_rotConsts(std::vector<double> &mmtsInt);
    void set_grnZpe(double grnZpe); // with respect to the minimum of all wells, default zero.
    virtual void get_VibFreq(std::vector<double>& vibFreq);

    double getInitPopulation() const { return m_initPopulation;};
    void setInitPopulation(double value) { m_initPopulation = value;};
    double getEqFraction() const { return m_eqFraction;};
    void setEqFraction(double value){ m_eqFraction = value;};

    virtual DensityOfStatesCalculator* get_DensityOfStatesCalculator(){
      return m_pDensityOfStatesCalculator;
    }
    void set_DensityOfStatesCalculator(DensityOfStatesCalculator* value){
      m_pDensityOfStatesCalculator = value;
    }

    virtual int getSpinMultiplicity();
    void   setSpinMultiplicity(int value);

    // Test the rovibrational density of states.
    virtual void testDensityOfStates() ;

    void set_grainValues(double relativeZPE);
    int get_cellOffset(void) {
      double modulus = fmod(get_zpe() - getEnv().EMin, getEnv().GrainSize);
      if(modulus < 0.0)  // presently modulus is only less than 0 for the excess reactant in an association rxn
        modulus = 0.0;   // however, this problem should become obsolete once supermolecule DOS is calculated on the fly
      return int(modulus) ;
    } ;


  } ;

  //**************************************************
  ///Transition states have no collisional properties
  class TransitionState : public ModelledMolecule
  {
    double m_ImFreq;            // Imaginary frequency of this barrier (For tunneling in QM calculations)
    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_ImFreq_chk;
    //================================================

  public:
    TransitionState(const MesmerEnv& Env, MesmerFlags& Flags);
    virtual ~TransitionState();

    // Initialize TransitionState.
    virtual bool InitializeMolecule(PersistPtr pp);
    double get_ImFreq();

  };

  //**************************************************
  ///For most molecules
  class CollidingMolecule : public ModelledMolecule
  {

  private:

    //-------------------------------
    // Variables:
    // Collision operator properties.
    //
    double              m_Sigma ;            // Lennard-Jones sigma.
    double              m_Epsilon ;          // Lennard-Jones epsilon.

    double              m_DeltaEdownExponent;  // Exponent of <Delta E down> according to the relation
                                               // <delta_E_down>(T) = <delta_E_down>_ref * (T / m_DeltaEdownRefTemp)^n
                                               // where m_DeltaEdownExponent is the exponent n
                                               // By default, n = 0, which means delta_E_down does not depend on temperature.
    double              m_DeltaEdownRefTemp;   // reference temperature of <Delta E down>, default 298.
    DPoint              m_DeltaEdown ;         // <Delta E down> for the exponential down model.
    double              m_collisionFrequency ; // Current value of collision frequency.
    int                 m_ncolloptrsize ;      // Size of the collision operator matrix.
    DistributionCalculator* m_pDistributionCalculator;

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_Sigma_chk;
    int m_Epsilon_chk;
    int m_DeltaEdown_chk;
    //================================================

    double m_grainFracBeta;                    // beta used to calculate grain distribution fraction
    std::vector<qd_real> m_grainDist ;          // Grain distribution (not normalized)
    qdMatrix             *m_egme ;              // Matrix containing the energy grained collision operator.

    //-------------------------------
    // Calculate collision frequency.
    double collisionFrequency(double beta, const double conc, BathGasMolecule *pBathGasMolecule) ;
    // Calculate collision operator.
    bool   collisionOperator (double beta) ;

  public:
    CollidingMolecule(const MesmerEnv& Env, MesmerFlags& Flags);
    virtual ~CollidingMolecule();

    // Initialize CollidingMolecule.
    virtual bool InitializeMolecule(PersistPtr ppp);

    // Initialize the Collision Operator.
    bool initCollisionOperator(double beta, BathGasMolecule *pBathGasMolecule) ;

    // Normalize the Collision Operator.
    void normalizeCollisionOperator();

    // Calculate a reaction matrix element.
    double matrixElement(int eigveci, int eigvecj, std::vector<double> &k, int ndim) ;

    void copyCollisionOperator(qdMatrix *CollOptr, const int size, const int locate, const double RducdOmega) const ;

    // Get Grain Boltzmann distribution.
    void grainDistribution(vector<qd_real> &grainFrac, const int numberOfGrains);
    void normalizedInitialDistribution(vector<qd_real> &grainFrac, const int numberOfGrains) ;
    void normalizedGrnBoltzmannDistribution(vector<qd_real> &grainFrac, const int numberOfGrains);
    void normalizedCellBoltzmannDistribution(vector<qd_real> &cellFrac, const int startingCell = 0);

    // Accessors.
    double getSigma() ;
    void   setSigma(double value);
    double getEpsilon() ;
    void   setEpsilon(double value);
    double get_collisionFrequency() const ;
    void set_colloptrsize(int ncolloptrsize) ;
    int  get_colloptrsize() const ;
    void   setDeltaEdown(const double value){ m_DeltaEdown = value; m_DeltaEdown_chk = 0;};
    void   setDeltaEdown(const double valueL, const double valueU, const double stepsize){
      m_DeltaEdown.set_range(valueL, valueU, stepsize);
      m_DeltaEdown_chk = 0;
    };
    void   setDeltaEdownRefTemp(const double value){m_DeltaEdownRefTemp = value;};
    void   setDeltaEdownExponent(const double value){m_DeltaEdownExponent = value;};
    double getDeltaEdownRefTemp (){return m_DeltaEdownRefTemp; }
    double getDeltaEdownExponent(){return m_DeltaEdownExponent;}
    const int get_grnZPE();

    double getDeltaEdown();
    void set_DistributionCalculator(DistributionCalculator* value);
    DistributionCalculator* get_DistributionCalculator();
  };


}//namespace

// CHECK FOR INPUTFILE PARAMETERS
// for these check values, initially the values are given -1. If not provided by user the value will remain as -1.
// If provided by the user the value will be increased to 0.
// During calcualtion, every inputfile parameter if called by any function, the class will check if the parameter
// is provided. If it is provided, the chk value will increase by 1. So if a chk value is 9, it is asked for 9 times.
// However, if the user did not provide the value and the values is asked. The program will stop and report error.
// Values with obvious default values are also accounted in this check but the program will not exit; only
// warning message will be given.

#endif // GUARD_Molecule_h
