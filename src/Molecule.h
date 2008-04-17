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
#include "XMLPersist.h"
#include "MesmerMath.h"
#include "DensityOfStates.h"

namespace mesmer
{
  //**************************************************
  /// Basic molecule: has name and some collision parameters.
  /// Used for bath gases and unmodelled product molecules.
  class Molecule
  {
  private:

    const MesmerEnv&     m_Env;
    int            m_flag;              // count the errors during initialization
    PersistPtr     m_ppPersist;         // Conduit for I/O
    std::string    m_Name ;             // Molecule name.
    std::string    m_Description;       // Longer description for the structure
    double         m_Mass ;             // Mass.
    double         m_Sigma ;            // Lennard-Jones sigma.
    double         m_Epsilon ;          // Lennard-Jones epsilon.

    // CHECK FOR INPUTFILE PARAMETERS
    // for these check values, initially the values are given -1. If not provided by user the value will remain as -1.
    // If provided by the user the value will be increased to 0.
    // During calcualtion, every inputfile parameter if called by any function, the class will check if the parameter
    // is provided. If it is provided, the chk value will increase by 1. So if a chk value is 9, it is asked for 9 times.
    // However, if the user did not provide the value and the values is asked. The program will stop and report error.
    // Values with obvious default values are also accounted in this check but the program will not exit; only
    // warning message will be given.
    //================================================
    int m_Mass_chk;
    int m_Sigma_chk;
    int m_Epsilon_chk;
    //================================================

    //Molecule(const Molecule&) ;
    //Molecule& operator=(const Molecule&) ;

  public:
    Molecule(const MesmerEnv& Env) ;
    virtual ~Molecule();

    // Initialize Molecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    PersistPtr  getPersistentPointer();
    void setPersistentPointer(PersistPtr value);
    std::string getName() const;
    std::string getDescription() const;
    const MesmerEnv& getEnv() const;

    int    getFlag() ;
    double getMass() ;
    double getSigma() ;
    double getEpsilon() ;

    void   setName(string name) ;
    void   setFlag(bool value) ;
    void   setMass(double value);
    void   setSigma(double value);
    void   setEpsilon(double value);
  };

  //**************************************************
  class BathGasMolecule : public Molecule
  {
  public:
    BathGasMolecule(const MesmerEnv& Env) : Molecule(Env){}
    virtual ~BathGasMolecule(){} ;
    // Initialize BathGasMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);
  };

  //**************************************************
  ///Molecule with multiple energy states. Probably not instantiated itself
  class ModelledMolecule : public Molecule
  {
    //
    // Molecular properties.
    //
    double              m_RotCstA ;          // Moment of inertia A.
    double              m_RotCstB ;          // Moment of inertia B.
    double              m_RotCstC ;          // Moment of inertia C.
    double              m_Sym ;              // Rotational symmetry number.
    double              m_ZPE ;              // Zero Point Energy. (Kcal/Mol)
    double              m_scaleFactor ;      // scale factor for input real/imaginary vibrational frequencies
    int                 m_SpinMultiplicity ; // spin multiplicity
    int                 m_grnZpe ;           // Zero point energy expressed in grains.
    DensityOfStatesCalculator *m_pDensityOfStatesCalculator ;

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_RC_chk;
    int m_Sym_chk;
    int m_ZPE_chk;
    int m_scaleFactor_chk;
    int m_SpinMultiplicity_chk;
    int m_VibFreq_chk;
    int m_grnZpe_chk;
    //================================================

    std::vector<double> m_eleExc; //Electronic excitation(for OH, NO, NS otherwise no memeber).
    std::vector<double> m_VibFreq ;          // Values of vibrational frequencies.

  public:
    //
    // Cell and grain averages.  Note: raw array used for efficiency,
    // but need to test validity of this. SHR 5/Jan/03.
    //
    std::vector<double> m_cellEne ;          // Cell energy array.
    std::vector<double> m_cellDOS ;          // Cell density of states array.
    std::vector<double> m_grainEne ;         // Grain average energy array.
    std::vector<double> m_grainDOS ;         // Grain density of states array.
    
    //----------------
    ModelledMolecule(const MesmerEnv& Env);
    virtual ~ModelledMolecule();

    // Initialize ModelledMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    // Get the density of states.
    void getCellDensityOfStates(std::vector<double> &cellDOS) ;

    // Get cell energies.
    void getCellEnergies(std::vector<double> &CellEne) ;

    // Get grain density of states.
    void getGrainDensityOfStates(std::vector<double> &grainDOS) ;

    // Get grain energies.
    void getGrainEnergies(std::vector<double> &grainEne) ;

    // Get Grain Boltzmann distribution.
    void grnBoltzDist(std::vector<double> &grainBoltzDist) ;

    // Get Electronic excitations
    void getEleExcitation(vector<double> &elecExci);

    // Get Grain canonical partition function.
    double grnCanPrtnFn() ;

    // Calculate Density of states
    bool calcDensityOfStates();

    // Calculate classical energy
    double getClassicalEnergy();
    
    // Accessors.
    virtual double get_zpe();
    double get_scaleFactor();
    void set_zpe(double value);
    void set_scaleFactor(double value);
    double get_Sym(void);
    int test_rotConsts(void);
    int  get_rotConsts(std::vector<double> &mmtsInt);
    void set_grnZpe(int grnZpe); // with respect to the minimum of all wells, default zero.
    virtual const int get_grnZpe();
    virtual void get_VibFreq(std::vector<double>& vibFreq);

    virtual DensityOfStatesCalculator* get_DensityOfStatesCalculator();
    void set_DensityOfStatesCalculator(DensityOfStatesCalculator* value);

    virtual int getSpinMultiplicity();
    void   setSpinMultiplicity(int value);

    // Calculate the average grain energy and then number of states per grain.
    void calcGrainAverages() ;

    // Test the rovibrational density of states.
    virtual void testDensityOfStates() ;
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
    TransitionState(const MesmerEnv& Env);
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
    double              m_DeltaEdown ;         // <Delta E down> for the exponential down model.
    double              m_collisionFrequency ; // Current value of collision frequency.
    int                 m_ncolloptrsize ;      // Size of the collision operator matrix.

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_DeltaEdown_chk;
    //================================================

    dMatrix             *m_egme ;              // Matrix containing the energy grained collision operator.

    //-------------------------------
    // Calculate collision frequency.
    double collisionFrequency(double beta, const double conc, Molecule *pBathGasMolecule) ;
    // Calculate collision operator.
    bool   collisionOperator (double beta) ;

  public:
    CollidingMolecule(const MesmerEnv& Env);
    virtual ~CollidingMolecule();

    // Initialize CollidingMolecule.
    virtual bool InitializeMolecule(PersistPtr ppp);

    // Initialize the Collision Operator.
    bool initCollisionOperator(double beta, Molecule *pBathGasMolecule) ;

    // Calculate a reaction matrix element.
    double matrixElement(int eigveci, int eigvecj, std::vector<double> &k, int ndim) ;

    void copyCollisionOperator(dMatrix *CollOptr, const int size, const int locate, const double RducdOmega) const ;

    // Accessors.
    double get_collisionFrequency() const ;
    void set_colloptrsize(int ncolloptrsize) ;
    int  get_colloptrsize() const ;
    void   setDeltaEdown(double value);
    double getDeltaEdown();
  };

  // class representing pair of molecules participating one association reaction
  // this class should account for their density of states as a convolution
  class SuperMolecule : public ModelledMolecule
  {
  public:
    ModelledMolecule* m_mol1;
    ModelledMolecule* m_mol2;

    SuperMolecule(const MesmerEnv& Env);
    virtual ~SuperMolecule();

    // Initialize SuperMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    // Accessors.

    // set composing member of the SuperMolecule, also copy necessary properties
    int getSpinMultiplicity();
    void get_VibFreq(std::vector<double>& vibFreq);
    double get_zpe();
    DensityOfStatesCalculator* get_DensityOfStatesCalculator();
    void setMembers(ModelledMolecule* mol1p, ModelledMolecule* mol2p);

    ModelledMolecule* getMember1();
    ModelledMolecule* getMember2();

    //virtual void testDensityOfStates(const MesmerEnv &m_Env) ;

  };

}//namespace
#endif // GUARD_Molecule_h
