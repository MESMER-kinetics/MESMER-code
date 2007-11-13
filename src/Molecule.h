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
#include "MesmerEnv.h"
#include "MesmerMath.h"
#include "formatfloat.h"

namespace mesmer
{
  //**************************************************
  /// Basic molecule: has name and some collision parameters.
  /// Used for bath gases and unmodelled product molecules.
  class Molecule
  {
  public:
    Molecule() ;
    virtual ~Molecule(){} ;

    // Initialize Molecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    PersistPtr  getPersistentPointer()     { return m_ppPersist;};
    std::string getName() const            { return m_Name ; } ;

    int    getFlag()                       { return m_flag; } ;
    double getMass() const                 { return m_Mass ; } ;
    double getSigma() const                { return m_Sigma ; } ;
    double getEpsilon() const              { return m_Epsilon ; } ;

    void   setFlag(bool value)             { if (value) ++m_flag; } ;
    void   setMass(double value)           { m_Mass = value; } ;
    void   setSigma(double value)          { m_Sigma = value; } ;
    void   setEpsilon(double value)        { m_Epsilon = value; } ;

  private:
    int            m_flag;              // count the errors during initialization
    PersistPtr     m_ppPersist;         // Conduit for I/O
    std::string    m_Name ;             // Molecule name.
    double         m_Mass ;             // Mass.
    double         m_Sigma ;            // Lennard-Jones sigma.
    double         m_Epsilon ;          // Lennard-Jones epsilon.

    //Molecule(const Molecule&) ;
    //Molecule& operator=(const Molecule&) ;
  };

  //**************************************************
  class BathGasMolecule : public Molecule
  {
  public:
    // Initialize BathGasMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);
  };

  //**************************************************
  ///Molecule with multiple energy states. Probably not instantiated itself
  class ModelledMolecule : public Molecule
  {
  public:

    ModelledMolecule();
    virtual ~ModelledMolecule();

    // Initialize ModelledMolecule.
    virtual bool InitializeMolecule(PersistPtr pp);

    // Get the density of states.
    void cellDensityOfStates(std::vector<double> &cellDOS, const MesmerEnv &mEnv) ;

    // Get cell energies.
    void cellEnergies(std::vector<double> &CellEne, const MesmerEnv &mEnv) ;

    // Get grain density of states.
    void grnDensityOfStates(std::vector<double> &grainDOS, const MesmerEnv &mEnv) ;

    // Get grain energies.
    void grnEnergies(std::vector<double> &grainEne, const MesmerEnv &mEnv) ;

    // Get Grain Boltzmann distribution.
    void grnBoltzDist(std::vector<double> &grainBoltzDist, const MesmerEnv &mEnv) ;

    // Get Grain canonical partition function.
    double grnCanPrtnFn(const MesmerEnv &mEnv) ;

    // Accessors.
    double get_zpe() const { return m_ZPE ; }

    // Get moments of inertia
    int get_rotConsts(std::vector<double> &mmtsInt);

    double get_Sym(void){return m_Sym;}

    void get_VibFreq(std::vector<double>& vibFreq){
      vibFreq.clear();
      vibFreq = m_VibFreq;
    }

    int getSpinMultiplicity() const        { return m_SpinMultiplicity; }
    void   setSpinMultiplicity(int value)  { m_SpinMultiplicity = value; }

  protected:

    // Calculate the rovibrational density of states for 1 cm-1 cells.
    void calcDensityOfStates(const MesmerEnv &mEnv) ;

    // Calculate the average grain energy and then number of states per grain.
    void calcGrainAverages(const MesmerEnv &mEnv) ;

    // Test the rovibrational density of states.
    void testDensityOfStates(const MesmerEnv &mEnv) ;

  private:
    //
    // Molecular properties.
    //
    std::vector<double> m_VibFreq ;          // Values of vibrational frequencies.
    double              m_RotCstA ;          // Moment of inertia A.
    double              m_RotCstB ;          // Moment of inertia B.
    double              m_RotCstC ;          // Moment of inertia C.
    double              m_Sym ;              // Rotational symmetry number.
    double              m_ZPE ;              // Zero Point Energy.
    int                 m_SpinMultiplicity ; // spin multiplicity

  public:
    //
    // Cell and grain averages.  Note: raw array used for efficiency,
    // but need to test validity of this. SHR 5/Jan/03.
    //
    std::vector<double>  m_cellEne ;         // Cell mid-point energy array.
    std::vector<double>  m_cellDOS ;         // Cell density of states array.
    std::vector<double>  m_grainEne ;        // Grain average energy array.
    std::vector<double>  m_grainDOS ;        // Grain density of states array.

  } ;

  //**************************************************
  ///Transition states have no collisional properties
  class TransitionState : public ModelledMolecule
  {
  public:
    // Initialize TransitionState.
    TransitionState();

  };

  //**************************************************
  ///For most molecules
  class CollidingMolecule : public ModelledMolecule
  {
  public:
    CollidingMolecule();
    ~CollidingMolecule();

    // Initialize CollidingMolecule.
    virtual bool InitializeMolecule(PersistPtr ppp);

    // Initialize the Collision Operator.
    bool initCollisionOperator(double beta, Molecule *pBathGasMolecule, const MesmerEnv &mEnv) ;

    // Diagonalize the Collision Operator. See ReactionManager::diagCollisionOperator()
    //void diagCollisionOperator() ;

    // Calculate a reaction matrix element.
    double matrixElement(int eigveci, int eigvecj, std::vector<double> &k, int ndim) ;

    void copyCollisionOperator(dMatrix *CollOptr, const int size, const int locate, const double RducdOmega, const MesmerEnv &mEnv) const ;

    // Accessors.
    void set_colloptrsize(int ncolloptrsize) {m_ncolloptrsize = ncolloptrsize ;} ;
    int  get_colloptrsize() const {return m_ncolloptrsize ;} ;
    void set_grnZpe(int grnZpe) {m_grnZpe = grnZpe ;} ;
    int  get_grnZpe() const {return m_grnZpe ;} ;
    double get_collisionFrequency() const { return m_collisionFrequency ;} ;

  private:
    //
    // Collision operator properties.
    //
    dMatrix *m_egme ;        // Matrix containing the energy grained collision operator.
    double   m_DeltaEdown ;  // <Delta E down> for the exponential down model.

    // Calculate collision frequency.
    double collisionFrequency(double beta, const double conc, Molecule *pBathGasMolecule) ;
    // Calculate collision operator.
    bool   collisionOperator (double beta, const MesmerEnv &mEnv) ;

    int    m_grnZpe ;             // Zero point energy expressed in grains.
    int    m_ncolloptrsize ;      // Size of the collision operator matrix.
    double m_collisionFrequency ; // Current value of collision frequency.
  };

  // class representing pair of molecules participating one association reaction
  // this class should account for their density of states as a convolution
  class PairedMolecules
  {
  public:
    PairedMolecules();
    ~PairedMolecules();

  private:
    //pointers to two independent molecules participating one association reaction
    std::string    m_Name ;             // Molecule pair name.
    double         m_Mass ;             // Mass.

    ModelledMolecule* p_mol1;
    ModelledMolecule* p_mol2;

    double               m_ZPE ;             // Zero Point Energy.
    std::vector<double>  m_VibFreq ;         // Values of vibrational frequencies.
    std::vector<double>  m_cellEne ;         // Pointer to cell mid-point energy array.
    std::vector<double>  m_cellDOS ;         // Pointer to cell density of states array.
    std::vector<double>  m_grainEne ;        // Pointer to grain average energy array.
    std::vector<double>  m_grainDOS ;        // Pointer to grain density of states array.

  };

}//namespace
#endif // GUARD_Molecule_h
