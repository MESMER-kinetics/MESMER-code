#ifndef GUARD_gDensityOfStates_h
#define GUARD_gDensityOfStates_h


#include "MolecularComponents.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{

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

    Rdouble m_ZPE ;             // Zero Point Energy. //wavenumbers
    double m_scaleFactor ;      // scale factor for input real/imaginary vibrational frequencies
    int    m_SpinMultiplicity ; // spin multiplicity

    //================================================
    // CHECK FOR INPUTFILE PARAMETERS
    int m_RC_chk;
    int m_Sym_chk;
    int m_ZPE_chk;
    int m_scaleFactor_chk;
    int m_SpinMultiplicity_chk;
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
    double get_zpe(); //cm-1
    void set_zpe(const double value){ m_ZPE = value; m_ZPE_chk = 0;};
    void set_zpe(const double valueL, const double valueU, const double stepsize){
      m_ZPE.set_range(valueL, valueU, stepsize, "ZPE");
      m_ZPE_chk = 0;
    }

    double get_Hf298Thermo(); //kJ/mol

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

    // This method tests if a rotor is heavy. It is a helper method
	// used to assess if a QM method will be expensive for calculating
	// the energy levels of an asymmetic top.

	bool IsHeavyTop(size_t n) ;

  private:

    bool initialization() ;

    bool ReadDOSMethods();

    bool ReadZeroPointEnergy(PersistPtr ppPropList) ;

    double ConvertEnergyConvention(
      const std::string& fromConvention, const std::string& toConvention, double fromValue);

    //Read <me:ZPE>, <me:Hf0>, <me:Hf298> or <me:HfAT0> from XML as specified by elName.
    //Returns true if found. Converts energy convention if the XML and nativeConvention are different .
    bool ReadEnergy(std::string elName, std::string nativeConvention);

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

  }//namespace

#endif //GUARD_gDensityOfStates_h

