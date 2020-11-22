#ifndef GUARD_gStructure_h
#define GUARD_gStructure_h

#include "MolecularComponents.h"

namespace mesmer
{
  class gStructure :public MolecularComponent
  {
    //-------------------------------------------------------------------------------------------------
    // Chemical Structure related properties
    //-------------------------------------------------------------------------------------------------

  private:
    double m_MolecularWeight;
    vector<double> m_PrincipalMI; // amuAng2
    dMatrix *m_AxisAlignment;

    struct atom
    {
      std::string id;
      std::string element;
      OpenBabel::vector3 coords;
      std::vector<std::string> connects; // Other atom ids.
    };
    std::map<std::string, atom> Atoms;

    std::map<std::string, std::pair<std::string, std::string> > Bonds;
    std::vector<std::string> m_atomicOrder;
    bool m_HasCoords;

    bool m_verbose; // Controls debug output.

    // Rotatable bonds
    vector<string> m_RotBondIDs;

    enum AxisLabel { X = 0, Y = 1, Z = 2 };

    // No default construction.
    gStructure();

    // Returns an ordered array of coordinates.
    void getAtomicCoords(vector<double> &coords, AxisLabel cartLabel) const;

    // Method to shift coordinates to the centre of mass/principal axis frame. 
    bool AlignCoords();

    //Calculates moment of inertia of a set of atoms about an axis define by at1 and at2.
    double CalcMomentAboutAxis(std::vector<std::string> atomset, OpenBabel::vector3 at1, OpenBabel::vector3 at2);

    // Calculates internal rotation eigenvector about an axis define by at1 and at2.
    bool CalcInternalRotVec(std::string bondID, std::vector<string> atomset, OpenBabel::vector3 at1, OpenBabel::vector3 at2, vector<double> &mode, bool ApplyMWeight);

    // Calculates bend eigenvector about an axis define by at1 and at2.
    double CalcBendVec(std::vector<string> atomset, OpenBabel::vector3 at1, OpenBabel::vector3 at2, vector<double>& mode);

    // Calculates stretch eigenvector along an axis define by at1 and at2.
    double CalcStretchVec(vector<string> atomset, OpenBabel::vector3 axis, vector<double>& mode);

    // Returns in atomset the IDs of all the atoms attached to atomID via bonds, but
    // does not include any atoms already in atomset or atoms beyond them.
    void GetAttachedAtoms(std::vector<std::string>& atomset, const std::string& atomID);

    // For a bond between atom at1 and atom at2 find all the atoms connected to at1
    // excluding those connect via at2.
    void findRotorConnectedAtoms(vector<string> &atomset, const string at1, const string at2);

    // Apply inertia weighting to the raw internal rotation velocity vector.
    void ApplyInertiaWeighting(vector<string> &atomset, vector<double> &velocity, double fctr) const;

		// Rotate one part of the molecule relative to the other by a given angle.
		void applyTorsionAngle(string bondID, double angle);

    // Calculates the GRIT for the current set of coordinates.
		void getGRIT(dMatrix& GRIT);

		// Calculates the effective reciprocal moment of inertia around
		// a given bond for a given configuration.
		double getReciprocalMOI(std::string bondID);
	
    // Helper functions for reduced mass methods:

    // Calculate vectors relative to the center of motion.
    OpenBabel::vector3 setVectorDirection(std::pair<std::string, std::string>& bondats, std::string& centre);

    // Orthogonalize a mode against the translationl and rotational modes.
    bool OrthogonalizeMode(vector<double>& mode);

    // Applies Gram-Schmit projection to remove component of second vector from the first.
    bool ProjectOutMode(vector<double>& mode, vector<double>& tmode);

    bool sequenceSearch(const string AtomID, vector<string> &stk, size_t idx);

  public:

    gStructure(Molecule* pMol);

    ~gStructure() { if (m_AxisAlignment) delete m_AxisAlignment; };

    void set_Verbose(bool verbose) { m_verbose = verbose; };

    //Returns true if atoms have coordinates
    bool ReadStructure();

    int NumAtoms() { return Atoms.size(); }

    bool IsAtom() {
      if (Atoms.empty())
        ReadStructure();
      return Atoms.size() == 1;
    }

    double CalcMW();

    std::map<std::string, int> GetElementalComposition() const;

    std::pair<std::string, std::string> GetAtomsOfBond(const std::string& bondID) {
      return Bonds[bondID];
    }

    OpenBabel::vector3 GetAtomCoords(const std::string atomID)
    {
      return Atoms[atomID].coords;
    }

    // Calculate the reduce moment of inertia about axis defined by specifed atoms.
    double reducedMomentInertia(pair<string, string>& bondats);

    // Calculate the angular dependent reduce moment of inertia about axis defined by specifed atoms.
    void reducedMomentInertiaAngular(string bondID, double phase, vector<double>& angles,
      vector<double>& redInvMOI, PersistPtr ppConfigData = NULL);

    // Calculate the internal rotation eigenvector. Based on the internal rotation 
    // mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010).
    // Typically this vector is used to project out an internal rotational mode
    // from a Hessian - in this situation mass weighting is required and this is
    // governed by the "ApplyWeight" flag. This vector is more recently been used 
    // in the calculation of internal rotation moments of inertia and in this 
    // case "ApplyWeight" should be false.
    void internalRotationVector(string bondID, vector<double>& mode, bool ApplyMWeight = true);

		// Calculate the square root of the determinant of the Generalized Rotation Inertia
		// Tensor. Used in the calculation of coupled rotor density of states and partition
		// functions.
		double getSqrtGRITDeterminant(vector<double>& angles);

    // Read librarymols.xml to obtain:
    //  the ab initio energy of the molecule, 
    //  the sums over the constituent atoms of:
    //  their enthalpy of formation at 0K and at 298K, 
    //  their enthalpy difference between 0K and 298K and
    //  their integrated Heat capacity between 0K and 298K.
    // All in kJ/mol. Return false on error.
    bool GetConstituentAtomThermo
    (double& ZPE, double& Hf0, double& Hf298, double& dH298, double& dStdH298);

    //Calculate moment of inertia matrix
    vector<double> CalcRotConsts();

    double getMass() const { return m_MolecularWeight; };

    void setMass(double value) { m_MolecularWeight = value; };

    int getAtomicOrder(std::string AtomID) const {
      size_t i(0);
      for (; i < m_atomicOrder.size() && AtomID != m_atomicOrder[i] ; i++);
      return (i < m_atomicOrder.size()) ? int(i) : -1;
    };

    // Returns an ordered array of masses.
    void getAtomicMasses(vector<double> &AtomicMasses) const;

    // Returns an ordered array of X coordinates.
    void getXCoords(vector<double> &coords) const;

    // Returns an ordered array of Y coordinates.
    void getYCoords(vector<double> &coords) const;

    // Returns an ordered array of Z coordinates.
    void getZCoords(vector<double> &coords) const;

    // Returns the alignment matix if it exists
    void getAlignmentMatrix(dMatrix &rAlignmentMatrix) const {
      if (m_AxisAlignment)
        rAlignmentMatrix = *m_AxisAlignment;
    };

    // Add rotatable bond ID (needed to calculate GRIT).
    void addRotBondID(std::string id) { m_RotBondIDs.push_back(id); };

		// Set rotatable bonds.
		void setRotBondID(std::vector<std::string>& BondIDs) { m_RotBondIDs = BondIDs ; };

		// Get rotatable bonds.
		void getRotBondID(std::vector<std::string>& BondIDs) { BondIDs = m_RotBondIDs ; };

    // Export to xmol and CML format.
    void exportToXYZ(const char* txt = NULL, bool last = false, PersistPtr ppConfigData = NULL);
    void exportToCML(const char* txt = NULL, bool last = false, PersistPtr ppConfigData = NULL);

    // Reduced mass methods:
    double bondStretchReducedMass(vector<string>& bondIDs, vector<double> &mode);

    double angleBendReducedMass(vector<string>& bondIDs, vector<double>& mode);

    double inversionReducedMass(vector<string>& bondIDs, vector<double>& mode);

    bool findFunctionalForm(string functionalForm) ;

  };

} //namespace
#endif //GUARD_gStructure_h
