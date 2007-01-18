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

#include <vector>
#include <memory>
#include "dMatrix.h"
#include "Persistence.h"

namespace mesmer
{
    static const int MAXCELL = 50000 ;
    static const int MAXGRN  = 500 ;

    //**************************************************
    /// Basic molecule: has name and some collision parameters.
    /// Used for bath gases and unmodelled product molecules.
    class Molecule : public IPersistObject
    {
    public:
        Molecule() ;
        virtual ~Molecule(){} ;

        // Initialize Molecule.
        virtual bool ReadFromXML(TiXmlElement* pnMol) ;

        // Get Molecule Name.
        std::string getName() { return m_Name ; } ;
        double getMass(){ return m_Mass ; } ;
        double getSigma(){ return m_Sigma ; } ;
        double getEpsilon(){ return m_Epsilon ; } ;

    protected:
        std::string   m_Name ;    // Molecule name.
        TiXmlElement* m_pXmlEl ;  // Its XML node address
        double        m_Mass ;    // Mass.
        double        m_Sigma ;   // Lennard-Jones sigma.
        double        m_Epsilon ; // Lennard-Jones epsilon.

    private:
        //   Molecule(const Molecule&) ;
        //   Molecule& operator=(const Molecule&) ;

    };

    //**************************************************
    class BathGasMolecule : public Molecule
    {
    public:
        // Initialize Molecule.
        virtual bool ReadFromXML(TiXmlElement* pnMol) ;
    };

    //**************************************************
    ///Molecule with multiple energy states. Probably not instantiated itself
    class ModelledMolecule : public Molecule
    {
    public:
        ModelledMolecule();
        virtual ~ModelledMolecule();
        // Initialize Molecule.
        virtual bool ReadFromXML(TiXmlElement* pnMol) ;

        // Get the density of states.
        void densityOfStates(double *) ;

        // Get cell energies.
        void cellEnergies(double *) ;

        // Accessors.
        double get_zpe() const { return m_ZPE ; } ;

    protected:

        // Calculate the rovibrational density of states for 1 cm-1 cells.
        void calcDensityOfStates() ;

        // Calculate the average grain energy and then number of states per grain.
        void calcGrainAverages() ;

        // Test the rovibrational density of states.
        void testDensityOfStates() ; 

        //
        // Memory management.
        //
        std::allocator<double>        m_alloc ;

    protected:
        //
        // Molecular properties.
        //
        std::vector<double> m_VibFreq ;     // Values of vibrational frequencies.
        double              m_MmtIntA ;     // Moment of inertia A.
        double              m_MmtIntB ;     // Moment of inertia B.
        double              m_MmtIntC ;     // Moment of inertia C.
        double              m_Sym ;         // Rotational symmetry number.
        double              m_ZPE ;         // Zero Point Energy.

        //
        // Cell and grain averages.  Note: raw array used for efficiency, 
        // but need to test validity of this. SHR 5/Jan/03.
        //
        int                  m_MaxCell ;     // Maximum number of cells to use.
        int                  m_MaxGrn  ;     // Maximum number of grains.
        int                  m_GrnSz ;       // Grain size in cm-1.
        double              *m_ecll ;        // Pointer to cell mid-point energy array.
        double              *m_cdos ;        // Pointer to cell density of states array.
        std::vector<double>  m_egrn ;        // Pointer to grain average energy array.
        std::vector<double>  m_gdos ;        // Pointer to grain density of states array.

    } ;

    //**************************************************
    ///Transition states have no collisional properties
    class TransitionState : public ModelledMolecule
    {
    public:
        TransitionState();

    };
    //**************************************************
    ///For most molecules
    class CollidingMolecule : public ModelledMolecule
    {
    public:
        CollidingMolecule();
        ~CollidingMolecule();

        // Initialize Molecule.
        virtual bool ReadFromXML(TiXmlElement* pnMol) ;

        // Initialize the Collision Operator.
        void initCollisionOperator(double beta) ;

        // Diagonalize the Collision Operator.
        void diagCollisionOperator() ;

        // Calculate collision frequency.
        double collisionFrequency(double temp, double conc, double bthMass, 
            double bthSigma, double bthEpsilon) ;

        // Calculate a reaction matrix element.
        double matrixElement(int eigveci, int eigvecj, std::vector<double> &k, int ndim) ;

        // Accessors.
        void copyCollisionOperator(dMatrix *CollOptr, const int size, const int locate) const ;

    protected:
        //
        // Collision operator properties.
        // 
        dMatrix *m_egme ;        // Matrix containing the energy grained collision operator.
        double   m_DeltaEdown ;  // <Delta E down> for the exponential down model.
    };


}//namespace
#endif // GUARD_Molecule_h
