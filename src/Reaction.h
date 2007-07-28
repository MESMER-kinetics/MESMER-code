#ifndef GUARD_Reaction_h
#define GUARD_Reaction_h

//-------------------------------------------------------------------------------------------
//
// Reaction.h 
//
// Author: Struan Robertson 
// Date:   1/Feb/2003
//
// This header file contains the declaration of the Reaction class.
//
//-------------------------------------------------------------------------------------------

#include <vector>
#include "Molecule.h"
#include "MoleculeManager.h"
#include "Persistence.h"
#include "MicroRate.h"

namespace mesmer
{
  class MicroRateCalculator;

    class Reaction
    {
        //
        // Declare System as a friend class. Not sure that this is the best or 
        // most OOP way to go as it clearly defeats the point of encapsulation
        // but the System needs to know a lot about the reactions it is
        // combining. Review latter. SHR 2/Apr/2003.
        //

        friend class System ;

    public:

		// Type of reaction.
        typedef enum ReactionType{ ASSOCIATION, 
            DISSOCIATION,
            ISOMERIZATION,
            ERROR_REACTION } ;

		typedef std::map<CollidingMolecule *, int> isomerMap ;

        // Constructors.
        Reaction(){} ;

        Reaction(MoleculeManager *pMoleculeManager);

        // Destructor.
        ~Reaction() ;

        // Copy constructor.
        //   Reaction(const Reaction& reaction) ;

        // Assignment operator.
        //   Reaction& operator=(const Reaction& reaction) ;

        // Initialize reaction.
        bool Initialize(PersistPtr ppReac) ;

        std::string& getName() { return m_Name ; } ;

        // Modifier for reaction type.
        void put_Reactiontype(ReactionType reactiontype) ;

        // Accessor for reaction type.
        ReactionType get_Reactiontype() const { } ;

        // Add microcanonical terms to collision operator
        void AddMicroRates(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

        // Determine the equilibrium constant.
        void CalcEquilConst() { } ;

        // Access microcanoincal rate coeffcients. 
        void get_MicroRateCoeffs(std::vector<double> &kmc) ;

        double get_PreExp() const           { return m_PreExp; }
        double get_ActivationEnergy()const  { return m_ActEne; }

        // Reactant information:

        int get_NumberOfReactants() const { return m_Reactant2 ? 2 : 1 ; } ;
        int get_NumberOfProducts()  const { return m_Product2 ? 2 : 1 ; } ;
        void get_unimolecularspecies(std::vector<CollidingMolecule *> &unimolecularspecies) const ;
        TransitionState* get_TransitionState() const {return m_TransitionState;}
        PersistPtr get_PersistPtr(){return m_ppPersist;}

    private:

        // Read a molecule name from the XML file and look it up
        Molecule* GetMolRef(PersistPtr pp);

        // Grain average microcanonical rate coefficients.
        bool grnAvrgMicroRateCoeffs();

        // Wrapper function to calculate and grain average microcanoincal rate coeffcients. 
        bool calcGrnAvrgMicroRateCoeffs() ;

        // Add isomer reaction terms to collision matrix.
        void AddIsomerReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

        // Add (reversible) association reaction terms to collision matrix.
        void AddAssocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

        // Add dissociation reaction terms to collision matrix.
        void AddDissocReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) ;

        std::string        m_Name ;             // Reaction name.
        MoleculeManager   *m_pMoleculeManager ; // Pointer to molecule manager.

        //
        // Reaction composition.
        //
        CollidingMolecule *m_Reactant ;         // Reactant Molecule.
        Molecule          *m_Reactant2 ;        // Subsidiary reactant molecule.
        CollidingMolecule *m_Product ;          // Product Molecule.
        Molecule          *m_Product2 ;         // Subsidiary product molecule.
        TransitionState   *m_TransitionState;   // TransitionState
        ReactionType        m_reactiontype ;     // Type of reaction.
        //
        // Reaction Rate data.
        //
        double              m_kfwd ;             // Forward canonical (high pressure) rate coefficient.
        std::vector<double> m_kfmc ;             // Forward microcanonical rate coefficients.
        std::vector<double> m_kfgrn ;            // Grained averaged forward microcanonical rates.

        double              m_ActEne ;           // Activation Energy
        double              m_PreExp ;           // Preexponetial factor

        // I/O and control
        PersistPtr          m_ppPersist;         // Conduit for I/O

		//
		// Point to microcanoical rate coeff. calculator.
		//
		MicroRateCalculator *m_pMicroRateCalculator ;
    } ;

}//namespace
#endif // GUARD_Reaction_h
