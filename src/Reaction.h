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

#include "MoleculeManager.h"
#include "MicroRate.h"
#include "Tunneling.h"
#include "MesmerConfig.h"

namespace mesmer
{
    class Reaction
    {
    public:

		typedef std::map<CollidingMolecule*, int> isomerMap ;
        typedef std::map<SuperMolecule*    , int> sourceMap ;

		// Constructors.

        Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id);

        // Destructor.
        virtual ~Reaction();

        // Initialize reaction.
        virtual bool InitializeReaction(PersistPtr ppReac) ;

        // Add microcanonical terms to collision operator
        void AddMicroRates(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega);

        // Access microcanoincal rate coeffcients.
        void get_MicroRateCoeffs(std::vector<double> &kmc) ;

        const std::string& getName() const    { return m_Name ; } ;
        double getExcessReactantConc() const  { return m_ERConc ; } ;
        double get_PreExp() const             { return m_PreExp ; } ;
        double get_NInf()const                { return m_NInf   ; } ;
        double getHeatOfReaction() const      { return m_HeatOfReaction ; };
        const MesmerEnv& getEnv() const       { return m_Env; } ;
        void resetCalcFlag()                  { restartCalc = true; };

        // Calculate association reaction coefficients
        virtual void grainRateCoeffDetailedBalance(const int dir);

        // Get activiation energy
        double get_ActivationEnergy(void);

        TransitionState* get_TransitionState()const { return m_TransitionState ; } ;

        // Get unimolecualr species information:
        virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const = 0 ;

        // Product information:
        virtual SuperMolecule* get_bi_molecularspecies(void) const {return NULL ; } ;

        // Get the principal source reactant (i.e. reactant not in excess) if it exists.
        virtual ModelledMolecule *get_pseudoIsomer(void) const {return NULL ; } ;

		// Get the imaginary frequency of the transitions state.
		double get_TSImFreq(void) const {return m_TransitionState->get_ImFreq() ; } ;

		bool thereIsTunnelling (void) const {return (m_pTunnelingCalculator) ? true : false ; } ;

		void calculateTunnelingCoeffs(std::vector<double>& TunnelingProbability) {m_pTunnelingCalculator->calculateTunnelingCoeffs(this, TunnelingProbability); } ;

    protected:

        // Read a molecule name from the XML file and look it up
        Molecule* GetMolRef(PersistPtr pp);

        // Read parameters requires to determine reaction heats and rates.
        bool ReadRateCoeffParameters(PersistPtr ppReac) ;

        // I/O and control
        PersistPtr           m_ppPersist;            // Conduit for I/O
												  
        /*
        Each of the backward/forward microcanonical rate coefficients are based on 
		the well bottom of its own well. That is, each of these vectors contains rate 
		coefficients starting from index 0 (the bottom of that well) to the top 
		(MaximumGrain - ZpeOfTheWell). However, each of the vectors has the same number 
		of members with some members be zero (either not calculated or no reaction ocurred).
        */

		std::vector<double>  m_CellKfmc ;            // Forward  microcanonical rate coefficients.
        std::vector<double>  m_CellKbmc ;            // Backward microcanonical rate coefficients.
        std::vector<double>  m_GrainKfmc ;           // Grained averaged forward  microcanonical rates.
        std::vector<double>  m_GrainKbmc ;           // Grained averaged backward microcanonical rates.

		//
        // Reaction composition.
        //
        SuperMolecule       *m_srct ;                 // Reactant molecules as a super-reactant
        CollidingMolecule   *m_rct1 ;                 // Reactant Molecule.
        ModelledMolecule    *m_rct2 ;                 // Subsidiary reactant molecule.
        CollidingMolecule   *m_pdt1 ;                 // Product Molecule.
        ModelledMolecule    *m_pdt2 ;                 // Subsidiary product molecule.
        TransitionState     *m_TransitionState;       // TransitionState
		
        MoleculeManager     *m_pMoleculeManager ;     // Pointer to molecule manager.
        MicroRateCalculator *m_pMicroRateCalculator ; // Pointer to microcanoical rate coeff. calculator.
        TunnelingCalculator *m_pTunnelingCalculator ; // Pointer to Tunneling Calculator

	private:

        //   Reaction();

        // Copy constructor.
        //   Reaction(const Reaction& reaction) ;

        // Assignment operator.
        //   Reaction& operator=(const Reaction& reaction) ;

        // Grain averaging shorthand function for microcanoical rate coefficients
        void rateConstantGrnAvg(const int _MG, const int _gsz, const std::vector<double> &DOS, const std::vector<double> &CellRC, std::vector<double> &GrainRC);

        // Grain average microcanonical rate coefficients.
        bool grnAvrgMicroRateCoeffs();

        // Wrapper function to calculate and grain average microcanoincal rate coeffcients.
        bool calcGrnAvrgMicroRateCoeffs() ;

        // Add reaction terms to collision matrix.
        virtual void AddReactionTerms(dMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) = 0 ;

        // Calculate reaction equilibrium constant.
        virtual double calcEquilibriumConstant() ;

		const MesmerEnv& m_Env;
        std::string m_Name ;        // Reaction name.

        bool restartCalc;           // restart calculation on DOS
        double m_PreExp ;           // Preexponetial factor
        double m_NInf ;             // Modified Arrhenius parameter
        double m_ERConc ;           // Concentration of the excess reactant
        double m_kfwd ;             // Forward canonical (high pressure) rate coefficient.

        //
        // Reaction Rate data.
        //
        double m_HeatOfReaction ;   // The heat of reaction corrected for zero point energies.

    } ;


}//namespace
#endif // GUARD_Reaction_h
