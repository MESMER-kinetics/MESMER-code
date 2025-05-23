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
#include "dMatrix.h"

namespace mesmer
{

  enum ReactionType {
    ISOMERIZATION,
    ASSOCIATION,
    DISSOCIATION,
    IRREVERSIBLE_ISOMERIZATION,
    IRREVERSIBLE_EXCHANGE,
    BIMOLECULAR_SINK,
    BIMOLECULAR_EXCHANGE,
    PSEUDOISOMERIZATION,
    SECONDORDERASSOCIATION,
    DIFFUSION,
    UNDEFINED_REACTION
  };

  class MicroRateCalculator;

  class TunnelingCalculator;

  class Reaction
  {
  public:
    //Orders Molecule pointers using the Molecule name.
    //Using the pointer itself was seen to give unpredictable ordering.
    //See Scott Meyers "Effective STL", Item 20
    struct MoleculePtrLess : public std::function<bool(const Molecule*, const Molecule*)>
    {
      bool operator()(const Molecule* mol1, const Molecule* mol2)const
      {
        return mol1->getName() < mol2->getName();
      }
    };
    struct ReactionPtrLess : public std::function<bool(const Reaction*, const Reaction*)>
    {
      bool operator()(const Reaction* r1, const Reaction* r2)const
      {
        return r1->getName() < r2->getName();
      }
    };

    typedef std::map<Molecule*, int, MoleculePtrLess> molMapType;

    // Constructors.

    Reaction(MoleculeManager* pMoleculeManager, const MesmerEnv& Env, MesmerFlags& Flags, const char* id);

    // Destructor.
    virtual ~Reaction();
    virtual void Finish() {} // Mostly does nothing except AssociationReaction restores m_ZPEs of reactants

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) = 0;
    PersistPtr get_PersistentPointer()const { return m_ppPersist; }

    const std::string& getName() const { return m_Name; }

    double getHeatOfReaction() const {
      const double pdtZPE = get_relative_pdtZPE();
      const double rctZPE = get_relative_rctZPE();
      return pdtZPE - rctZPE;
    };
    int getHeatOfReactionInt() const { return int(getHeatOfReaction()); }
    const MesmerEnv& getEnv() const { return m_Env; };
    MesmerFlags& getFlags() { return m_Flags; };
    void resetCalcFlag() { m_reCalcMicroRateCoeffs = true; };

    // return reactant and product zero-point energy
    virtual double get_relative_rctZPE(void) const = 0;
    virtual double get_relative_pdtZPE(void) const = 0;
    virtual double get_relative_TSZPE(void) const = 0;

    // Get threshold energy
    virtual double get_ThresholdEnergy(void);
    /* This function should be considered as a function to get Einf.
    In ILT, not the theoretical threshold energy but the experimental Einf is used.
    This function returns user defined m_EInf, otherwise zero.
    ILT can be used in all reaction types if necessary. */

    // get products and reactants
    virtual size_t get_products(std::vector<Molecule*>& product) const = 0;
    virtual size_t get_reactants(std::vector<Molecule*>& reactants) const = 0;

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual Molecule* get_reactant(void) const = 0;

    Molecule* get_TransitionState() const { return m_TransitionState; };

    // Get unimolecular species information:
    virtual void get_unimolecularspecies(std::vector<Molecule*>& unimolecularspecies) const = 0;

    enum reactionType { all, rev, reactantsOnly, productsOnly };
    std::string getReactionString(reactionType = all);

    // Get the imaginary frequency of the transitions state.
    double get_TSImFreq(void) const;

    bool thereIsTunnelling(void) const { return (m_pTunnelingCalculator) ? true : false; };

    // Get tunnelling probabilities if they are defined.
    void calculateCellTunnelingCoeffs(std::vector<double>& TunnelingProbability);

    // calculate flux in grains
    void fluxCellToGrain();

    // returns the flux in cells for foreign modifications
    std::vector<double>& get_CellFlux(void) { return m_CellFlux; };

    // returns the forward grain microcanoincal rate coefficients for foreign modifications
    const std::vector<double>& get_GrainKfmc(void) { return m_GrainKfmc; };

    // returns the forward grain microcanoincal rate coefficients for foreign modifications
    const std::vector<double>& get_MtxGrnKf(void) { return m_MtxGrnKf; };

    // get canonical pseudo first order irreversible loss rate coefficient
    virtual double GetCanonicalIrreversibleLossRate(void) { return 0.0; };

    // set the bottom energy of m_CellFlux
    void setCellFluxBottom(const double energyValue);

    // return the grain idx in flux where the forward & reverse kofEs begin, respectively
    void calcFluxFirstNonZeroIdx(void);

    // get the grain in flux vector which corresponds to the threshold energy
    // normally this is the first grain, except for cases where the threshold energy is negative
    int get_fluxFirstNonZeroIdx(void) { return int(m_GrnFluxFirstNonZeroIdx); };

    // set & get flux Start Idx for calculating k(e)s from flux
    void set_EffGrnFwdThreshold(int idx) { m_EffGrainedFwdThreshold = idx; };
    int get_EffGrnFwdThreshold(void) { return int(m_EffGrainedFwdThreshold); };

    // set & get the forward threshold energy for calculating backward k(e)s from flux
    void set_EffGrnRvsThreshold(int idx) { m_EffGrainedRvsThreshold = idx; };
    int get_EffGrnRvsThreshold(void) { return int(m_EffGrainedRvsThreshold); };

    // get the backward threshold energy for calculating backward k(e)s from flux
    int get_fluxGrnZPE(void) { return int(m_FluxGrainZPE); };
    int get_fluxZPE(void) { return int(m_FluxCellZPE); };

    // calculate the effective threshold energy for utilizing in k(E) calculations, necessary for cases
    // with a negative threshold energy
    virtual void calcEffGrnThresholds(void) = 0;

    // Get the bottom cell offset of m_CellFlux.
    const size_t getFluxCellOffset(void) { return m_FluxCellOffset; };

    // Wrapper function to calculate and grain average microcanoincal rate coeffcients.
    virtual bool calcGrnAvrgMicroRateCoeffs();

    // Calculate high pressure rate coefficients.
    virtual bool HighPresRateCoeffTest(PersistPtr ppbase);

    static void setTestInterval(double TMin, double TMax, double dTemp) {
      m_dTemp = dTemp;
      m_TMin = TMin;
      m_TMax = TMax;
    };

    static void getTestInterval(double& TMin, double& TMax, double& dTemp) {
      dTemp = m_dTemp;
      TMin = m_TMin;
      TMax = m_TMax;
    };

    // Add reaction terms to the reaction matrix.
    virtual void AddReactionTerms(qdMatrix* CollOptr, molMapType& isomermap, const double rMeanOmega) = 0;

    // Add contracted basis set reaction terms to the reaction matrix.
    virtual void AddContractedBasisReactionTerms(qdMatrix* CollOptr, molMapType& isomermap) = 0;

    // Is reaction equilibrating and therefore contributes
    // to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double& Keq, Molecule** rct, Molecule** pdt) { return false; };

    // returns the reaction type
    virtual ReactionType getReactionType() = 0;

    // Calculate high pressure rate coefficients at current T.
    virtual void HighPresRateCoeffs(vector<double>* pCoeffs) = 0;

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() = 0;

    // For reactions involving a source update pseudoisomer map.
    virtual void updateSourceMap(molMapType& sourcemap) {
      /* For reactions without source terms this is a NULL operation. */
    };

    // Get the concentration of the excess reactant. 
    double get_concExcessReactant() const { return (m_bERConcPercent) ? m_ERConc*m_Env.conc : m_ERConc; }
    void   set_concExcessReactant(double conc, bool bPercentConc = false) {
      m_ERConc = conc; m_bERConcPercent = bPercentConc;
    }
    Molecule* getExcessReactant() { return m_ExcessReactant; }

    // Test if excess concentration is greater than bath gas concentration.
    // This method will be over ridden only by association reactions and
    // derivatives thereof.
    virtual bool IsExcessGreaterThanBathConc(double bathConc) const { return false; }

    // The following method takes an effective unimolecular rate coefficient 
    // and if required, normalizes it by, a concentration and/or any other 
    // factors in order to obtain a second order rate coefficient. For a
    // unimolecular reaction this function does nothing. This function is 
    // used in the calculation of Chi^2 values and analytical representation.
    virtual void normalizeRateCoefficient(double& rateCoefficient, std::string ref = "") const = 0;

    void setUsesProductProperties(bool b = true);
    bool UsesProductProperties() const { return m_UsesProductProperties; }

  protected:

    // Read a molecule name from the XML file and look it up
    // The defaultType is used if there is no me:type or role attribute
    Molecule* GetMolRef(PersistPtr pp, const char* defaultType = NULL);

    // Read parameters requires to determine reaction heats and rates.
    bool ReadRateCoeffParameters(PersistPtr ppReac);

    // Read excess reactant concentration
    bool ReadExcessReactantConcentration(PersistPtr ppReac);

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs() = 0;

    // I/O and control
    PersistPtr           m_ppPersist;            // Conduit for I/O

    //
    // Reaction composition.
    //

    Molecule* m_TransitionState;       // Transition State
    Molecule* m_ExcessReactant;
    MoleculeManager* m_pMoleculeManager;     // Pointer to molecule manager.
    MicroRateCalculator* m_pMicroRateCalculator; // Pointer to microcanoical rate coeff. calculator.
    TunnelingCalculator* m_pTunnelingCalculator; // Pointer to Tunneling Calculator

    /*
    Each of the backward/forward microcanonical rate coefficients are based on
    the bottom of the relevant well. The cell and grain vectors for each well all have
    the same number of elements; although the number of trailing 0 elements differ
    by the quantity (MaximumGrain - ZpeOfTheWell).
    */

    //
    // Reaction Rate data.
    //

    // _2008_04_24__12_35_40_  <- Please search for this string in the current file for further description.
    double m_FluxCellZPE;              // cell ZPE of m_GrainFlux
    double m_FluxGrainZPE;             // grain ZPE of m_GrainFlux
    size_t m_FluxCellOffset;           // cell Offset when converting m_CellFlux to m_GrainFlux

    std::vector<double>  m_CellFlux;  // Microcanonical transition state fluxes. (QM or classical)
    std::vector<double>  m_GrainFlux; // Grain summed microcanonical transition state fluxes..

    std::vector<double>  m_GrainKfmc; // Grained averaged forward  microcanonical rates.
    std::vector<double>  m_MtxGrnKf;  // Grained averaged forward  microcanonical rates as used in collision operator.

    // Temperature range and interval for rate coefficient test.

    static double m_dTemp;
    static double m_TMin;
    static double m_TMax;

    // Previously private but needed in IrreversibleUnimolecularReaction::calcFluxFirstNonZeroIdx(void)
    bool m_UsesProductProperties;
    int m_GrnFluxFirstNonZeroIdx;  // idx of the starting grain for calculating forward/backward k(E)s from flux
    int m_EffGrainedFwdThreshold;  // effective threshold energy (in grains) for forward flux calculations
    int m_EffGrainedRvsThreshold;  // effective threshold energy (in grains) for backward flux calculations

  private:

    // Grain average microcanonical rate coefficients.
    bool grnAvrgMicroFluxCoeffs();

    const MesmerEnv& m_Env;
    MesmerFlags& m_Flags;
    std::string m_Name;            // Reaction name.

    bool   m_reCalcMicroRateCoeffs; // re-calculation on DOS

    double m_ERConc;       // Concentration of the excess reactant (This is a complement to reactions with
                           // excess species. This value is not used in unimolecular reactions.)
    bool m_bERConcPercent; // If true, excess species concentration is a percentage of bath gas concentration. 
                           // (Also a complement to reactions with excess species, so not used in unimolecular reactions.)

  };

  // _2008_04_24__12_35_40_
  //
  //  Transition state flux construction:
  //
  //  The horizontal dashes on the graph below represents the flux in cell level of a reaction. It can start from
  //  the bottom of the higher well of the reaction as there will be no flux for any energy lower than this point.
  //
  //  It is user's taste to choose where to start a Flux vector, as long as the user specifies the ZPE of the Flux
  //  by setCellFluxBottom()
  //
  //  It is important to calculate the flux so that it is based at least higher than the higher well so that the derivation
  //  of forward/reverse microcanoincal rate coefficients can be processed by Mesmer.
  //
  //                       flux                            kfmc                           kbmc
  //                       ___                             --->                           <---
  //                       ___                             --->                           <---
  //                       ___                             --->                           <---
  //                       ___                             --->             ___           <---
  //                      /___\                            --->            /   \          <---
  //                     / ___ \                           --->           /     \         <---
  //                    /  ___  \                          --->          /       \        <---
  //                   /   ___   \                         --->         /         \       <---
  //                  /    ___    \______ higher well      --->        /           \______<---
  //                 /                                     ---x = 0.0 /
  //                /                                      ---x = 0.0/
  // lower well ___/                                       ---x .___/
  //
  // ------------------
  //
  //                     E        DOS     kf(E)       flux           kr(E)
  //                                  --->             ___           <---
  //                 /  12          6 ---> 10       60 ___           <---
  //         grain  |   11          5 --->  9       45 ___           <---
  //                 \  10          4 --->  8       32 ___           <---
  //                                  --->            /___\          <---
  //                                  --->           / ___ \         <---
  //                                  --->          /  ___  \        <---
  //                                  --->         /   ___   \       <---
  //                                  --->        /    ___    \______<---
  //                                  ---x = 0.0 /
  //                                  ---x = 0.0/
  //                                  ---x .___/
  //
  //  Let's say if we have a grain has three cells (wavenumbers), and their energy are 10, 11, 12, respectively.
  //  Their cell DOS are listed above, which are 4, 5, 6, respectively; also the forward microcanoincal rate coefficients.
  //  Then, we have flux 6*10=60 for cell 12, 5*9=45 for cell 11, 4*8=32 for cell 10.
  //
  //  Conversely, when we first calculating kfmc, we put k(E) = W(E)/ [h * rho(E)]. Flux of cell 12 is simply [W(E)/ h] = 60
  //  without having to divide by rho(E), which is 6 for this cell.

}//namespace
#endif // GUARD_Reaction_h
