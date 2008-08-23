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
#include "Tunneling.h"

namespace mesmer
{
  class Reaction
  {
  public:
    //Orders Molecule pointers using the Molecule name.
    //Using the pointer itself was seen to give unpredictable ordering.
    //See Scott Meyers "Effective STL", Item 20
    struct MoleculePtrLess : public binary_function<const Molecule*, const Molecule*, bool>
    {
      bool operator()(const Molecule* mol1, const Molecule* mol2)const
      { return mol1->getName() < mol2->getName(); }
    };
    struct ReactionPtrLess : public binary_function<const Reaction*, const Reaction*, bool>
    {
      bool operator()(const Reaction* r1, const Reaction* r2)const
      { return r1->getName() < r2->getName(); }
    };

    typedef std::map<CollidingMolecule*, int, MoleculePtrLess> isomerMap ;
    typedef std::map<ModelledMolecule* , int, MoleculePtrLess> sourceMap ;

    // Constructors.

    Reaction(MoleculeManager *pMoleculeManager, const MesmerEnv& Env, const char *id);

    // Destructor.
    virtual ~Reaction();

    // Initialize reaction.
    virtual bool InitializeReaction(PersistPtr ppReac) = 0 ;

    // Access microcanoincal rate coeffcients.
    void get_MicroRateCoeffs(std::vector<double> &kmc) ;

    const std::string& getName() const    { return m_Name ; } ;
    double get_PreExp()                   { return m_PreExp.get_value() ; } ;
    void set_PreExp(double value)         { m_PreExp = value;}
    void set_PreExp(double valueL, double valueU, double stepsize){ m_PreExp.set_range(valueL, valueU, stepsize); };
    double get_NInf()                     { return m_NInf.get_value() ; } ;
    void set_NInf(double value)           { m_NInf = value;}
    void set_NInf(double valueL, double valueU, double stepsize)  { m_NInf.set_range(valueL, valueU, stepsize); }
    double getHeatOfReaction() const      {
      const double pdtZPE = get_relative_pdtZPE();
      const double rctZPE = get_relative_rctZPE();
      return pdtZPE - rctZPE;
    };
    int getHeatOfReactionInt() const      { 
      const int pdtZPE = int(get_relative_pdtZPE());
      const int rctZPE = int(get_relative_rctZPE());
      return pdtZPE - rctZPE; 
    };
    const MesmerEnv& getEnv() const       { return m_Env; } ;
    void resetCalcFlag()                  { reCalcDOS = true; };

    // return reactant and product zero-point energy
    virtual double get_relative_rctZPE(void) const = 0;
    virtual double get_relative_pdtZPE(void) const = 0;
    virtual double get_relative_TSZPE(void) const = 0;

    // Get activiation energy
    double get_ThresholdEnergy(void);
    /* This function should be considered as a function to get the activiation energy. 
    In ILT, not the theoretical threshold energy but the experimental activation energy is used. 
    This function returns user defined m_ActivationEnergy, otherwise zero. 
    ILT can be used in all reaction types if necessary. */

    // get products
    virtual int get_products(std::vector<ModelledMolecule *> &product) const = 0;

    // get the reactant, which reacts in a first order or pseudo first order process
    virtual ModelledMolecule *get_reactant(void) const = 0;

    // get reactant collisionoperator size
    virtual int getRctColloptrsize() = 0;

    // Set activiation energy
    void set_ThresholdEnergy(const double value){m_ActivationEnergy = value; };
    void set_ThresholdEnergy(const double value, const double valueL, const double valueU, const double stepsize){
      m_ActivationEnergy.set_range(valueL, valueU, stepsize);
    };

    TransitionState* get_TransitionState() const { return m_TransitionState ; } ;

    // Get unimolecualr species information:
    virtual int get_unimolecularspecies(std::vector<ModelledMolecule *> &unimolecularspecies) const = 0 ;

    // Get the imaginary frequency of the transitions state.
    double get_TSImFreq(void) const {return m_TransitionState->get_ImFreq() ; } ;

    bool thereIsTunnelling (void) const {return (m_pTunnelingCalculator) ? true : false ; } ;

    void calculateCellTunnelingCoeffs(std::vector<double>& TunnelingProbability) {m_pTunnelingCalculator->calculateCellTunnelingCoeffs(this, TunnelingProbability); } ;

    void calculateGrainTunnelingCoeffs(std::vector<double>& TunnelingProbability);

    // calculate TSFlux in grains
    void TSFluxCellToGrain(const std::vector<double>& shiftedTScellFlux);

    // shift transitions state cell flux
    void shiftTScellFlux(std::vector<double>& shiftedTScellFlux);

    // returns the flux in cells for foreign modifications
    std::vector<double>& get_CellFlux(void) {return m_CellTSFlux; };

    // returns the forward grain microcanoincal rate coefficients for foreign modifications
    const std::vector<double>& get_GrainKfmc(void) {return m_GrainKfmc; };

    // returns the backward grain microcanoincal rate coefficients for foreign modifications
    std::vector<double>& get_GrainKbmc(void) {return m_GrainKbmc; };

    // returns the forward grain microcanoincal rate coefficients for foreign modifications
    const std::vector<double>& get_IrreversibleRateCoefficients(void) {return m_GrainKfmc; };

    // get canonical psuedo first order irreversible loss rate coefficient
    virtual double GetCanonicalIrreversibleLossRate(void){return 0.0;};

    // set the bottom energy of m_CellTSFlux
    void setCellFluxBottom(const double energyValue);

    // get the bottom grain ZPE of m_GrainTSFlux
    const int getTSFluxGrnZPE(void){return int(m_FluxGrainZPE);};

    // get the bottom cell offset of m_CellTSFlux
    const int getTSFluxCellOffset(void){return m_FluxCellOffset;};

    // Wrapper function to calculate and grain average microcanoincal rate coeffcients.
    bool calcGrnAvrgMicroRateCoeffs() ;

    // Add reaction terms to collision matrix.
    virtual void AddReactionTerms(qdMatrix *CollOptr, isomerMap &isomermap, const double rMeanOmega) = 0 ;

    // Is reaction equilibrating and therefore contributes
    // to the calculation of equilibrium fractions.
    virtual bool isEquilibratingReaction(double &Keq, ModelledMolecule **rct, ModelledMolecule **pdt) { return false ; } ;

// is the reaction an irreversible reaction
    virtual bool isIrreversible(){return false;};

    // Calculate reaction equilibrium constant.
    virtual double calcEquilibriumConstant() = 0 ;

  protected:

    // Read a molecule name from the XML file and look it up
    // The defaultType is used if there is no me:type attribute
    Molecule* GetMolRef(PersistPtr pp, const char* defaultType = NULL);

    // I/O and control
    PersistPtr           m_ppPersist;            // Conduit for I/O

    //
    // Reaction composition.
    //

    TransitionState     *m_TransitionState;       // TransitionState

    MoleculeManager     *m_pMoleculeManager ;     // Pointer to molecule manager.
    MicroRateCalculator *m_pMicroRateCalculator ; // Pointer to microcanoical rate coeff. calculator.
    TunnelingCalculator *m_pTunnelingCalculator ; // Pointer to Tunneling Calculator

    /*
    Each of the backward/forward microcanonical rate coefficients are based on
    the bottom of the relevant well. The cell and grain vectors for each well all have 
    the same number of elements; although the number of trailing 0 elements differ
    by the quantity (MaximumGrain - ZpeOfTheWell).
    */

    //
    // Reaction Rate data.
    //
    double m_forwardCanonicalRate;
    double m_backwardCanonicalRate;

    // _2008_04_24__12_35_40_  <- Please search for this string in the current file for further description.
    double m_FluxGrainZPE;                       // grain ZPE of m_GrainTSFlux
    int m_FluxCellOffset;                        // cell Offset when converting m_CellTSFlux to m_GrainTSFlux

    std::vector<double>  m_CellTSFlux ;          // Microcanonical transition state fluxes. (QM or classical)
    std::vector<double>  m_GrainTSFlux ;         // Grain summed microcanonical transition state fluxes..

    std::vector<double>  m_GrainKfmc ;           // Grained averaged forward  microcanonical rates.
    std::vector<double>  m_GrainKbmc ;           // Grained averaged backward microcanonical rates.

  private:

    //   Reaction();

    // Copy constructor.
    //   Reaction(const Reaction& reaction) ;

    // Assignment operator.
    //   Reaction& operator=(const Reaction& reaction) ;

    // Grain average microcanonical rate coefficients.
    bool grnAvrgMicroRateCoeffs();

    // Read parameters requires to determine reaction heats and rates.
    virtual bool ReadRateCoeffParameters(PersistPtr ppReac) = 0;

    // Grain averaged microcanonical rate coefficients.
    virtual void calcGrainRateCoeffs() = 0;

    // Test k(T)
    virtual void testRateConstant() = 0;

    const MesmerEnv& m_Env;
    std::string m_Name ;        // Reaction name.

    bool reCalcDOS;             // re-calculation on DOS
    DPoint m_PreExp ;           // Preexponetial factor
    DPoint m_NInf ;             // Modified Arrhenius parameter
    DPoint m_ActivationEnergy;  // Activation Energy
    double m_kfwd ;             // Forward canonical (high pressure) rate coefficient.
  } ;

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
  //                     TS flux                           kfmc                           kbmc
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
  //                     E        DOS     kf(E)     TS flux         kr(E)
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
