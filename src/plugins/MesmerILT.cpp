
#include "../AssociationReaction.h"
#include "../gStructure.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class MesmerILT : public MicroRateCalculator
  {
  public:

  /********************************************************************************
  Constructor which registers this class with the list of plugins, initializes the ID
  and also does some initialization specific to this class.
  ********************************************************************************/
    MesmerILT(const char* id) : m_id(id),
      m_PreExp(0.0),
      m_NInf(0.0),
      m_TInf(298.0),
      m_EInf(0.0),
      m_isRvsILTpara(false)
    { Register(); }
  
    virtual const char* getID()  { return m_id; }
    virtual ~MesmerILT() {}

  /******************************************************************************
  Because the class can be used for more than one reaction, a new instance is made
  for each use (in addition to the instance made at startup). This is done by
  ParseForPlugin()(Reaction.cpp near line 250), after finding MesmerILT requested
  in the XML input file, calling Clone(); a function of the following form is
  required in each derived class.
  ******************************************************************************/
    virtual MesmerILT* Clone() { return new MesmerILT(*this); }

  /********************************************************************************
  Optional description which will appear in a verbose listing of plugins.
  Put two spaces after \n  in Description() for good formatting in thelisting.
  ********************************************************************************/
  virtual const char* Description()  { return
    "Uses Inverse Laplace Transform to calculate of microcanonical rates.\n  "
    "Standard ILT, unimolecular ILT, and reverse ILT are implemented.\n  "
    "\n";
  };

  /*********************************************************************
  Called by ParseForPlugin to input the plugin's data from the XML file.
  Can be omitted if the plugins does not need to input its own data.
  *********************************************************************/
  virtual bool ParseData(PersistPtr pp);

  /*********************************************************************
  This is the function which does most of the real work of the plugin.
  See Reaction.cpp around line 112;
  *********************************************************************/
  virtual bool calculateMicroCnlFlux(Reaction* pReact) ;

  virtual double get_ThresholdEnergy(Reaction* pReac) ;

  private:

    bool calculateAssociationMicroRates(Reaction* pReact);
    bool calculateUnimolecularMicroRates(Reaction* pReact);
    
    bool UnimolecularConvolution(Reaction* pReact);
    bool BimolecularConvolution(Reaction* pReact, vector<double>& ConvolvedCellDOS, double ma, double mb) const ;
    
  private:
    const char* m_id;

    // All the parameters that follow are for an Arrhenius expression of the type:
    // k(T) = Ainf*(T/Tinf)^ninf * exp(-Einf/(RT))

    Rdouble m_PreExp ;           // Pre-exponential factor
    Rdouble m_NInf ;             // Modified Arrhenius parameter
    double  m_TInf ;             // T infinity
    Rdouble m_EInf ;             // E infinity
    bool    m_isRvsILTpara;      // The ILT parameters provided are for reverse direction.

  };

  /*********************************************************************************
  Global instance, defining the ID of the plugin.
  IDs should be XML compatible so must not contain any of these characters ><"'&
  **********************************************************************************/
  MesmerILT theMesmerILT("MesmerILT");

  /*********************************************************************
  The ParseData funtion here is more complicated than normal because it
  will handle older alternative formats in addition to the preferred
  "New mesmer form". Here is an example:

  <reaction>
  ...
    <me:MCRCMethod xsi:type="MesmerILT">
      <me:preExponential>6.00e-12</me:preExponential>
      <me:activationEnergy  units="cm-1" reverse="true">546.0</me:activationEnergy>
      <me:TInfinity>1.0</me:TInfinity>
      <me:nInfinity>0.0</me:nInfinity>
    </me:MCRCMethod>

  *********************************************************************/
  bool MesmerILT::ParseData(PersistPtr pp)
  {
    Reaction* pReact = m_parent; //use old var name
    PersistPtr ppReac = pReact->get_PersistentPointer();

    /***************************************************************************
    The input data can be provided in a number of ways (for historical reasons).
    OpenBabel outputs <rateParameters> <A> <n> <E>
    Attempt to read these first, and if not present read the mesmer version.
    ***************************************************************************/
    PersistPtr ppActEne, ppPreExponential;//@, pp;
    const char* pActEnetxt=NULL, *pPreExptxt=NULL;
    bool rangeSet(false) ;
    PersistPtr ppRateParams = ppReac->XmlMoveTo("rateParameters") ;
    if(ppRateParams) {
      //OpenBabel form
      ppActEne = ppRateParams->XmlMoveTo("E") ;
      pActEnetxt = ppRateParams->XmlReadValue("E", optional);
      ppPreExponential = ppRateParams->XmlMoveTo("A") ;
      pPreExptxt = ppRateParams->XmlReadValue("A");
    } 
    else {
      /***********************************************************************************
      Read the Mesmer form of the data
      New form has <me:MCRCMethod name="MesmerILT" xsi:type="me:MesmerILT">
      and parameters as subelement of this. In old form they are siblings.
      Now handled in ParseForPlugin
      pp = ppReac->XmlMoveTo("me:MCRCMethod");
      if(pp && !pp->XmlReadValue("xsi:type", optional))

      The lines below that read the values of"me:activationEnergy" and "me:preExponential"
      will use values taken from defaults.xml if the data XML file does not contain them.
      The appropriate line is added to the internal XML tree, so that when the output XML
      file is used subsequently as input values are explicit. This mechanism, which
      applies for most XmlRead operations unless there is a 'optional' parameter, is the
      recommended way to handle default values. It allows the default to be changed by the
      user, records the use of the default in mesmer.log, and provides error messages,
      including optional exhortations for the user to check the default (see the manual).
      ***********************************************************************************/
      ppActEne = pp->XmlMoveTo("me:activationEnergy") ;
      pActEnetxt = pp->XmlReadValue("me:activationEnergy");
      ppPreExponential = pp->XmlMoveTo("me:preExponential") ;
      pPreExptxt = pp->XmlReadValue("me:preExponential");
    }

    // Specify the direction of the following ILT parameters.
    m_isRvsILTpara = ppActEne->XmlReadBoolean("reverse") ;
   
    // Activation energy details.    
    if (pActEnetxt) {
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput = (unitsTxt) ? unitsTxt : "kJ/mol" ;
      double value(getConvertedEnergy(unitsInput, tmpvalue));
      
      // The output streams cerr, cwarn and cinfo are for progressively less
      // important messages. cerr and cwarn probably appear on the command line
      // and all arewritten to the log file (usually mesmer.log).
      // All messages should have a std::endl.
      if (value<0.0) {
        cerr << "Activation energy should not be negative when used with ILT." << endl;
        return false;
      }

      /******************************************************************************
      m_EInf, the activation energy, behaves most of the time like a normal variable
      of type double.  But it can be a "range variable", taking a range of values
      when used in grid search and fitting routines. When me:activationEnergy hash
      "lower", "upper" and "stepsize" attributes, the following code sets this up.
      Several other data variables can also have this optional feature.
      The first parameter of ReadRdoubleRange is the variable name that will appear
      in the log file to confirm the setup of this feature.
      ******************************************************************************/
      ReadRdoubleRange(string(pReact->getName()+":activationEnergy"), ppActEne, m_EInf,
        rangeSet, getConvertedEnergy(unitsInput, 1.0)) ;
      m_EInf = value ;
      if (rangeSet) {
        double valueL, valueU, stepsize ;
        m_EInf.get_range(valueL,valueU,stepsize) ;

        /******************************************************************************
        Issue an error message and abandon the parsing of the plugin if the activation
        energy is negative. When using cerr or cwarn, the message is displayed on the
        console and also added to the log file (usually mesmer.log). Using cinfo outputs
        to the log file only.
        The message may be prefixed by a "context" like "In R1:". Making an ErrorContext
        sets the current context, and the previous context is restored when the object goes
        out of scope.
        ******************************************************************************/
        if(valueL<0.0){
          cerr << "Lower bound of activation energy should not be negative when used with ILT.";
          return false;
        }
      }
    } else {
      cerr << "No activation energy specified for ILT method in reaction " << this->getID() << ". Please correct input file.";
      return false;
    }

    // Pre-exponential factor details.
    if (pPreExptxt) {
      istringstream s2(pPreExptxt); 
      s2 >> m_PreExp ;
      if (m_PreExp<0.0) {
        cerr << "Pre-exponential factor should not be negative when used with ILT." << endl;
        return false;
      }
      ReadRdoubleRange(string(pReact->getName()+":preExponential"), ppPreExponential, m_PreExp, rangeSet) ;  
      if (rangeSet) {
        double valueL, valueU, stepsize ;
        m_PreExp.get_range(valueL,valueU,stepsize) ;
        if(valueL<0.0){
          cerr << "Lower bound of pre-exponential factor should not be negative when used with ILT.";
          return false;
        }
      }
    } else {
      cerr << "No pre-exponential factor specified for ILT method in reaction " << this->getID() << ". Please correct input file.";
      return false;
    }

    const char* pNInftxt = pp->XmlReadValue("me:nInfinity", optional);
    if (pNInftxt)
    {
      PersistPtr ppNInf = pp->XmlMoveTo("me:nInfinity") ;
      istringstream s2(pNInftxt); s2 >> m_NInf ;
      ReadRdoubleRange(string(pReact->getName()+":nInfinity"), ppNInf, m_NInf, rangeSet) ;  
    }

    double TInf = ppReac->XmlReadDouble("me:TInfinity");
    if(TInf <= 0) {
      cinfo << "Tinfinity is less than or equal to 0; set to the default value of 298 K" << endl;
      m_TInf = 298.0 ;
    } else {
      m_TInf = TInf ;   
    }

    if (m_isRvsILTpara)
      // Read product properties and set flag in Reaction to indicate
      // that even irreversible reactions may use these.
      pReact->setUsesProductProperties();

    // Read <me:excessReactantConc> as a child of <me:MCRCMethod>.
    // If it is not found m_ERConc is set to NaN. This causes
    // <me:excessReactantConc> to be read a child of <Reaction> in Reaction.cpp.
    // A default value from defaults.xml will be used if still not found.
    pReact->set_concExcessReactant(pp->XmlReadDouble("me:excessReactantConc",optional));

    return ILTCheck(pReact, ppReac) ; 
  }

  bool MesmerILT::calculateMicroCnlFlux(Reaction* pReact)
  {
    // Check to see what type of reaction we have
    if (pReact->getReactionType() == ISOMERIZATION || 
      pReact->getReactionType() == IRREVERSIBLE_ISOMERIZATION ||
      pReact->getReactionType() == DISSOCIATION){// if it's unimolecular 
        if(!calculateUnimolecularMicroRates(pReact))  // and the microrate calculation is unsuccessful return false
          return false;
    } else {            // if it's not unimolecular
      if(!calculateAssociationMicroRates(pReact))  // and the microrate calculation is unsuccessful return false
        return false;
    }
    return true;
  }

  bool MesmerILT::calculateUnimolecularMicroRates(Reaction* pReact)
  {
    //
    // Need to determine if the supplied Arrhenius parameters are for the association 
    // or dissociation direction and then invoke the appropriate algorithm.
    //
    double relative_ZPE(0.0) ;
    if (m_isRvsILTpara){
      vector<Molecule *> products ; 
      int numberOfProducts = pReact->get_products(products) ;

      if (numberOfProducts != 2) 
        return false ;

      Molecule* p_rct1 = pReact->get_reactant();
      Molecule* p_pdt1 = products[0];
      Molecule* p_pdt2 = products[1];

      const double ma = p_pdt1->getStruc().getMass();
      const double mb = p_pdt2->getStruc().getMass();
 
      // Allocate some work space for density of states and extract densities of states from molecules.
      vector<double> pdtsCellDOS; // Convoluted cell density of states of reactants.
      
      if(!countDimerCellDOS(p_pdt1->getDOS(), p_pdt2->getDOS(), pdtsCellDOS))
        return false;

      BimolecularConvolution(pReact, pdtsCellDOS, ma, mb) ;

      relative_ZPE = pReact->get_relative_pdtZPE();

    } else {

      UnimolecularConvolution(pReact) ;

      relative_ZPE = pReact->get_relative_rctZPE() ;

    }
    pReact->setCellFluxBottom(relative_ZPE + m_EInf);

    return true;
  }

  bool MesmerILT::calculateAssociationMicroRates(Reaction* pReact)
  {
    AssociationReaction *pAssocReaction = dynamic_cast<AssociationReaction*>(pReact) ;
    if(!pAssocReaction){
      cerr << "The MesmerILT method is not available for Irreversible Exchange Reactions"<< endl;
      return false ;
    }

    vector<Molecule *> unimolecularspecies;
    pAssocReaction->get_unimolecularspecies(unimolecularspecies);

    Molecule*  p_pdt1 = unimolecularspecies[0];
    Molecule*  p_rct1 = pAssocReaction->get_pseudoIsomer();
    Molecule*  p_rct2 = pAssocReaction->get_excessReactant();

    const double ma = p_rct1->getStruc().getMass();
    const double mb = p_rct2->getStruc().getMass();
    const double mc = ma + mb ;

    // Allocate some work space for density of states and extract densities of states from molecules.
    vector<double> rctsCellDOS; // Convoluted cell density of states of reactants.

    pAssocReaction->getRctsCellDensityOfStates(rctsCellDOS) ;

    // Perform convolution.

    BimolecularConvolution(pReact, rctsCellDOS, ma, mb) ;

    // the flux bottom energy is equal to the well bottom of the reactant
    pAssocReaction->setCellFluxBottom(pReact->get_relative_rctZPE() + m_EInf);

    cinfo << "Association ILT calculation completed" << endl;

    return true;
  }

  bool MesmerILT::UnimolecularConvolution(Reaction* pReact)
  {
    //
    // Obtain Arrhenius parameters. Note constraint: Ninf >= 0.0
    //
    const double Ninf = m_NInf ;     
    const double Tinf = m_TInf ;
    const double Ainf = m_PreExp ;

    Molecule* p_rct = pReact->get_reactant();
    size_t MaximumCell = pReact->getEnv().MaxCell;

    // Allocate some work space for density of states and extract densities of states from reactant.
    vector<double> rctCellDOS; 
    if(!p_rct->getDOS().getCellDensityOfStates(rctCellDOS))
      return false;

    //
    // Initialize reaction flux vector.
    //
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    const double gammaValue = MesmerGamma(Ninf);
    const double beta0      = 1.0/(boltzmann_RCpK*Tinf);
    const double constant   = Ainf * pow(beta0,Ninf)/gammaValue;

    vector<double> conv;
    if (Ninf > 0.0 ) {
      //
      // The expression held in the elements of the vector work has been altered from the
      // simple mean to analytic integral_x_to_y{E^(Ninf-1)dE}, where x and y are lower and
      // upper energy limits of the cell respectively.
      //
      vector<double> work(MaximumCell, 0.0);
      double cellSize = pReact->getEnv().CellSize ;
      for (size_t i(0); i < MaximumCell; ++i){
		double ene = double(i)*cellSize ;
        work[i] = (pow((ene+cellSize),Ninf)-pow(ene,Ninf))/Ninf;
      }
      FastLaplaceConvolution(work, rctCellDOS, conv);    // FFT convolution replaces the standard convolution
    } else {
		cerr << "nInfinity for unimolecular ILT must be greater than zero... if you want zero, respecify as a small number, e.g., 0.0001" << endl;
		exit(1);
    }

    for (size_t i(0); i < MaximumCell; ++i)
      rxnFlux[i] = constant * conv[i];

    return true;
  }

  bool MesmerILT::BimolecularConvolution(Reaction* pReact, vector<double>& ConvolvedCellDOS, double ma, double mb) const
  {
    //
    // Initialize reaction flux vector.
    //
    size_t MaximumCell = pReact->getEnv().MaxCell;
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    const double gammaValue = MesmerGamma(m_NInf + 1.5);

    // Note electronic degeneracies were already accounted for in DOS calculations.
    // tp_C = 3.24331e+20: defined in Constant.h, constant used in the translational
    // partition function.

    double _ant = m_PreExp * tp_C * pow( ( ma * mb / (ma + mb)), 1.5 ) / gammaValue;
    _ant /= (pow((m_TInf * boltzmann_RCpK), m_NInf));

    //
    // The expression held in the elements of the vector work has been altered from
	// the simple power of the mean value to analytic integral_x_to_y{E^(Ninf+0.5)dE},
	// where x and y are lower and upper energy limits of the cell respectively. Note
    // constraint: Ninf > -1.5
    //
    const double NinfTrans = m_NInf + 1.5; 
    vector<double> work(MaximumCell, 0.0);
    double cellSize = pReact->getEnv().CellSize ;
    for (size_t i(0); i < MaximumCell; ++i){
	  double ene = double(i)*cellSize ;
      work[i] = (pow((ene+cellSize),NinfTrans)-pow(ene,NinfTrans))/NinfTrans ;
    }

    vector<double> conv;
    FastLaplaceConvolution(work, ConvolvedCellDOS, conv);

    for (size_t i(0); i < MaximumCell; ++i)
      rxnFlux[i] = _ant * conv[i];

    return true;
  }

  //
  // This function the activation energy as the threshold energy. This is not stricitly correct as 
  // the activation energy also includes tunnelling effects and temperature dependencies. However,
  // in terms of getting mircocanonical rates it is functionally appropriate.
  //
  double MesmerILT::get_ThresholdEnergy(Reaction* pReac) {

    //if (m_EInf < 0.0) now checked during parsing
    //  cerr << "Providing negative E_infinity in Reaction " << getName() << " is invalid.";

    double RxnHeat = pReac->getHeatOfReaction(); 
    double threshold(0.0) ;

    if (m_isRvsILTpara){
      threshold = (m_EInf > 0.0) ? m_EInf + RxnHeat : RxnHeat ;
    } else {
      threshold = m_EInf ;
    }

    if (threshold < RxnHeat){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }

    return threshold ;

  }

  //-------------------------------------------------
  //   Short note for variables & abbreviations in Mesmer: (REVERSIBLE association reaction)
  //
  //   zpe_react:           zero-point energy of the reactant
  //   zpe_prodt:           zero-point energy of the product
  //   barri_hgt:           Barrier height (equals to theoretical calculated threshold energy)
  //   Einf:           activation energy (experimental value)
  //   TS:                  transition state
  //   PES:                 potential energy surface
  //   A+B:                 molecules A and B
  //   A-B:                 complex formed by association of molecules A and B
  //            |
  //           /|\          TS
  //            |         *****       -\ barri_hgt                        -\
  //  potential |  A+B ***     *      -/               -\         activation\
  //   energy   |  (+)          *                        \          Energy   \
  //            |                *                        \                  /
  //            |                 *                       /                 /
  //            |                  *                     / zpe_react       /
  //            |               /-  **** A-B            /                -/
  //            |   zpe_prodt  /         (-)           /
  //           O|              \-                    -/
  //              ------------------------------------------------------------->
  //                             reaction coordinate
  //  PES
  //
  //   Definition of a REVERSIBLE association reaction in Mesmer:
  //
  //   1. A REVERSIBLE association reaction is going forward when the reaction is going from left to right in this
  //      potential energy surface.
  //   2. A reaction PES can change in different temperature, caused by rotational contribution to the total energy.
  //-------------------------------------------------

}//namespace
