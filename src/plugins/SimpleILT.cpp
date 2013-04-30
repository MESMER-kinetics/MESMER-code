#include <vector>
#include <string>
#include "../System.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class SimpleILT : public MicroRateCalculator
  {
  public:
    SimpleILT(const char* id) : m_id(id),
      m_PreExp(0.0), m_EInf(0.0)  { Register(); }
  
    virtual const char* getID()  { return m_id; }
    virtual ~SimpleILT() {}
    virtual SimpleILT* Clone() { return new SimpleILT(*this); }

    virtual bool calculateMicroRateCoeffs(Reaction* pReac);

    virtual double get_ThresholdEnergy(Reaction* pReac) ;

    virtual bool ReadParameters(Reaction* pReac) ;
    
  private:
    const char* m_id;

    // All the parameters that follow are for an Arrhenius expression of the type:
    // k(T) = Ainf * exp(-Einf/(RT))

    Rdouble m_PreExp ; // Preexponetial factor
    Rdouble m_EInf ;   // E infinity

  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleILT theSimpleILT("SimpleILT");
  //************************************************************

  bool SimpleILT::ReadParameters(Reaction* pReact) {

    PersistPtr ppReac = pReact->get_PersistentPointer();

    // OpenBabel outputs <rateParameters> <A> <n> <E>
    // Attempt to read these first and if not present read the mesmer 
    // version which will add the default if necessary.
    
    PersistPtr ppActEne, ppPreExponential, pp;
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
      //Mesmer forms
      //New form has <me:MCRCMethod name="MesmerILT" xsi:type="MesmerILT">
      //and parameters as subelement of this. In old form they are siblings.
      pp = ppReac->XmlMoveTo("me:MCRCMethod");
      if(pp && !pp->XmlReadValue("xsi:type",  optional))
        pp = ppReac;
      ppActEne = pp->XmlMoveTo("me:activationEnergy") ;
      pActEnetxt = pp->XmlReadValue("me:activationEnergy");
      ppPreExponential = pp->XmlMoveTo("me:preExponential") ;
      pPreExptxt = pp->XmlReadValue("me:preExponential");
    }

    // Activation energy details.    
    if (pActEnetxt) {
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput = (unitsTxt) ? unitsTxt : "kJ/mol" ;
      double value(getConvertedEnergy(unitsInput, tmpvalue));
      
      if (value<0.0) {
        cerr << "Activation energy should not be negative when used with ILT." << endl;
        return false;
      }
      ReadRdoubleRange(string(pReact->getName()+":activationEnergy"), ppActEne, m_EInf,
        rangeSet, getConvertedEnergy(unitsInput, 1.0)) ;  
      m_EInf = value ;
      if (rangeSet) {
        double valueL, valueU, stepsize ;
        m_EInf.get_range(valueL,valueU,stepsize) ;
        if(valueL<0.0){
          cerr << "Lower bound of activation energy should not be negative when used with ILT.";
          return false;
        }
      }
    } else {
      cerr << "No activation energy specified for ILT method in reaction " << this->getID() << ". Please correct input file.";
      return false;
    }

    if (pPreExptxt)
    {
      double value(0.0) ;
      stringstream s2(pPreExptxt); s2 >> value ;
      ReadRdoubleRange(string(pReact->getName()+":preExp"), ppPreExponential, m_PreExp, rangeSet) ;  
      m_PreExp = value ;
    } else {
      cerr << "Specifying ILT without pre-exponential term provided in reaction " << this->getID() << ". Please correct input file.";
      return false;
    }

    return ILTCheck(pReact, ppReac) ; 
  }

  //
  // This method calculates the reaction flux by Laplace inversion
  // of the Arrhenius equation for a reaction proceeding in the 
  // unimolecular direction.
  //

  bool SimpleILT::calculateMicroRateCoeffs(Reaction* pReact)
  {
    vector<Molecule *> Isomers ;
    pReact->get_unimolecularspecies(Isomers) ;  
    Molecule *p_rcts = Isomers[0] ;

    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    // Allocate some work space for and obtain density of states of the unimolecuar reactant.

    vector<double> rctsCellDOS; 
    if(!p_rcts->getDOS().getCellDensityOfStates(rctsCellDOS))
      return false;

    // Obtain the Arrhenius parameters.

    const int nEinf = int(m_EInf) ; 
    const double preExp = m_PreExp ;

    // Calculate microcanonical rate coefficients using simple ILT expression.

    for (int i = nEinf; i < MaximumCell ; ++i ) {
      rxnFlux[i] = preExp * rctsCellDOS[i-nEinf];
    }

    // the flux bottom energy is equal to the well bottom of the source term
    pReact->setCellFluxBottom(p_rcts->getDOS().get_zpe() - pReact->getEnv().EMin);

    return true;
  }

  //
  // This function the activation energy as the threshold energy. This is not stricitly correct as 
  // the activation energy alos includes tunnelling effects and temperature dependencies. However,
  // in terms of getting mircocanonical rates it is functionally appropriate.
  //
  double SimpleILT::get_ThresholdEnergy(Reaction* pReac) {

    double RxnHeat = pReac->getHeatOfReaction(); 

    if (m_EInf < RxnHeat){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }

    return m_EInf ;

  }
  
    
}//namespace
