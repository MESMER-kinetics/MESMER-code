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
      m_PreExp(0.0), m_EInf(0.0), m_isRvsILTpara(false)  {
      Register();
    }

    virtual const char* getID()  { return m_id; }
    virtual ~SimpleILT() {}
    virtual SimpleILT* Clone() { return new SimpleILT(*this); }

    virtual bool calculateMicroCnlFlux(Reaction* pReac);

    virtual double get_ThresholdEnergy(Reaction* pReac);
    virtual bool ParseData(PersistPtr pp);

  private:
    const char* m_id;

    // All the parameters that follow are for an Arrhenius expression of the type:
    // k(T) = Ainf * exp(-Einf/(RT))

    Rdouble m_PreExp;      // Preexponetial factor
    Rdouble m_EInf;        // E infinity
    bool    m_isRvsILTpara; // The ILT parameters provided are for reverse direction.

  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  SimpleILT theSimpleILT("SimpleILT");
  //************************************************************

  // Read ILT parameters
  bool SimpleILT::ParseData(PersistPtr pp)
  {
    Reaction* pReact = m_parent; //use old var name
    PersistPtr ppReac = pReact->get_PersistentPointer();

    // OpenBabel outputs <rateParameters> <A> <n> <E>
    // Attempt to read these first and if not present read the mesmer 
    // version which will add the default if necessary.

    PersistPtr ppActEne, ppPreExponential;//@, pp;
    const char* pActEnetxt = NULL, *pPreExptxt = NULL;
    bool rangeSet(false);
    PersistPtr ppRateParams = ppReac->XmlMoveTo("rateParameters");
    if (ppRateParams) {
      //OpenBabel form
      ppActEne = ppRateParams->XmlMoveTo("E");
      pActEnetxt = ppRateParams->XmlReadValue("E", optional);
      ppPreExponential = ppRateParams->XmlMoveTo("A");
      pPreExptxt = ppRateParams->XmlReadValue("A");
    }
    else {
      //Mesmer forms
      //New form has <me:MCRCMethod name="SimpleILT" xsi:type="SimpleILT">
      //and parameters as subelement of this. In old form they are siblings.
      //@Now handled in ParseForPlugin
      //@pp = ppReac->XmlMoveTo("me:MCRCMethod");
      //@if(pp && !pp->XmlReadValue("xsi:type", optional))
      //@pp = ppReac;
      ppActEne = pp->XmlMoveTo("me:activationEnergy");
      pActEnetxt = pp->XmlReadValue("me:activationEnergy");
      ppPreExponential = pp->XmlMoveTo("me:preExponential");
      pPreExptxt = pp->XmlReadValue("me:preExponential");
    }

    // Specify the direction of the following ILT parameters.
    m_isRvsILTpara = ppActEne->XmlReadBoolean("reverse");

    // Activation energy details.    
    if (pActEnetxt) {
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput = (unitsTxt) ? unitsTxt : "kJ/mol";
      double value(getConvertedEnergy(unitsInput, tmpvalue));

      if (value < 0.0) {
        cerr << "Activation energy should not be negative when used with ILT." << endl;
        return false;
      }
      ReadRdoubleRange(string(pReact->getName() + ":activationEnergy"), ppActEne, m_EInf,
        rangeSet, getConvertedEnergy(unitsInput, 1.0));
      m_EInf = value;
      if (rangeSet) {
        double valueL, valueU, stepsize;
        m_EInf.get_range(valueL, valueU, stepsize);
        if (valueL < 0.0){
          cerr << "Lower bound of activation energy should not be negative when used with ILT.";
          return false;
        }
      }
    }
    else {
      cerr << "No activation energy specified for ILT method in reaction " << this->getID() << ". Please correct input file.";
      return false;
    }

    if (pPreExptxt)
    {
      double value(0.0);
      stringstream s2(pPreExptxt); s2 >> value;
      ReadRdoubleRange(string(pReact->getName() + ":preExp"), ppPreExponential, m_PreExp, rangeSet);
      m_PreExp = value;
    }
    else {
      cerr << "Specifying ILT without pre-exponential term provided in reaction " << this->getID() << ". Please correct input file.";
      return false;
    }

    return ILTCheck(pReact, ppReac);
  }

  //
  // This method calculates the reaction flux by Laplace inversion
  // of the Arrhenius equation for a reaction proceeding in the 
  // unimolecular direction.
  //

  bool SimpleILT::calculateMicroCnlFlux(Reaction* pReact)
  {
    double relative_ZPE(0.0);
    vector<Molecule *> Isomers;
    if (m_isRvsILTpara) {
      pReact->get_products(Isomers);
      relative_ZPE = pReact->get_relative_pdtZPE();
    }
    else {
      pReact->get_reactants(Isomers);
      relative_ZPE = pReact->get_relative_rctZPE();
    }
    Molecule *pMol = Isomers[0];

    const int MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    // Allocate some work space for and obtain density of states of the unimolecuar reactant.

    vector<double> rctsCellDOS;
    if (!pMol->getDOS().getCellDensityOfStates(rctsCellDOS))
      return false;

    // Calculate microcanonical rate coefficients using simple ILT expression.

    for (size_t i(0); i < MaximumCell; ++i) {
      rxnFlux[i] = m_PreExp * rctsCellDOS[i];
    }

    // The flux bottom energy is equal to the well bottom of the source term
    pReact->setCellFluxBottom(relative_ZPE + m_EInf);

    return true;
  }

  //
  // This function returns the activation energy as the threshold energy. This is not stricitly
  // correct as the activation energy also includes tunnelling effects and temperature dependencies.
  // However, in terms of getting mircocanonical rates it is functionally appropriate.
  //
  double SimpleILT::get_ThresholdEnergy(Reaction* pReac) {

    double RxnHeat = pReac->getHeatOfReaction();

    if (m_EInf < RxnHeat){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }

    return m_EInf;

  }


}//namespace
