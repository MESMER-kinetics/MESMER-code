#include "SimpleILT.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  SimpleILT theSimpleILT("SimpleILT");
  SimpleILT oldSimpleILT("Simple ILT");
  //************************************************************

  bool SimpleILT::ReadParameters(Reaction* pReact) {

    PersistPtr ppReac = pReact->get_PersistentPointer();

    // OpenBabel outputs <rateParameters> <A> <n> <E>
    // Attempt to read these first; 
    // if not present read the mesmer version which will add the default if necessary.
    PersistPtr ppActEne, ppPreExponential;
    const char* pActEnetxt=NULL, *pPreExptxt=NULL;
    bool rangeSet(false) ;
    PersistPtr ppRateParams = ppReac->XmlMoveTo("rateParameters") ;

    if(ppRateParams) {
      ppActEne = ppRateParams->XmlMoveTo("E") ;
      pActEnetxt = ppRateParams->XmlReadValue("E", optional);
      ppPreExponential = ppRateParams->XmlMoveTo("A") ;
      pPreExptxt = ppRateParams->XmlReadValue("A");
    }
    else {
      ppActEne = ppReac->XmlMoveTo("me:activationEnergy") ;
      pActEnetxt = ppReac->XmlReadValue("me:activationEnergy");
      ppPreExponential = ppReac->XmlMoveTo("me:preExponential") ;
      pPreExptxt = ppReac->XmlReadValue("me:preExponential");
    }

    if (pActEnetxt)
    {
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput = (unitsTxt) ? unitsTxt : "kJ/mol" ;

      double value(getConvertedEnergy(unitsInput, tmpvalue));
      if(value<0.0)
      {
        cerr << "activation energy should not be negative when used with ILT" << endl;
        return false;
      }
      ReadRange(string("activationEnergy"), ppPreExponential, m_EInf, value/tmpvalue, rangeSet) ;
      m_EInf = value ;
      if (rangeSet) {
        double valueL, valueU, stepsize ;
        m_EInf.get_range(valueL,valueU,stepsize) ;
        if(valueL<0.0){
          cerr << "lower bound of activation energy should not be negative when used with ILT";
          return false;
        }
      }
    }
    else{
      cerr << "Specifying ILT without activation energy provided in reaction "
           << this->getName() << ". Please correct input file.";
      return false;
    }

    if (pPreExptxt)
    {
      double value(0.0) ;
      stringstream s2(pPreExptxt); s2 >> value ;
      ReadRange(string("me:preExponential"), ppPreExponential, m_PreExp, 1.0, rangeSet) ;
      m_PreExp = value ;
    } else {
      cerr << "Specifying ILT without pre-exponential term provided in reaction " << this->getName() << ". Please correct input file.";
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
  
//-----------------------------------------------------------------------------------------------
//
// ILT Utility methods
//

  //
  // Utility function to read parameter range. 
  //   The parameter cnvrsnFctr applies any conversion factor that is required.
  //   The parameter rangeSet indicates if a range has actually been set for the rdouble variable.
  //
  bool SimpleILT::ReadRange(const string& name, PersistPtr ppbase, Rdouble& rdouble, double cnvrsnFctr, bool& rangeSet)
  {
    PersistPtr pp = ppbase ;

    const char* pLowertxt = pp->XmlReadValue("lower", optional);
    const char* pUppertxt = pp->XmlReadValue("upper", optional);
    const char* pStepStxt = pp->XmlReadValue("stepsize", optional);

    if (pLowertxt && pUppertxt){
      rangeSet = true ;
      double valueL(0.0), valueU(0.0), stepsize(0.0);
      stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt);
      s3 >> valueL; s4 >> valueU; s5 >> stepsize;
      rdouble.set_range(valueL*cnvrsnFctr,valueU*cnvrsnFctr,stepsize*cnvrsnFctr, name.c_str());
      RangeXmlPtrs.push_back(pp);
    } else {
      rangeSet = false ;
    }

    return true;
  }

  //
  // Utility function to check for inconsistencies. 
  //
  bool SimpleILT::ILTCheck(Reaction* pReac, PersistPtr ppReac)
  {
    // A few checks on features not allowed in ILT methods.
    
    if (pReac->get_TransitionState())
    {
      cerr << "Reaction " << pReac->getName() 
        << " uses ILT method, which should not have transition state."<<endl;
      return false;
    }
    const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling", optional) ;
    if(pTunnelingtxt)
    {
      cerr << "Tunneling parameter in Reaction " << pReac->getName() << " is invalid in ILT."<<endl;
      return false;
    }

    const char* pCrossingtxt = ppReac->XmlReadValue("me:crossing", optional) ;
    if(pCrossingtxt)
    {
      cerr << "Crossing parameter in Reaction " << pReac->getName() << " is invalid in ILT."<<endl;
      return false;
    }
    
    return true ;

  }
    
}//namespace
