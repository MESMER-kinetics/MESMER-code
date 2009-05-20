#include "Reaction.h"
//
// Test the forward microcanonical rate coefficients.
//
using namespace std;
using namespace Constants ;

namespace mesmer
{

  bool MicroRateCalculator::testMicroRateCoeffs(
    Reaction*  pReact,
    PersistPtr ppbase) const
  {
    vector<Molecule *> unimolecularspecies;
    pReact->get_unimolecularspecies(unimolecularspecies);
    Molecule * pReactant = unimolecularspecies[0];

    string comment("Microcanonical rate coefficients");
    PersistPtr ppList = ppbase->XmlWriteMainElement("me:microRateList", comment );
    int MaximumGrain = (pReact->getEnv().MaxGrn - pReact->get_fluxFirstNonZeroIdx());

    // Allocate some work space for density of states.

    vector<double> grainEne;
    vector<double> grainDOS;
    const vector<double>& grainKfmc = pReact->get_GrainKfmc();
    pReactant->getDOS().getGrainEnergies(grainEne) ;
    pReactant->getDOS().getGrainDensityOfStates(grainDOS) ;

    ctest << "\nCanonical rate coefficients for " << pReact->getName() << ", calculated from microcanonical rates\n{\n";
    for(int i = 0 ; i < 29 ; ++i)
    {
      double Temperature = double(i+2)*100.0 ;
      double beta = 1.0/(boltzmann_RCpK*Temperature) ;

      double sm1 = 0.0, sm2 = 0.0, tmp = 0.0;

      for ( int i = 0 ; i < MaximumGrain ; ++i ) {
        tmp  = grainDOS[i] * exp(-beta * grainEne[i]) ;
        sm1 += grainKfmc[i] * tmp ;
        sm2 += tmp ;
      }
      sm1 /= sm2 ;
      formatFloat(ctest, Temperature, 6,  7) ;
      formatFloat(ctest, sm1,         6, 15) ;
      ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:microRate");
      ppItem->XmlWriteValueElement("me:T",   Temperature, 6) ;
      ppItem->XmlWriteValueElement("me:val", sm1,         6) ;
    }
    ctest << "}\n";

    return true;
  }
  
  std::string MicroRateCalculator::getName() {return name;}

    // Read ILT parameters
  bool MicroRateCalculator::ReadParameters(Reaction* pReact) {

    pReact->setUsesILT();

    PersistPtr ppReac = pReact->get_PersistentPointer();

    // OpenBabel outputs <rateParameters> <A> <n> <E>
    // Attempt to read these first; 
    // if not present read the mesmer version which will add the default if necessary.
    //**TODO preExponential
    PersistPtr ppActEne, ppPreExponential;
    const char* pActEnetxt=NULL, *pPreExptxt=NULL;
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
      pReact->set_revILT(ppActEne->XmlReadBoolean("reverse")); //specify the direction of the following ILT parameters
      double tmpvalue = 0.0;
      stringstream s2(pActEnetxt); s2 >> tmpvalue ;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput;
      if (unitsTxt){
        unitsInput = unitsTxt;
      }
      else{
        unitsInput = "kJ/mol";
      }
      const char* pLowertxt = ppActEne->XmlReadValue("lower", optional);
      const char* pUppertxt = ppActEne->XmlReadValue("upper", optional);
      const char* pStepStxt = ppActEne->XmlReadValue("stepsize", optional);
      double value(getConvertedEnergy(unitsInput, tmpvalue));
      if (pLowertxt && pUppertxt){
        double tmpvalueL(0.0), tmpvalueU(0.0), tmpstepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt);
        s3 >> tmpvalueL; s4 >> tmpvalueU; s5 >> tmpstepsize;
        double valueL(getConvertedEnergy(unitsInput, tmpvalueL));
        double valueU(getConvertedEnergy(unitsInput, tmpvalueL));
        double stepsize(getConvertedEnergy(unitsInput, tmpstepsize));
        pReact->set_EInf(NaN);
        Rdouble::set_range_indirect(valueL,valueU,stepsize, "EInf");
      }
      else{
        pReact->set_EInf(value);
      }
    }
    else{
      cerr << "Specifying ILT without activation energy provided in reaction "
           << this->getName() << ". Please correct input file.";
      return false;
    }

    if (pPreExptxt)
    {
      double value = 0.0;
      stringstream s2(pPreExptxt); s2 >> value ;
      const char* pLowertxt = ppPreExponential->XmlReadValue("lower", optional);
      const char* pUppertxt = ppPreExponential->XmlReadValue("upper", optional);
      const char* pStepStxt = ppPreExponential->XmlReadValue("stepsize", optional);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt);
        s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        pReact->set_PreExp(NaN);
        Rdouble::set_range_indirect(valueL,valueU,stepsize, "PreExp");

      }
      else{
        pReact->set_PreExp(value);
      }
    }
    else{
      cerr << "Specifying ILT without pre-exponential term provided in reaction " << this->getName() << ". Please correct input file.";
      return false;
    }

    const char* pNInftxt = ppReac->XmlReadValue("me:nInfinity", optional);
    if (pNInftxt)
    {
      PersistPtr ppNInf = ppReac->XmlMoveTo("me:nInfinity") ;
      double value = 0.0;
      stringstream s2(pNInftxt); s2 >> value ;
      const char* pLowertxt = ppNInf->XmlReadValue("lower", false);
      const char* pUppertxt = ppNInf->XmlReadValue("upper", false);
      const char* pStepStxt = ppNInf->XmlReadValue("stepsize", false);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        pReact->set_NInf(NaN);
        Rdouble::set_range_indirect(valueL,valueU,stepsize, "NInf");
      }
      else{
        pReact->set_NInf(value);
      }
    }

    double TInf = ppReac->XmlReadDouble("me:TInfinity");
    if(TInf <= 0) {
      cinfo << "Tinfinity is less than or equal to 0; set to the default value of 298 K" << endl;
      TInf = 298;
    }
    else
      pReact->set_TInf(TInf);         // else set Tinf to the value in the input

    //A few checks on features not allowed in ILT methods
    if (pReact->get_TransitionState())
    {
      cerr << "Reaction " << pReact->getName() 
        << " uses ILT method, which should not have transition state."<<endl;
      return false;
    }
    const char* pTunnelingtxt = ppReac->XmlReadValue("me:tunneling", optional) ;
    if(pTunnelingtxt)
    {
      cerr << "Tunneling parameter in Reaction " << getName() << " is invalid in ILT."<<endl;
      return false;
    }

    return true;
  }

}//namespace
