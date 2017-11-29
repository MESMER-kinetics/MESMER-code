//-------------------------------------------------------------------------------------------
//
// CanonicalRateCoefficient.cpp
//
// Author: Struan Robertson
// Date:   29/Nov/2017
//
// This file contains the implementation of the CanonicalRateCoefficient class. This class 
// is a dummy MicroRateCalculator that it is used in conjunction with the irreversible
// exchange reaction class to give a loss rate coefficient. This is usually used in the 
// situations where there are competing association and bimolecular reactions that remove 
// the same, deficient, species.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"

using namespace std;

namespace mesmer
{
  class CanonicalRateCoefficient : public MicroRateCalculator
  {
  public:
    CanonicalRateCoefficient(const char* id) : m_id(id),
      m_PreExp(0.0),
      m_NInf(0.0),
      m_TInf(298.0),
      m_EInf(0.0)
    {
      Register();
    }

    virtual const char* getID() { return m_id; }
    virtual ~CanonicalRateCoefficient() {}
    virtual CanonicalRateCoefficient* Clone() { return new CanonicalRateCoefficient(*this); }

    virtual bool calculateMicroCnlFlux(Reaction* pReac);

    virtual double get_ThresholdEnergy(Reaction* pReac);
    virtual bool ParseData(PersistPtr pp);

  private:
    const char* m_id;

    // All the parameters that follow are for an Arrhenius expression of the type:
    // k(T) = Ainf * (T/Tinf)^ninf * exp(-Einf/(RT))

    Rdouble m_PreExp;  // Pre-exponential factor
    Rdouble m_NInf;    // Modified Arrhenius parameter
    double  m_TInf;    // T infinity
    Rdouble m_EInf;    // E infinity

  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  CanonicalRateCoefficient theCanonicalRateCoefficient("CanonicalRateCoefficient");
  //************************************************************

  // Read Arrhenius parameters
  bool CanonicalRateCoefficient::ParseData(PersistPtr pp)
  {
    Reaction* pReact = m_parent; //use old var name
    PersistPtr ppReac = pReact->get_PersistentPointer();

    // Read in Arrhenius parameters if they are define.
    bool rangeSet(false);
    // Activation energy details. 
    PersistPtr ppActEne = pp->XmlMoveTo("me:activationEnergy");
    if (ppActEne) {
      const char* ppActEnetxt = pp->XmlReadValue("me:activationEnergy");
      double tmpvalue = 0.0;
      stringstream s2(ppActEnetxt); s2 >> tmpvalue;
      const char* unitsTxt = ppActEne->XmlReadValue("units", false);
      string unitsInput = (unitsTxt) ? unitsTxt : "kJ/mol";
      double value(getConvertedEnergy(unitsInput, tmpvalue));

      if (value < 0.0)
        throw(std::runtime_error(string("Reaction " + pReact->getName() + " definition has negative activation energy.")));

      ReadRdoubleRange(string(pReact->getName() + ":activationEnergy"), ppActEne, m_EInf, rangeSet, getConvertedEnergy(unitsInput, 1.0));
      m_EInf = value;
      if (rangeSet) {
        double valueL, valueU, stepsize;
        m_EInf.get_range(valueL, valueU, stepsize);

        // Issue an error message and abandon the parsing of the plugin if the activation
        // energy is negative. 
        if (valueL < 0.0)
          throw(std::runtime_error(string("Reaction " + pReact->getName() + " definition has negative lower bound for activation energy.")));
      }
    }
    else
      throw(std::runtime_error(string("No activation energy specified for reaction " + pReact->getName() + ". Please correct input file.")));

    // Pre-exponential factor.
    PersistPtr ppPreExp = pp->XmlMoveTo("me:preExponential");
    if (ppPreExp) {
      const char* ppPreExptxt = pp->XmlReadValue("me:preExponential");
      istringstream s2(ppPreExptxt);
      s2 >> m_PreExp;
      if (m_PreExp < 0.0)
        throw(std::runtime_error(string("Reaction " + pReact->getName() + " definition has negative pre-exponential factor.")));

      ReadRdoubleRange(string(pReact->getName() + ":preExponential"), ppPreExp, m_PreExp, rangeSet);
      if (rangeSet) {
        double valueL, valueU, stepsize;
        m_PreExp.get_range(valueL, valueU, stepsize);
        if (valueL < 0.0)
          throw(std::runtime_error(string("Reaction " + pReact->getName() + " definition has negative lower bound for pre-exponential factor.")));
      }
    }
    else
      throw(std::runtime_error(string("No pre-exponential factor specified for reaction " + pReact->getName() + ". Please correct input file.")));

    // Modified Arrhenius exponent.
    const char* pNInftxt = pp->XmlReadValue("me:nInfinity", optional);
    if (pNInftxt)
    {
      PersistPtr ppNInf = pp->XmlMoveTo("me:nInfinity");
      istringstream s2(pNInftxt); s2 >> m_NInf;
      ReadRdoubleRange(string(pReact->getName() + ":nInfinity"), ppNInf, m_NInf, rangeSet);
    }

    // Reference temperature.
    double TInf = pp->XmlReadDouble("me:TInfinity");
    if (TInf <= 0) {
      cinfo << "Tinfinity is less than or equal to 0; set to the default value of 298 K" << endl;
      m_TInf = 298.0;
    }
    else {
      m_TInf = TInf;
    }

		return true;
  }

  //
  // This method calculates the dummy reaction flux .
  //
  bool CanonicalRateCoefficient::calculateMicroCnlFlux(Reaction* pReact)
  {
    // Allocate space to hold flux.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(1, 0.0);

		const double beta = pReact->getEnv().beta;
    rxnFlux[0] = m_PreExp*pow((1.0 / (beta*boltzmann_RCpK*m_TInf)), m_NInf)*exp(-beta*m_EInf);

    return true;
  }

  //
  // This function returns the activation energy as the threshold energy. This is not stricitly
  // correct as the activation energy also includes tunnelling effects and temperature dependencies.
  // However, in terms of getting mircocanonical rates it is functionally appropriate.
  //
  double CanonicalRateCoefficient::get_ThresholdEnergy(Reaction* pReac) {
    return m_EInf;
  }

}//namespace
