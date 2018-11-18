//-------------------------------------------------------------------------------------------
//
// HighPressureRateCoefficients.cpp
//
// Author: Struan Robertson
// Date:   18/Nov/2018
//
// This class implements the methods to calculate the high pressure (canonical) rate
// coefficients for a given range of tempertures using the microcanonical rates obtained
// from input data.
//
//-------------------------------------------------------------------------------------------

#include "../System.h"
#include "../calcmethod.h"

namespace mesmer
{
class HighPressureRateCoefficients : public CalcMethod
{
public:
  HighPressureRateCoefficients(const char* id) : m_id(id),
    m_pp(),
    m_TMax(200.0),
    m_TMin(2200.0),
    m_NTpt(20)
  { Register(); }

  virtual ~HighPressureRateCoefficients() {}
  virtual const char* getID()  { return m_id; }
  virtual HighPressureRateCoefficients* Clone() { return new HighPressureRateCoefficients(*this); }

  // Read in data for this method from XML.
  virtual bool ParseData(PersistPtr pp);

  // Function to do the work.
  virtual bool DoCalculation(System* pSys);

private:

  PersistPtr m_pp;
  const char* m_id;

  double m_TMin;
  double m_TMax;
  size_t m_NTpt;

};

////////////////////////////////////////////////
//Global instance
HighPressureRateCoefficients theHighPressureRateCoefficients("HighPressureRateCoefficients");
///////////////////////////////////////////////

bool HighPressureRateCoefficients::ParseData(PersistPtr pp) {
  //save to write XML
  m_pp = pp;

  bool status(true);

  //Read in range parameters, or use values from defaults.xml.
  m_NTpt = pp->XmlReadInteger("me:NumTemp");

  m_TMax = pp->XmlReadDouble("me:MaxTemp");
  m_TMin = pp->XmlReadDouble("me:MinTemp");
  if (m_TMin > m_TMax)
    throw(std::runtime_error("Hig hPressure Rate Coefficients: Max. Temp. less than Min. Temp."));

  return status;
};

bool HighPressureRateCoefficients::DoCalculation(System* pSys) {

  bool status(true);

  // Set some environment variables, even though some of them are not used directly.

  MesmerEnv& Env = pSys->getEnv();

  if (Env.GrainSize <= 0)
    Env.GrainSize = 100;
  Env.MaxGrn = Env.MaxCell / Env.GrainSize;

  bool bMicroRateEnabled = pSys->m_Flags.microRateEnabled;
  pSys->m_Flags.microRateEnabled = true;

  // Loop over reactions

  ReactionManager *pReactionManager = pSys->getReactionManager();

  for (size_t i(0); i < pReactionManager->size(); ++i) {

    Reaction *pReaction = (*pReactionManager)[i];

    pReaction->calcGrnAvrgMicroRateCoeffs();

  }

  pSys->m_Flags.microRateEnabled = bMicroRateEnabled;

  return status;
}

}//namespace

