#ifndef GUARD_SimpleILT_h
#define GUARD_SimpleILT_h

#include <vector>
#include <string>
#include "../System.h"

namespace mesmer
{
  class SimpleILT : public MicroRateCalculator
  {
  public:

    // Constructor which registers with the list of MicroRateCalculators in the base class.
    SimpleILT(const std::string& id) : 
      MicroRateCalculator(id), 
      m_PreExp(0.0),
      m_EInf(0.0) {}

    virtual ~SimpleILT() {}
    virtual SimpleILT* Clone() { return new SimpleILT(*this); }

    virtual bool calculateMicroRateCoeffs(Reaction* pReac);

    virtual double get_ThresholdEnergy(Reaction* pReac) ;

    virtual bool ReadParameters(Reaction* pReac) ;
    
    // Utility function to read parameter range. 
    static bool ReadRange(const string& name, PersistPtr ppbase, Rdouble& rdouble, double cnvrsnFctr, bool& rangeSet) ;

    // Utility function to check for inconsistencies. 
    static bool ILTCheck(Reaction* pReac, PersistPtr ppReac) ;
    
  private:
   
    // All the parameters that follow are for an Arrhenius expression of the type:
    // k(T) = Ainf * exp(-Einf/(RT))

    Rdouble m_PreExp ; // Preexponetial factor
    Rdouble m_EInf ;   // E infinity

  };
}//namespace

#endif // GUARD_SimpleILT_h
