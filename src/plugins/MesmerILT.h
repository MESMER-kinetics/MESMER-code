#ifndef GUARD_MesmerILT_h
#define GUARD_MesmerILT_h

#include "System.h"

namespace mesmer
{
  class MesmerILT : public MicroRateCalculator
  {
  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    MesmerILT(const std::string& id) : 
      MicroRateCalculator(id),
      m_PreExp(0.0),
      m_NInf(0.0),
      m_TInf(298.0),
      m_EInf(0.0),
      m_isRvsILTpara(false) {}

    virtual ~MesmerILT() {}
    virtual MesmerILT* Clone() { return new MesmerILT(*this); }

    virtual bool calculateMicroRateCoeffs(Reaction* pReact) ;

    virtual double get_ThresholdEnergy(Reaction* pReac) ;
    
    virtual bool ReadParameters(Reaction* pReac); 

  private:

    bool calculateAssociationMicroRates(Reaction* pReact);
    bool calculateUnimolecularMicroRates(Reaction* pReact);
    
    bool UnimolecularConvolution(Reaction* pReact);
    bool BimolecularConvolution(Reaction* pReact, vector<double>& ConvolvedCellDOS, double ma, double mb, double mc);
    
    // All the parameters that follow are for an Arrhenius expression of the type:
    // k(T) = Ainf*(T/Tinf)^ninf * exp(-Einf/(RT))

    Rdouble m_PreExp ;           // Preexponetial factor
    Rdouble m_NInf ;             // Modified Arrhenius parameter
    double  m_TInf ;             // T infinity
    Rdouble m_EInf ;             // E infinity
    bool    m_isRvsILTpara;      // The ILT parameters provided are for reverse direction.

  };
}//namespace

#endif // GUARD_MesmerILT_h
