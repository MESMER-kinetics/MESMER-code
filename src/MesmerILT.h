#ifndef GUARD_MesmerILT_h
#define GUARD_MesmerILT_h

#include "System.h"

namespace mesmer
{
  class MesmerILT : public MicroRateCalculator
  {
  public:

    ///Constructor which registers with the list of MicroRateCalculators in the base class
    MesmerILT(const std::string& id) : MicroRateCalculator(id) {}

    virtual ~MesmerILT() {}

    virtual bool calculateMicroRateCoeffs(Reaction* pReact) ;

  private:

    bool calculateAssociationMicroRates(Reaction* pReact);
    bool calculateUnimolecularMicroRates(Reaction* pReact);
    
    bool UnimolecularConvolution(Reaction* pReact);
    bool BimolecularConvolution(Reaction* pReact, vector<double>& ConvolvedCellDOS, double ma, double mb, double mc);

  };
}//namespace

#endif // GUARD_MesmerILT_h
