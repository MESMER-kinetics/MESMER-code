//-------------------------------------------------------------------------------------------
//
// EckartCoefficients.h
//
// Author: Chi-Hsiu Liang
// Date:   _2008_03_19__10_02_55_
//
// Produces Eckart tunneling coefficients
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "System.h"
#include "AssociationReaction.h"
#include "IrreversibleExchangeReaction.h"

namespace mesmer
{
  class EckartCoefficients : public TunnelingCalculator
  {
  public:
  
    ///Constructor which registers with the list of TunnelingCalculators in the base class
    EckartCoefficients(const std::string& id) : TunnelingCalculator(id){}
  
    virtual ~EckartCoefficients() {}
  
    virtual bool calculateCellTunnelingCoeffs(Reaction* pReact, std::vector<double>& TunnelingProbability);
  };
}//namespace


