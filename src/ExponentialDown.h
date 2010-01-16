#ifndef GUARD_ExponentialDown_h
#define GUARD_ExponentialDown_h

#include "EnergyTransferModel.h"
namespace mesmer
{
  class ExponentialDown : public EnergyTransferModel
  {
  public:

    ///Constructor which registers with the list of energy transfer models in the base class
    ExponentialDown(const std::string& id) : EnergyTransferModel(id) {}

    virtual ~ExponentialDown() {}

    virtual double calculateTransitionProbability(double Ei, double Ej);

    virtual bool ReadParameters(const gWellProperties* pgWellProperties) ; 

  private:

    double m_deltaEdown ;
    double m_alpha ;

  };
}//namespace

#endif // GUARD_ExponentialDown_h
