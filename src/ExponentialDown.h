#ifndef GUARD_ExponentialDown_h
#define GUARD_ExponentialDown_h

//#include "System.h"
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

    virtual bool ReadParameters(); 

  private:

    double m_deltaEdown ;

  };
}//namespace

#endif // GUARD_ExponentialDown_h
