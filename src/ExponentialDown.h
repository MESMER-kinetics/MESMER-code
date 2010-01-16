#ifndef GUARD_ExponentialDown_h
#define GUARD_ExponentialDown_h

//-------------------------------------------------------------------------------------------
//
// ExponentialDown.h
//
// Author: Struan Robertson
// Date:   16/Jan/2010
//
// This header file contains the declaration of the ExponentialDown class.
//
// The temperature dependence of <delta_E_down> is accounted for as:
//
// <delta_E_down>(T) = <delta_E_down>_ref * (T / m_DeltaEdownRefTemp)^n
//
//-------------------------------------------------------------------------------------------

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
