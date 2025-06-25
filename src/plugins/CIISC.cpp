//-------------------------------------------------------------------------------------------
//
// CIISC.cpp
//
// Authors: Struan Robertson
// Date:    9/Feb/2025
//
// Calculates microcanonical rate coefficients for a collision induced inter-system crossing.
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "../System.h"
#include "../gWellProperties.h"
#include "../gDensityOfStates.h"
#include "../MicroRate.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  class CIISC : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    CIISC(const char* id) : m_id(id) { Register(); }

    virtual ~CIISC() {}
    virtual const char* getID() { return m_id; }

    virtual CIISC* Clone() { return new CIISC(*this); }

    virtual bool ParseData(PersistPtr);

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual double get_ThresholdEnergy(Reaction* pReac) { return max(0.0,(pReac->get_relative_pdtZPE() - pReac->get_relative_rctZPE())) ; };

  private:
    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  CIISC theCIISC("CIISC");
  //************************************************************

  bool CIISC::ParseData(PersistPtr)
  {
    return true;
  }

  bool CIISC::calculateMicroCnlFlux(Reaction* pReact)
  {
    // get MaxCell from MesmerEnv structure via Reaction class
    const size_t MaximumCell = pReact->getEnv().MaxCell;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    double coll_frq = pReact->get_reactant()->getColl().get_collisionFrequency();
    for (size_t i(0); i < MaximumCell; ++i) {
      // Calculate microcanonical rate coefficients using CIISC expression.
      rxnFlux[i] = 0.01 * coll_frq;
    }

    // The flux bottom energy is equal to the ZPE of the transition state.
    pReact->setCellFluxBottom(max(pReact->get_relative_rctZPE(), pReact->get_relative_pdtZPE()));

    return true;
  }

}//namespace
