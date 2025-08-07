//-------------------------------------------------------------------------------------------
//
// InternalConversion.cpp
//
// Authors: Struan Robertson
// Date:    21/Jul/2025
//
// Calculates microcanonical rate coeffcients for internal conversion.
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
  class InternalConversion : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    InternalConversion(const char* id) : m_format(), m_id(id) { Register(); }

    virtual ~InternalConversion() {}
    virtual const char* getID() { return m_id; }

    virtual InternalConversion* Clone() { return new InternalConversion(*this); }

    virtual bool ParseData(PersistPtr);

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual double get_ThresholdEnergy(Reaction* pReac) { return max(0.0, (pReac->get_relative_pdtZPE() - pReac->get_relative_rctZPE())); };

  private:

    string m_format;

    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  InternalConversion theInternalConversion("InternalConversion");
  //************************************************************

  bool InternalConversion::ParseData(PersistPtr pp)
  {
    const char* pFormat = pp->XmlReadValue("Format", optional);

    m_format = pFormat;

    return true;
  }

  bool InternalConversion::calculateMicroCnlFlux(Reaction* pReact)
  {
    // get MaxCell from MesmerEnv structure via Reaction class
    const size_t MaximumCell = pReact->getEnv().MaxCell;

    Molecule* p_rct = pReact->get_reactant();

    // Allocate some work space for density of states and extract densities of states from reactant.
    vector<double> rctCellDOS;
    if (!p_rct->getDOS().getCellDensityOfStates(rctCellDOS))
      return false;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    double cellSize = pReact->getEnv().CellSize;
    size_t threshold = size_t(get_ThresholdEnergy(pReact) / cellSize);

    // Calculate microcanonical rate coefficients using InternalConversion expression.
    if (m_format == "Analytic") {
      for (size_t i(0), j(threshold); i < MaximumCell; ++i, ++j) {
        double tmp = 1.11e06 - 450.0 * double(i) * cellSize;
        rxnFlux[i] = (tmp > 0.0) ? tmp : 0.0;
        rxnFlux[i] *= rctCellDOS[j];
      }
    }
    else {
      // Unknown format.
      throw(std::runtime_error("Unknown format for internal conversion."));
    }

    // The flux bottom energy is equal to the ZPE of the transition state.
    pReact->setCellFluxBottom(max(pReact->get_relative_rctZPE(), pReact->get_relative_pdtZPE()));

    return true;
  }

}//namespace
