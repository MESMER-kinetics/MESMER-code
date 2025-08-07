//-------------------------------------------------------------------------------------------
//
// GeneralISC.cpp
//
// Authors: Struan Robertson
// Date:    21/Jul/2025
//
// Calculates microcanonical rate coefficients based on an analytical expression.
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
  class GeneralISC : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    GeneralISC(const char* id) : m_format(), m_threshold(0.0), m_id(id) { Register(); }

    virtual ~GeneralISC() {}
    virtual const char* getID() { return m_id; }

    virtual GeneralISC* Clone() { return new GeneralISC(*this); }

    virtual bool ParseData(PersistPtr);

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual double get_ThresholdEnergy(Reaction* pReac) { return m_threshold; };

  private:

    string m_format;

    double m_threshold;
        
    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  GeneralISC theGeneralISC("GeneralISC");
  //************************************************************

  bool GeneralISC::ParseData(PersistPtr pp)
  {
    const char* pFormat = pp->XmlReadValue("Format", optional);

    m_format = pFormat;

    m_threshold = 1700.0;

    return true;
  }

  bool GeneralISC::calculateMicroCnlFlux(Reaction* pReact)
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

    // Calculate microcanonical rate coefficients using GeneralISC expression.
    if (m_format == "Analytic") {
      size_t i(0), j(0);
      double ene = m_threshold;
      j = nint(m_threshold / cellSize);
      size_t upl = min(MaximumCell, size_t((2950.0 - m_threshold)/cellSize));
      for (; i < upl; ++i, ++j, ene += cellSize) {
        rxnFlux[i]  = 2.e03 * exp((ene - m_threshold) / 105.0);
        rxnFlux[i] *= rctCellDOS[j];
      }
      upl = min(MaximumCell, size_t((10000.0 - m_threshold)/ cellSize));
      for (; i < upl; ++i, ++j, ene += cellSize) {
        rxnFlux[i]  = 2.96e08 * exp((ene - 2950.0) / 12300.0);
        rxnFlux[i] *= rctCellDOS[j];
      }
    }
    else {
      // Unknown format.
      throw(std::runtime_error("Unknown format for General ISC."));
    }

    // The flux bottom energy is equal to the ZPE of the transition state.
    pReact->setCellFluxBottom(pReact->get_relative_rctZPE() + pReact->get_ThresholdEnergy());

    return true;
  }

}//namespace
