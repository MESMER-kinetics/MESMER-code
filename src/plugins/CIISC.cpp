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
#include "../MicroRate.h"

using namespace std;
using namespace Constants;
namespace mesmer
{
  class CIISC : public MicroRateCalculator
  {
  public:
    ///Constructor which registers with the list of MicroRateCalculators in the base class
    CIISC(const char* id) : m_CIISCProb(1.0), m_ExpTransDecay(0.0), m_BaseTransDecay(0.0), m_id(id) { Register(); }

    virtual ~CIISC() {}
    virtual const char* getID() { return m_id; }

    virtual CIISC* Clone() { return new CIISC(*this); }

    virtual bool ParseData(PersistPtr);

    virtual bool calculateMicroCnlFlux(Reaction* pReact);

    virtual double get_ThresholdEnergy(Reaction* pReac) { return max(0.0,(pReac->get_relative_pdtZPE() - pReac->get_relative_rctZPE())) ; };

  private:
    
    Rdouble m_CIISCProb; // Collisional transition probability.

    Rdouble m_ExpTransDecay;
    Rdouble m_BaseTransDecay;
    
    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  CIISC theCIISC("CIISC");
  //************************************************************

  bool CIISC::ParseData(PersistPtr pp)
  {
    double CIISCProb = pp->XmlReadDouble("me:CIISCProbability");
    m_CIISCProb = (IsNan(CIISCProb)) ? 1.0 : max(CIISCProb,0.0);

    bool rangeSet;
    PersistPtr ppCIISCProb = pp->XmlMoveTo("me:CIISCProbability");
    ReadRdoubleRange(string(m_parent->getName() + ":CIISCProbability"), ppCIISCProb, m_CIISCProb, rangeSet);

    double ExpTransDecay = pp->XmlReadDouble("me:ExpTransDecay");
    m_ExpTransDecay = (IsNan(ExpTransDecay)) ? 0.0 : max(ExpTransDecay, 0.0);

    PersistPtr ppExpTransDecay = pp->XmlMoveTo("me:ExpTransDecay");
    ReadRdoubleRange(string(m_parent->getName() + ":ExpTransDecay"), ppExpTransDecay, m_ExpTransDecay, rangeSet);

    double BaseTransDecay = pp->XmlReadDouble("me:BaseTransDecay");
    m_BaseTransDecay = (IsNan(BaseTransDecay)) ? 0.0 : max(BaseTransDecay, 0.0);

    PersistPtr ppBaseTransDecay = pp->XmlMoveTo("me:BaseTransDecay");
    ReadRdoubleRange(string(m_parent->getName() + ":BaseTransDecay"), ppBaseTransDecay, m_BaseTransDecay, rangeSet);

    return true;
  }

  bool CIISC::calculateMicroCnlFlux(Reaction* pReact)
  {
    // Allocate some work space for density of states and extract densities of states from reactant.
    vector<double> rctCellDOS;
    if (!pReact->get_reactant()->getDOS().getCellDensityOfStates(rctCellDOS))
      return false;

    // get MaxCell from MesmerEnv structure via Reaction class
    const size_t MaximumCell = pReact->getEnv().MaxCell;
    const double cellSize = pReact->getEnv().CellSize;
    const double coeff = (m_ExpTransDecay > 0.0) ? cellSize / m_ExpTransDecay : 0.0 ;

    // Allocate space to hold transition state flux and initialize elements to zero.
    vector<double>& rxnFlux = pReact->get_CellFlux();
    rxnFlux.clear();
    rxnFlux.resize(MaximumCell, 0.0);

    size_t threshold = size_t(get_ThresholdEnergy(pReact) / cellSize);

    double coll_frq = pReact->get_reactant()->getColl().get_collisionFrequency();
    for (size_t i(0), j(threshold) ; i < MaximumCell; ++i) {
      // Calculate microcanonical rate coefficients using CIISC expression.
      rxnFlux[i] = coll_frq * (m_CIISCProb * exp(-double(i)*coeff) + m_BaseTransDecay);
      rxnFlux[i] *= rctCellDOS[j];
    }

    // The flux bottom energy is equal to the ZPE of the transition state.
    pReact->setCellFluxBottom(max(pReact->get_relative_rctZPE(), pReact->get_relative_pdtZPE()));

    return true;
  }

}//namespace
