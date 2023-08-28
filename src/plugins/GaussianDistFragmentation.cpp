//-------------------------------------------------------------------------------------------
//
// GaussianDistFragmentation.cpp
//
// Author: Robin Shannon
// Date:   16th June 2023
//
// This file contains the implemetation of an abstract class to read a distribution of energies
// post fragmentation from dynamics simulations or some other means
//
//-------------------------------------------------------------------------------------------

#include "../Molecule.h"
#include "../Fragmentation.h"
#include "../gWellProperties.h"
#include "../Reaction.h"

using namespace Constants;
using namespace std;

namespace mesmer
{

  //
  // Implementation class for the calculation of fragment distribution on dissocistion
  // for the a dynamical distribution.
  //

  class GaussianDist : public FragDist
  {
  public:

    ///Constructor which registers with the list of DistributionCalculators in the base class
    GaussianDist(const char* id) : m_id(id),
      m_muIntercept(),
      m_muUnits("cm-1"),
      m_muGradient(0.0),
      m_sigmaIntercept(),
      m_sigmaUnits("cm-1"),
      m_sigmaGradient(0.0),
      m_sigmaTwoIntercept(),
      m_sigmaTwoUnits("cm-1"),
      m_sigmaTwoGradient(0.0),
      m_pReaction(NULL),
      m_rctDOS()
    {
      Register();
    };

    // Destructor.
    virtual ~GaussianDist() {};

    virtual GaussianDist* Clone() { return new GaussianDist(*this); }

    // Read any data from XML and store in this instance. Default is do nothing.
    virtual bool ReadParameters(PersistPtr ppFragDist, std::string name) {
      m_muIntercept = ppFragDist->XmlReadDouble("me:muIntercept");
      // m_muUnits = ppFragDist->XmlReadPropertyAttribute("me:muIntercept", "units", optional);
      m_muGradient = ppFragDist->XmlReadDouble("me:muGradient", optional);
      m_sigmaIntercept = ppFragDist->XmlReadDouble("me:sigmaIntercept");
      // m_sigmaUnits = ppFragDist->XmlReadPropertyAttribute("me:sigmaIntercept", "units", optional);
      m_sigmaGradient = ppFragDist->XmlReadDouble("me:sigmaGradient", optional);
      m_sigmaTwoIntercept = ppFragDist->XmlReadDouble("me:sigmaTwoIntercept");
      // m_sigmaTwoUnits = ppFragDist->XmlReadPropertyAttribute("me:sigmaTwoIntercept", "units", optional);
      m_sigmaTwoGradient = ppFragDist->XmlReadDouble("me:sigmaTwoGradient", optional);
      string out3 = "Reading dynPriorDist params";
      std::cout << out3 << std::endl;
      return true;
    };

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction);

    virtual void calculate(double Energy, std::vector<double>& dist);

    virtual const char* getID() { return m_id; }

    // Return resources
    virtual void clear() {
      m_rctDOS.clear();
    };

  private:
    const char* m_id;
    double m_muIntercept;
    string m_muUnits;
    double m_muGradient;
    double m_sigmaIntercept;
    string m_sigmaUnits;
    double m_sigmaGradient;
    double m_sigmaTwoIntercept;
    string m_sigmaTwoUnits;
    double m_sigmaTwoGradient;

  protected:
    Reaction* m_pReaction;
    vector<double> m_rctDOS;
  };

  //************************************************************
  // Global instance, defining its id
  GaussianDist theGaussianDist("GaussianFrag");
  //************************************************************

  void GaussianDist::initialize(Reaction* pReaction) {

    m_pReaction = pReaction;

    ReactionType reactionType = m_pReaction->getReactionType();
    Molecule* pSpcs(NULL);
    if (reactionType == PSEUDOISOMERIZATION) {
      pSpcs = m_pReaction->get_reactant();
    }
    else if (reactionType == BIMOLECULAR_EXCHANGE) {
      vector<Molecule*> products;
      m_pReaction->get_products(products);
      pSpcs = products[0];
    }

    pSpcs->getDOS().getCellDensityOfStates(m_rctDOS);


  };

  // Calculate dissociation distribution

  void GaussianDist::calculate(double Energy, std::vector<double>& dist) {

    // Get the difference in zero point energies between the well and the adduct.
    const double DeltaH = m_pReaction->getHeatOfReaction();

    // Calcualte threshold for reverse reaction.
    const double rvsThreshold = m_pReaction->get_ThresholdEnergy() - DeltaH;

    // Get the excess energy available for redistribution among bimolecular species. 
    double XsE = max((Energy - rvsThreshold), 0.0);

    const size_t excessEnergy = static_cast<size_t>(XsE);

    const size_t cellOffSet = m_pReaction->getFluxCellOffset();

    // Determine gaussian parameters for the energy distribution at current energy.
    // Currently assumes linear dependence specified by an intercept and gradient. 

    const double tempmu = m_muIntercept + Energy * m_muGradient;
    const double tempsigma = m_sigmaIntercept + Energy * m_sigmaGradient;
    const double tempsigmaTwo = m_sigmaTwoIntercept + Energy * m_sigmaTwoGradient;
    const double m_mu = getConvertedEnergy(m_muUnits, tempmu);
    const double m_sigma = getConvertedEnergy(m_sigmaUnits, tempsigma);
    const double m_sigmaTwo = getConvertedEnergy(m_sigmaTwoUnits, tempsigmaTwo);
    std::cout << m_sigmaTwo << std::endl;

    if (excessEnergy > 0) {

      // Calculate cell distribution vector. The distribution is shifted by the cell-to-grain
      // off set so as to match the shift applied to the reaction flux above.
      // First determine pre-exponetial factor for an asymetric gaussian curve.
      double pre = sqrt(2.0 / M_PI) * (1.0 / (m_sigma + m_sigmaTwo));
      vector<double> mDist(m_rctDOS.size(), 0.0);
      for (size_t i(0), j(cellOffSet); i < excessEnergy && j < mDist.size(); i++, j++) {
        double ene = static_cast<double>(i) / excessEnergy;
        if (ene > m_mu) {
          mDist[j] = pre * exp((-1.0 * pow((ene - m_mu), 2.0)) / (2.0 * pow(m_sigmaTwo, 2.0)));
        }
        else {
          mDist[j] = pre * exp((-1.0 * pow((ene - m_mu), 2.0)) / (2.0 * pow(m_sigma, 2.0)));
        }
      }

      // Average cell distribution over grains.
      // Get grain size for grain averaging.
      const size_t GrainSize = m_pReaction->get_reactant()->getEnv().GrainSize;
      double sum(0.0);
      for (size_t i(0), index(0); i < dist.size(); i++) {
        for (size_t j = 0; j < (GrainSize) && index < mDist.size(); j++, index++) {
          dist[i] += mDist[index];
        }
        sum += dist[i];
      }

      // Normalize fragment distribution.
      double rSum = 1.0 / sum;
      for (size_t i(0); i < dist.size(); i++) {
        dist[i] *= rSum;
      }

    }
    return;
  }

}