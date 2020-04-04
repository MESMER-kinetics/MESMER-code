//-------------------------------------------------------------------------------------------
//
// BoltzmannDistribution.cpp
//
// Author: Chi-Hsiu Liang
// Date:   _2008_05_15_
//
// Produces Boltzmann distribution
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include "../Distribution.h"
#include "../Molecule.h"
#include "../MesmerFlags.h"
#include "../gDensityOfStates.h"
#include "../gWellProperties.h"

namespace mesmer
{
  class BoltzmannDistribution : public DistributionCalculator
  {
  public:

    ///Constructor which registers with the list of DistributionCalculators in the base class
    BoltzmannDistribution(const char* id) :m_id(id), m_temp(-1.0), m_excitation(0.0) { Register(); }

    virtual ~BoltzmannDistribution() {}
    virtual const char* getID() { return m_id; }
    virtual BoltzmannDistribution* Clone() { return new BoltzmannDistribution(*this); }

    virtual bool ParseData(PersistPtr pp);

    virtual bool calculateDistribution(Molecule* m_host, std::vector<double>& dist);

  private:
    const char* m_id;

    double m_temp;
    double m_excitation;
  };

  //************************************************************
  //Global instance, defining its id
  BoltzmannDistribution theBoltzmannDistribution("Boltzmann");
  //************************************************************

  bool BoltzmannDistribution::ParseData(PersistPtr pp) {

    // Read in in optional temperature (Mostly used for shock tube simulations).

    const char* ptxt = pp->XmlReadValue("me:Temperature", optional);
    if (ptxt) {
      stringstream ss(ptxt);
      ss >> m_temp;
    }

    // Read in in optional excitation (mostly for photo-excitation studies).

    PersistPtr ppExcite = pp->XmlMoveTo("me:Excitation");
    if (ppExcite) {

      const char* p = ppExcite->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";

      const char* ptxt = pp->XmlReadValue("me:Excitation", optional);
      if (ptxt) {
        stringstream ss(ptxt);
        ss >> m_excitation;
      }
      m_excitation = getConvertedEnergy(units, m_excitation);
      if (m_excitation < 0.0) {
        throw std::runtime_error("Initial distribution of" + getParent()->getName() + "specifed with a negative excitation energy.");
      }
    }

    return true;
  }

  bool BoltzmannDistribution::calculateDistribution(Molecule* m_host, std::vector<double>& dist)
  {

    // Note any adjustments need to account for a reservoir are performed by the calling method.

    // Get the grain energies.
    vector<double> Ene;
    m_host->getDOS().getGrainEnergies(Ene);

    //Get the cell energies and rovibrational Densities of States.

    vector<double> cellEne, tmpDOS;
    getCellEnergies(m_host->getEnv().MaxCell, m_host->getEnv().CellSize, cellEne);
    m_host->getDOS().getCellDensityOfStates(tmpDOS);

    // Get the value of beta. 

    const double beta = (m_temp > 0.0) ? 1.0 / (boltzmann_RCpK * m_temp) : m_host->getEnv().beta;

    // Calculate unnormalized Boltzmann dist., accounting for excitation if needed.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    size_t cellOffSet = m_host->getDOS().get_cellOffset();
    size_t nExcitation = size_t(m_excitation) + cellOffSet;
    vector<double> tmp(tmpDOS.size(), 0.0);
    for (size_t i(0), jj(nExcitation); jj < tmp.size(); ++i, ++jj) {
      tmp[jj] = exp(log(tmpDOS[i]) - beta * cellEne[i] + 10.0);
    }

    // Now form the grain sums accounting for excitation and cell off set.
    dist.clear();
    size_t cellPerGrain = m_host->getEnv().cellPerGrain();
    double tsum(0.0);
    for (size_t i(0), ii(0) ; ii< Ene.size() ; ii++ ) {
      double sum(0.0);
      for (size_t j(0) ; j < cellPerGrain && i < tmp.size() ; j++, i++) {
        sum += tmp[i];
      }
      tsum += sum;
      dist.push_back(sum);
    }

    // Normalize distribution. 
    for (size_t i(0); i < dist.size(); ++i) {
      dist[i] /= tsum;
    }

    // Print distribution if requested.
    if (m_host->getFlags().InitialDistEnabled) {
      ctest << "\nInitial distribution vector" << endl;
      for (size_t i(0); i < dist.size(); i++) {
        formatFloat(ctest, dist[i], 6, 15);
        ctest << endl;
      }
    }

    return true;
  }
}//namespace

