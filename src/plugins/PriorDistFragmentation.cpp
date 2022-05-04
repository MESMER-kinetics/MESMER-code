//-------------------------------------------------------------------------------------------
//
// PriorDistFragmentation.cpp
//
// Author: Struan Robertson
// Date:   8th October 2019
//
// This file contains the implemetation of the fragmentation abstract base class and related 
// implementation classes. 
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
  // for the prior model.
  //
  class priorDist : public FragDist
  {
  public:

    ///Constructor which registers with the list of DistributionCalculators in the base class
    priorDist(const char* id) : m_id(id),
      m_pReaction(NULL),
      m_rctDOS(),
      m_upperConv(),
      m_lowerConv(),
      m_WriteDist(false),
      m_XSEne(),
      m_DistTable()
    {
      Register();
    };

    // Destructor.
    virtual ~priorDist() {};

    virtual const char* getID() { return m_id; }
    virtual priorDist* Clone() { return new priorDist(*this); }

    // Read any data from XML and store in this instance. Default is do nothing.
    virtual bool ReadParameters(PersistPtr ppFragDist, std::string name) {

      // Check if fragmentation distribution is to be written.

      PersistPtr pp = ppFragDist->XmlMoveTo("me:WriteDistribution");
      if (pp) {
        m_WriteDist = true;
      }

      return true;
    }

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction);

    // Calculate distribution
    virtual void calculate(double excessEnergy, std::vector<double>& dist);

    // Return resources
    virtual void clear();

  private:

    const char* m_id;

  protected:

    Reaction* m_pReaction;

    vector<double> m_rctDOS;

    vector<double> m_upperConv;

    vector<double> m_lowerConv;

    // Varaibles for holding table data.

    bool m_WriteDist;

    vector<double> m_XSEne;

    vector<vector<double> > m_DistTable;

  };

  //************************************************************
  //Global instance, defining its id
  priorDist thePriorDist("Prior");
  //************************************************************

  // Initialize the fragment distribution.
  void priorDist::initialize(Reaction* pReaction) {

    m_pReaction = pReaction;

    ReactionType reactionType = m_pReaction->getReactionType();

    Molecule* pXsSpcs(NULL), * pSpcs(NULL);
    if (reactionType == PSEUDOISOMERIZATION) {
      pXsSpcs = m_pReaction->getExcessReactant();
      pSpcs = m_pReaction->get_reactant();
    }
    else if (reactionType == BIMOLECULAR_EXCHANGE) {
      vector<Molecule*> products;
      m_pReaction->get_products(products);
      pSpcs = products[0];
      pXsSpcs = products[1];
    }

    vector<double> xsDOS;
    pXsSpcs->getDOS().getCellDensityOfStates(xsDOS);

    pSpcs->getDOS().getCellDensityOfStates(m_rctDOS);

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    size_t Size = xsDOS.size();
    vector<double> Trans_DOS;
    getCellEnergies(Size, pSpcs->getEnv().CellSize, Trans_DOS);
    for (size_t i(0); i < Trans_DOS.size(); i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]);
    }

    FastLaplaceConvolution(xsDOS, Trans_DOS, m_upperConv);

    FastLaplaceConvolution(m_upperConv, m_rctDOS, m_lowerConv);

  };

  // Calculate dissociation distribution

  void priorDist::calculate(double Energy, std::vector<double>& dist) {

    // Get the difference in zero point energies between the well and the adduct.
    const double DeltaH = m_pReaction->getHeatOfReaction();

    // Calcualte threshold for reverse reaction.
    const double rvsThreshold = m_pReaction->get_ThresholdEnergy() - DeltaH;

    // Get the excess energy available for redistribution among bimolecular species. 
    double XsE = max((Energy - rvsThreshold), 0.0);

    const size_t excessEnergy = static_cast<size_t>(XsE);

    const size_t cellOffSet = m_pReaction->getFluxCellOffset();

    if (excessEnergy > 0) {

      // Calculate cell distribution vector. The distribution is shifted by the cell-to-grain
      // off set so as to match the shift applied to the reaction flux above.

      vector<double> mDist(m_rctDOS.size(), 0.0);
      for (size_t i(0), j(cellOffSet); i < excessEnergy && j < mDist.size(); i++, j++) {
        mDist[j] = m_rctDOS[i] * m_upperConv[excessEnergy - i] / m_lowerConv[excessEnergy];
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

    if (m_WriteDist) {
      m_XSEne.push_back(XsE);
      m_DistTable.push_back(dist);
    }

    }

    return;
  }

  // Return resources and write distribution.

  void priorDist::clear() {

    m_rctDOS.clear();
    m_upperConv.clear();
    m_lowerConv.clear();

    // If required, write out the fragmentation distribution.

    if (m_WriteDist && m_XSEne.size() > 0) {
      stest << endl << "Fragmenation distribtuon for reaction " << m_pReaction->getName() << ":" << endl << "{";
      string padding = "               ";
      size_t maxSize(0);
      stest << endl << padding;
      for (size_t i(0); i < m_XSEne.size(); i++) {
        stest << formatFloat(m_XSEne[i], 5, 15);
        maxSize = max(maxSize, m_DistTable[i].size());
      }
      stest << endl << endl;
      stringstream ss;
      size_t forwardThreshE = m_pReaction->get_EffGrnFwdThreshold();
      vector<double> rctEne;
      m_pReaction->get_reactant()->getDOS().getGrainEnergies(rctEne);
      for (size_t j(0), jj(forwardThreshE); j < maxSize; j++, jj++) {
        ss << formatFloat(rctEne[jj], 5, 15);
        for (size_t i(0); i < m_DistTable.size(); i++) {
          vector<double>& tmp = m_DistTable[i];
          if (j < tmp.size())
            ss << formatFloat(tmp[j], 5, 15);
          else
            ss << padding;
        }
        ss << endl;
      }
      stest << ss.str() <<"}" << endl;

      m_XSEne.clear();
      for (size_t i(0); i < m_DistTable.size(); i++)
        m_DistTable[i].clear();
      m_DistTable.clear();
    }
  };


  class modPriorDist : public priorDist
  {
  public:

    ///Constructor which registers with the list of DistributionCalculators in the base class
    modPriorDist(const char* id) : priorDist(id),
      m_order(),
      m_nexp(),
      m_Tref(0.0)
    {};

    // Destructor.
    virtual ~modPriorDist() {};

    virtual modPriorDist* Clone() { return new modPriorDist(*this); }

    // Read any data from XML and store in this instance. Default is do nothing.
    virtual bool ReadParameters(PersistPtr ppFragDist, std::string name) {

      priorDist::ReadParameters(ppFragDist, name);

      m_order = ppFragDist->XmlReadDouble("me:modPriorOrder");
      m_nexp = ppFragDist->XmlReadDouble("me:modPriorNexp");
      m_Tref = ppFragDist->XmlReadDouble("me:modPriorTref");
      bool rangeSet(false);
      PersistPtr ppOrder = ppFragDist->XmlMoveTo("me:modPriorOrder");
      ReadRdoubleRange(name + std::string(":modPriorOrder"), ppOrder, m_order, rangeSet);
      PersistPtr ppNexp = ppFragDist->XmlMoveTo("me:modPriorNexp");
      ReadRdoubleRange(name + std::string(":modPriorNexp"), ppNexp, m_nexp, rangeSet);

      return true;
    };

    // Initialize the fragment distribution.
    virtual void initialize(Reaction* pReaction);

  private:

    Rdouble m_order;
    Rdouble m_nexp;
    double m_Tref;

  };

  //************************************************************
  // Global instance, defining its id
  modPriorDist theModPriorDist("modPrior");
  //************************************************************

  void modPriorDist::initialize(Reaction* pReaction) {

    m_pReaction = pReaction;

    ReactionType reactionType = m_pReaction->getReactionType();

    Molecule* pXsSpcs(NULL), * pSpcs(NULL);
    if (reactionType == PSEUDOISOMERIZATION) {
      pXsSpcs = m_pReaction->getExcessReactant();
      pSpcs = m_pReaction->get_reactant();
    }
    else if (reactionType == BIMOLECULAR_EXCHANGE) {
      vector<Molecule*> products;
      m_pReaction->get_products(products);
      pSpcs = products[0];
      pXsSpcs = products[1];
    }

    vector<double> xsDOS;
    pXsSpcs->getDOS().getCellDensityOfStates(xsDOS);

    pSpcs->getDOS().getCellDensityOfStates(m_rctDOS);

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    size_t Size = xsDOS.size();
    vector<double> Trans_DOS;
    getCellEnergies(Size, pSpcs->getEnv().CellSize, Trans_DOS);
    for (size_t i(0); i < Trans_DOS.size(); i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]);
    }

    const double beta = m_pReaction->getEnv().beta;
    double order = m_order * pow((1.0 / (beta * m_Tref * boltzmann_RCpK)), m_nexp);
    for (size_t i(0); i < m_rctDOS.size(); i++) {
      m_rctDOS[i] = pow(m_rctDOS[i], order);
    }

    FastLaplaceConvolution(xsDOS, Trans_DOS, m_upperConv);

    FastLaplaceConvolution(m_upperConv, m_rctDOS, m_lowerConv);

  };

}