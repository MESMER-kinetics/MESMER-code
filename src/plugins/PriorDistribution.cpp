//-------------------------------------------------------------------------------------------
//
// PriorDistribution.cpp
//
// Author: Robin Shannon (Refactored: SHR)
// Date:   26/Nov/2012
//
// Produces Prior Distribution
//
//-------------------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>

#include "../Distribution.h"
#include "../MolecularComponents.h"
#include "../Molecule.h"
#include "../MoleculeManager.h"

namespace mesmer
{
  class PriorDistribution : public DistributionCalculator
  {
  public:

    ///Constructor which registers with the list of DistributionCalculators in the base class
    PriorDistribution(const char* id) : m_id(id){ Register(); }

    virtual ~PriorDistribution() {}
    virtual const char* getID()  { return m_id; }

    virtual bool ParseData(PersistPtr pp);
    virtual bool calculateDistribution( Molecule* m_host, std::vector<double>& distribution);

  private:

    void GetNormalisedDist(const vector<double>& DOS1, const vector<double>&  DOS2, const vector<double>&  DOS3, vector<double>& Dist) ;

    void GetGrainAveragedDistribution(const vector<double>& DOS, vector<double>& dist,  Molecule* m_host);

    const char* m_id;
    Molecule* m_CoFragment;
    int m_XsEne;
  };

  //************************************************************
  //Global instance, defining its id
  PriorDistribution thePriorDistribution("Prior");
  //************************************************************

  bool PriorDistribution::ParseData(PersistPtr pp)
  {
    string CoFragmentName = pp->XmlReadValue("me:CoFragment",optional);
    if (!CoFragmentName.length())
    {
      cerr << "No Co-fragment name specified for Prior distribution." << endl;
      return false;
    }

    double Xs = pp->XmlReadDouble("me:EnergyExcess");
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "cm-1";
    m_XsEne = static_cast<int> (getConvertedEnergy(units, Xs));

    Molecule* host = getParent();
    MesmerFlags& Flags = const_cast<MesmerFlags&>(host->getFlags());
    m_CoFragment = 
      host->getMoleculeManager()->addmol(CoFragmentName, "PriorCoFragment", host->getEnv(), Flags) ;
    m_CoFragment->getDOS();
    return true;
  }

  bool PriorDistribution::calculateDistribution(Molecule* m_host, std::vector<double>& dist) {

    // Get the rovibrational Density of states for the primary species in the prior distribution
    vector<double> DOS1;
    m_host->getDOS().getCellDensityOfStates(DOS1);

    // Get the rovibrational Density of states for the Cofragment in the prior distribution
    vector<double> DOS2; 
    m_CoFragment->getDOS().getCellDensityOfStates(DOS2);

    // Get average cell energies
    const int MaximumGrain = m_host->getEnv().MaxGrn;
    const int MaximumCell  = m_host->getEnv().MaxCell;

    // Make sure Excess energy is not larger that the energy of the highest cell.
    if (m_XsEne > MaximumCell) {
      cwarn << "Excess energy in prior distribution greater that highest cell energy in master equation";
      m_XsEne = MaximumCell;
    }

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    vector<double> Trans_DOS ;
    getCellEnergies(m_XsEne, Trans_DOS) ;
    for (int i(0) ; i < m_XsEne ; i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]) ;
    }

    // Resize rovibrational DOS vectors so densities so energies greater than the XsEne are not considered.

    DOS1.resize(m_XsEne);
    DOS2.resize(m_XsEne);

    // Get cell prior distribution for Reactant.

    vector<double> ReactCellDist(m_XsEne, 0.0) ;
    GetNormalisedDist(DOS1, DOS2, Trans_DOS, ReactCellDist);

    // Print cell distribution if Flag present

    if (m_host->getFlags().InitialDistEnabled){

      vector<double> CoReactCellDist(m_XsEne, 0.0) ;
      GetNormalisedDist(DOS2, DOS1, Trans_DOS, CoReactCellDist);

      vector<double> TransCellDist(m_XsEne, 0.0) ;
      GetNormalisedDist(Trans_DOS, DOS2, DOS1, TransCellDist);

      ctest << "\nInitial distribution vector" << endl ;
      ctest << "\nReactant\tCoProduct\tTranslational" << endl ;
      for (int i=0; i < m_XsEne; i++){
        formatFloat(ctest, ReactCellDist[i],    6, 15) ;
        formatFloat(ctest,  CoReactCellDist[i], 6, 15) ;
        formatFloat(ctest, TransCellDist[i],    6, 15) ;
        ctest << endl ;
      }
    }

    dist.clear();
    dist = vector<double>(MaximumGrain, 0.0) ;
    GetGrainAveragedDistribution(ReactCellDist, dist, m_host);

    return true;
  }

  //  Function to perform convolutions required to obtain the prior distribution
  void PriorDistribution::GetNormalisedDist(const vector<double>& DOS1, const vector<double>&  DOS2, const vector<double>&  DOS3, vector<double>& Dist) {
    size_t Size = DOS1.size();
    vector<double> FirstConv(Size, 0.0);
    FastLaplaceConvolution(DOS2, DOS3, FirstConv);

    // Calculate distribution and normalisation constant.
    double Norm(0.0) ; // Normalization constant.
    for (size_t i(0), j(Size-1) ; i < Size ; i++, j--) {
      Dist[i] = DOS1[i] * FirstConv[j];
      Norm += Dist[i] ;
    }

    // Normalize the distribution vector.
    for (size_t i(0) ; i < Size ; i++) {
      Dist[i] /= Norm ;
    }
  }

  void PriorDistribution::GetGrainAveragedDistribution(const vector<double>& DOS, vector<double>& dist,  Molecule* m_host) {
    const int GrainSize = m_host->getEnv().GrainSize;
    const size_t Size2 = DOS.size();
    size_t index(0) ;
    for (size_t i(0) ; i < dist.size() ; i++) {       
      for (int j=0; j<(GrainSize) && index < Size2; j++, index++ ) {
        dist[i] += DOS[index];
      }
    }
  }

}//namespace


