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
#include "../Molecule.h"
#include "../MoleculeManager.h"
#include "../gDensityOfStates.h"

namespace mesmer
{
  class PriorDistribution : public DistributionCalculator
  {
  public:

    ///Constructor which registers with the list of DistributionCalculators in the base class
    PriorDistribution(const char* id) : m_id(id){ Register(); }

    virtual ~PriorDistribution() {}
    virtual const char* getID()  { return m_id; }
    virtual PriorDistribution* Clone() { return new PriorDistribution(*this); }

    virtual bool ParseData(PersistPtr pp);
    virtual bool calculateDistribution( Molecule* m_host, std::vector<double>& distribution);

  private:

    void GetNormalisedDist(const vector<double>& DOS1, const vector<double>&  DOS2, const vector<double>&  DOS3, vector<double>& Dist) ;

    void GetGrainAveragedDistribution(const vector<double>& DOS, vector<double>& dist,  Molecule* m_host);

    const char* m_id;
    Molecule* m_CoFragment;
    Rdouble m_XsEne;
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

    PersistPtr ppXsEne = pp->XmlMoveTo("me:EnergyExcess");
    if(!ppXsEne)
    {
      cerr << "No Co-fragment excessEnergy specified in Prior distribution." << endl;
      return false;
    }
    const char* p = ppXsEne->XmlReadValue("units", optional);
    string units = p ? p : "cm-1";
    m_XsEne = getConvertedEnergy(units, pp->XmlReadDouble("me:EnergyExcess"));
    bool rangeSet;
    ReadRdoubleRange("Prior:ExcessEnergy", ppXsEne, m_XsEne, rangeSet, getConvertedEnergy(units, 1.0));

    Molecule* host = getParent();
    MesmerFlags& Flags = const_cast<MesmerFlags&>(host->getFlags());
    m_CoFragment = 
      host->getMoleculeManager()->addmol(CoFragmentName, "PriorCoFragment", host->getEnv(), Flags) ;
//    m_CoFragment->getDOS();

    return true;
  }

  bool PriorDistribution::calculateDistribution(Molecule* m_host, std::vector<double>& dist) {
    
    size_t XsEne = static_cast<size_t>(m_XsEne);

    // Get the rovibrational Density of states for the primary species in the prior distribution
    vector<double> DOS1;
    m_host->getDOS().getCellDensityOfStates(DOS1);

    // Get the rovibrational Density of states for the Cofragment in the prior distribution
    vector<double> DOS2; 
    m_CoFragment->getDOS().getCellDensityOfStates(DOS2);

    // Get average cell energies
    const size_t MaximumGrain = m_host->getEnv().MaxGrn;
    const size_t MaximumCell  = m_host->getEnv().MaxCell;

    // Make sure Excess energy is not larger that the energy of the highest cell.
    if (XsEne > MaximumCell) {
      cwarn << "Excess energy in prior distribution greater that highest cell energy in master equation";
      XsEne = MaximumCell;
    }

    // The (classical) translational density of states. Prefactors are not included 
    // because they cancel on normalization.

    vector<double> ene, Trans_DOS ;
    getCellEnergies(XsEne, m_host->getEnv().CellSize, ene) ;
    Trans_DOS = ene;
    for (size_t i(0) ; i < Trans_DOS.size() ; i++) {
      Trans_DOS[i] = sqrt(Trans_DOS[i]) ;
    }

    // Resize rovibrational DOS vectors so densities and energies greater than the XsEne are not considered.

    DOS1.resize(XsEne);
    DOS2.resize(XsEne);

    // Get cell prior distribution for Reactant.

    vector<double> ReactCellDist(XsEne, 0.0) ;
    GetNormalisedDist(DOS1, DOS2, Trans_DOS, ReactCellDist);

    // Print cell distribution if Flag present

    if (m_host->getFlags().InitialDistEnabled){

      vector<double> CoReactCellDist(XsEne, 0.0) ;
      GetNormalisedDist(DOS2, DOS1, Trans_DOS, CoReactCellDist);

      vector<double> TransCellDist(XsEne, 0.0) ;
      GetNormalisedDist(Trans_DOS, DOS2, DOS1, TransCellDist);

      ctest << endl << "Initial distribution vector" << endl ;
      ctest << endl << setw(16) << "Energy/cm-1," << setw(17) << "Reactant," << setw(17) << "CoProduct," << setw(16) << "Translational" << endl ;
      for (int i=0; i < XsEne; i++){
        ctest << formatFloat(ene[i],            6, 15) << ", ";
        ctest << formatFloat(ReactCellDist[i],  6, 15) << ", ";
        ctest << formatFloat(CoReactCellDist[i],6, 15) << ", ";
        ctest << formatFloat(TransCellDist[i],  6, 15) << endl ;
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
    const size_t GrainSize = m_host->getEnv().GrainSize;
    const size_t Size2 = DOS.size();
    size_t index(0) ;
    for (size_t i(0) ; i < dist.size() ; i++) {       
      for (size_t j(0) ; j < GrainSize && index < Size2; j++, index++ ) {
        dist[i] += DOS[index];
      }
    }
  }

}//namespace


