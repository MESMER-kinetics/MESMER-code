//-------------------------------------------------------------------------------------------
//
// PriorDistribution.cpp
//
// Author: Robin Shannon 
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

    virtual bool calculateDistribution( Molecule* m_host, std::vector<double>& distribution);

  private:

    void CoFragDOS(PersistPtr pp,  Molecule* m_host, vector<double>& DOS) ;

    void GetNormalisedDist(const vector<double>& DOS1, const vector<double>&  DOS2, const vector<double>&  DOS3, vector<double>& Dist) ;

    void GetGrainAveragedDistribution(const vector<double>& DOS, vector<double>& dist,  Molecule* m_host);

    const char* m_id;

  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  PriorDistribution thePriorDistribution("Prior");
  //************************************************************

  bool PriorDistribution::calculateDistribution(Molecule* m_host, std::vector<double>& dist) {

    dist.clear();

    // Get the rovibrational Density of states for the primary species in the prior distribution

    vector<double> DOS1;
    m_host->getDOS().getCellDensityOfStates(DOS1);

    // Get the rovibrational Density of states for the Cofragment in the prior distribution

    PersistPtr pp = m_host->get_PersistentPointer();

    vector<double> DOS2; 
    CoFragDOS(pp, m_host, DOS2); 

    // Get the excess energy

    double Xs = pp->XmlReadDouble("me:EnergyExcess");
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "cm-1";
    int XsEne = static_cast<int> (getConvertedEnergy(units, Xs));

    // Get average cell energies

    const int MaximumGrain = m_host->getEnv().MaxGrn;
    const int MaximumCell  = m_host->getEnv().MaxCell;
    dist.resize(MaximumGrain);

    //Make sure Excess energy is not larger that the energy of the highest cell.

    if (XsEne > MaximumCell) {
      cwarn << "Excess energy in prior distribution greater that highest cell energy in master equation";
      XsEne = MaximumCell;
    }

    // Get the translational density of states.

    vector<double> Trans_DOS(XsEne, 1.0);

    //Resize rovibrational DOS vectors so densities so energies greater than the XsEne are not considered

    DOS1.resize(XsEne);
    DOS2.resize(XsEne);

    // Get cell prior distribution for Reactant.

    vector<double> ReactCellDist(XsEne, 0.0) ;
    GetNormalisedDist(DOS1, DOS2, Trans_DOS, ReactCellDist);

    vector<double> CoReactCellDist(XsEne, 0.0) ;
    GetNormalisedDist(DOS2, DOS1, Trans_DOS, CoReactCellDist);

    vector<double> TransCellDist(XsEne, 0.0) ;
    GetNormalisedDist(Trans_DOS, DOS2, DOS1, TransCellDist);

    // Print cell distribution if Flag present

    if (m_host->getFlags().InitialDistEnabled){
      ctest << "\nInitial distribution vector" << endl ;
      ctest << "\nReactant\tCoProduct\tTranslational" << endl ;
      for (int i=0; i < XsEne; i++){
        formatFloat(ctest, ReactCellDist[i],   6, 15) ;
        formatFloat(ctest, TransCellDist[i],   6, 15) ;
        formatFloat(ctest, CoReactCellDist[i], 6, 15) ;
        ctest << endl ;
      }
    }

    GetGrainAveragedDistribution(ReactCellDist, dist, m_host);

    return true;
  }

  // Function to get DOS of state in the cofragment species 
  void PriorDistribution::CoFragDOS(PersistPtr pp,  Molecule* m_host, vector<double>& DOS) {

    PersistPtr pp2 = pp->XmlMoveTo("CoFragment");
    pp2 = pp2->XmlMoveTo("molecule");

    MesmerFlags Flag = m_host->getFlags();

    const char* typetxt = "PriorCoFragment";

    //Construct a new Molecule corresponding to Cofragment 
    Molecule *CoFrag = new Molecule(m_host->getEnv(), Flag, typetxt);
    CoFrag->InitializeMolecule(pp2);
    CoFrag->activateRole(typetxt);
    CoFrag->getDOS().getCellDensityOfStates(DOS);
  }

  //  Function to perform convolutions required to obtain the prior distribution
  // vector<double> GetNormalisedDist(vector<double> DOS1,  vector<double> DOS2, vector<double> DOSTrans )
  // {
  //int Size = DOS1.size();
  //vector<double> FirstConv(Size);
  //FastLaplaceConvolution(DOS2, DOSTrans, FirstConv);
  //vector<double> seccondConv(Size);
  //FastLaplaceConvolution(DOS1, DOSTrans, seccondConv);
  //
  //vector<double> TripleConv(Size);
  //FastLaplaceConvolution(FirstConv, seccondConv, TripleConv);
  //
  //// Sum elements of triple convolution to get normalisation constant
  //
  // 
  // double Norm=0;
  // for(std::vector<double>::iterator j=TripleConv.begin();j!=TripleConv.end();++j){
  //    Norm += *j;
  // }
  // 
  // // Get Distribution for fragment with DOS1
  // vector<double> DOSNorm(Size);
  // vector<double> Prior(Size);
  //
  // for(int i=1; i<=(Size-1); i++) {
  //	 DOSNorm[i]= (DOS1[i] / Norm);
  // }
  // FastLaplaceConvolution(DOSNorm, FirstConv, Prior);
  // return Prior;
  // }

  //  Function to perform convolutions required to obtain the prior distribution.
  void PriorDistribution::GetNormalisedDist(const vector<double>& DOS1, const vector<double>&  DOS2, const vector<double>&  DOS3,  vector<double>& Dist)
  {
    size_t SizeTemp = DOS1.size();
    double Norm(0.0) ; // Normalization constant.
    for (size_t i(0) ; i < SizeTemp ; i++) {
      for (size_t j(0); j < SizeTemp-i ; j++) {
        size_t jj = SizeTemp- i - j - 1;
        Dist[i] += (DOS1[i]) * DOS2[j] * DOS3[jj];
      }
      Norm += Dist[i] ;
    }

    // Calculate the distribution vector for species with DOS1
    for (size_t i(0) ; i<(SizeTemp); i++) {
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


