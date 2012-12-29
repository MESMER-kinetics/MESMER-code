//-------------------------------------------------------------------------------------------
//
// PriorDistribution.cpp
//
// Author: Robin_Shannon_2012_11_26_
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
    PriorDistribution(const char* id) :m_id(id){ Register(); }

    virtual ~PriorDistribution() {}
    virtual const char* getID()  { return m_id; }

    virtual bool calculateDistribution( Molecule* m_host, std::vector<double>& distribution);
	vector<double> CoFragDOS(PersistPtr pp,  Molecule* m_host);
	vector<double> GetNormalisedDist(vector<double> DOS1,  vector<double> DOS2, vector<double> DOSTrans );
    vector<double> GetGrainAveragedDistribution(vector<double> DOS, vector<double> dist,  Molecule* m_host);
    private:
    const char* m_id;
  };

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  PriorDistribution thePriorDistribution("Prior");
  //************************************************************

  bool PriorDistribution::calculateDistribution(Molecule* m_host,
                                                std::vector<double>& dist)
  {

    dist.clear();



    // Get the rovibrational Density of states for the primary species in the prior distribution

	vector<double> DOS1;

	m_host->getDOS().getCellDensityOfStates(DOS1);
   
	// Get the rovibrational Density of states for the Cofragment in the prior distribution

    PersistPtr pp = m_host->get_PersistentPointer();
	
    vector<double> DOS2 = CoFragDOS(pp, m_host); 
	
	// Get the excess energy

    double Xs = pp->XmlReadDouble("me:EnergyExcess");
    const char* p = pp->XmlReadValue("units", optional);
    string units = p ? p : "cm-1";
	double XsE = getConvertedEnergy(units, Xs);
    int XsEne = static_cast<int> (XsE);
	
	// Get average cell energies

	const int MaximumGrain = m_host->getEnv().MaxGrn;
    const int MaximumCell  = m_host->getEnv().MaxCell;
    std::vector<double> Ene;
    getCellEnergies(MaximumCell, Ene);
	dist.resize(MaximumGrain);
	
	//Make sure Excess energy is not larger that the energy of the highest cell.

	if (XsE > MaximumCell) {
	cerr << "Excess energy in prior distribution greater that highest cell energy in master equation";
	XsE = MaximumCell;
	}
	
	// Get the translational density of states.

    vector<double> Trans_DOS(XsEne);

    for(int i=0; i < XsEne; i++) {
        Trans_DOS[i] = 1.0;
    }
    
	//Resize rovibrational DOS vectors so densities so energies greater than the XsEne are not considered

	DOS1.resize(XsEne);
	DOS2.resize(XsEne);

	// Get cell prior distribution for Reactant.

	vector<double> ReactCellDist = GetNormalisedDist(DOS1, DOS2, Trans_DOS);
	
	vector<double> CoReactCellDist = GetNormalisedDist(DOS2, DOS1, Trans_DOS);

	vector<double> TransCellDist = GetNormalisedDist(Trans_DOS, DOS2, DOS1);

	// Print cell distribution if Flag present

    if (m_host->getFlags().InitialDistEnabled){
      ctest << "\nInitial distribution vector" << endl ;
	 ctest << "\nReactant\tCoProduct\tTranslational" << endl ;
	  for (int i=0; i < XsEne; i++){
		formatFloat(ctest, ReactCellDist[i],  6,  15) ;
		formatFloat(ctest, TransCellDist[i],  6,  15) ;
		formatFloat(ctest, CoReactCellDist[i],  6,  15) ;
		ctest << endl ;
	  }
    }

  dist = GetGrainAveragedDistribution(ReactCellDist, dist, m_host);

    return true;
  }

 // Function to get DOS of state in the cofragment species 
 vector<double> PriorDistribution::CoFragDOS(PersistPtr pp,  Molecule* m_host)
  {
     
	PersistPtr pp2 = pp->XmlMoveTo("CoFragment");
    pp2 = pp2->XmlMoveTo("molecule");
    
	MesmerFlags Flag = m_host->getFlags();

	const char* typetxt = "PriorCoFragment";


	 //Construct a new Molecule corresponding to Cofragment 
    Molecule *CoFrag = new Molecule(m_host->getEnv(), Flag, typetxt);
	CoFrag->InitializeMolecule(pp2);
	CoFrag->activateRole(typetxt);
	vector<double> FragDOS;
    CoFrag->getDOS().getCellDensityOfStates(FragDOS);
    return FragDOS;
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

//  // Function to perform convolutions required to obtain the prior distribution
 vector<double> PriorDistribution::GetNormalisedDist(vector<double> DOS1,  vector<double> DOS2, vector<double> DOS3 )
 {
int SizeTemp = DOS1.size();


	vector<double> TempDist(SizeTemp);
    for(int i=0; i<(SizeTemp); i++) {
        for(int j=0; j<(SizeTemp-i); j++) {
            int jj = SizeTemp- i - j - 1;
            TempDist[i] = TempDist[i] + (DOS1[i]) * DOS2[j] * DOS3[jj];
        }
    }

    /* calculate the normalization constant */

    double Norm=0.0;
 for(std::vector<double>::iterator j=TempDist.begin();j!=TempDist.end();++j){
    Norm += *j;
 }

 // Calculate the distribution vector for species with DOS1
 vector<double> Dist(SizeTemp);
 for( int i=0; i<(SizeTemp); i++) {
        for(int j=0; j<(SizeTemp-i); j++) {
           int jj = SizeTemp - i - j - 1 ;
            Dist[i] = Dist[i] + (DOS1[i] / Norm) * DOS2[j] * DOS3[jj];
        }
    }
 
 return Dist;
 }

 vector<double> PriorDistribution::GetGrainAveragedDistribution(vector<double> DOS, vector<double> dist,  Molecule* m_host)
 {
  const int GrainSize = m_host->getEnv().GrainSize;
  int Size = dist.size();
  int Size2 = DOS.size();
  for( int i=0; i<(Size); i++) {       
	   for(int j=0; j<(GrainSize); j++) {
	   int index = i  * GrainSize;
	   if ( index + j < Size2){
		dist[i]+= DOS[index + j];
		}
		else{
		dist[i]+= 0.0;
        }
	   }
  }
  return dist;
 }

}//namespace


