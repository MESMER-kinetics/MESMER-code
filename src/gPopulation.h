#ifndef GUARD_gPopulation_h
#define GUARD_gPopulation_h

#include "MolecularComponents.h"

namespace mesmer
{
  class gPopulation:public MolecularComponent
  {

    //-------------------------------------------------------------------------------------------------
    // Population and equilibrium fraction
    //-------------------------------------------------------------------------------------------------

  private:

    double m_initPopulation ;   // initial population of the molecule.
    double m_eqFraction ;       // equilibrium fraction of the species
    map <int,double> grainPopulations; // a map which holds any initial grain populations that have been specified

  public:

    //
    // Constructor, destructor and initialization
    //
    gPopulation(Molecule* pMol);

    double getInitPopulation() const { return m_initPopulation;};

    void setInitPopulation(double value) {
      if(grainPopulations.size()==0){
        m_initPopulation = value;
      } // only let the population be specified if there's no grain pop specified
      else{cerr << "initial grain population in this isomer has already been defined... ignoring population specifications " << endl;}
    };

    void getInitGrainPopulation(map<int,double>& inputMap) { 
      inputMap = grainPopulations;
    };

    //  note: any given molecule should have EITHER a total population OR a grain population, BUT NOT BOTH
    void setInitGrainPopulation(int grain, double value) { 
      map<int,double>::iterator it;
      it = grainPopulations.find(grain);
      if(it==grainPopulations.end() && m_initPopulation==0){
        grainPopulations[grain] = value;
        m_initPopulation += value;
      }
      else if(it!=grainPopulations.end()){  // ignore redefinitions of grain populations
        cerr << "initial population of grain " << grain << " has been defined twice... ignoring redefinition " << endl;
      }
      else if(m_initPopulation!=0){ // only let the grain population be specified if there's no total population specified
        cerr << "initial population in this isomer has already been defined... ignoring grain population specifications " << endl;
      }
    };
    double getEqFraction() const { return m_eqFraction;};
    void setEqFraction(double value){ m_eqFraction = value;};

  };
} //namespace
#endif //GUARD_gPopulation_h
