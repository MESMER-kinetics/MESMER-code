//
// SuperMolecule.cpp
//
// Author: Struan Robertson
//-------------------------------------------------------------------------------------------
#include "Molecule.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  //
  //Constructor
  //
  SuperMolecule::SuperMolecule(const MesmerEnv& Env) : ModelledMolecule(Env),
    m_mol1(NULL),
    m_mol2(NULL)
  {}

    SuperMolecule::~SuperMolecule()
  {}

  //
  //Initialization
  //
  bool SuperMolecule::InitializeMolecule(PersistPtr pp)
  {
    //the construction of SuperMolecule should always rely on the components, there is therefore no need to initialize it
    //apart from name and m_ppPersist

    if (pp) {
      setPersistentPointer(pp);

      const char* id = pp->XmlReadValue("id");
      if (id) setName(id);
      else{
        cinfo << "Molecular name is absent. Default name <source> is used." << endl;
        string tempName = "source"; setName(tempName);
        //setFlag(true);
      }
      return true;
    }
    else{
      cerr << "Invalid PersistPtr.\n";
      return false;
    }
  }

  int SuperMolecule::getSpinMultiplicity(){
    return (m_mol1->getSpinMultiplicity() * m_mol2->getSpinMultiplicity());
  }

  void SuperMolecule::get_VibFreq(std::vector<double>& vibFreq){
    vibFreq.clear();
    std::vector<double> vmol2;
    m_mol1->get_VibFreq(vibFreq);
    m_mol1->get_VibFreq(vmol2);
    for (std::vector<double>::size_type i = 0; i < vmol2.size(); ++i)
      vibFreq.push_back(vmol2[i]);
  }

  double SuperMolecule::get_zpe() {
    return (m_mol1->get_zpe() + m_mol2->get_zpe());
  }

  DensityOfStatesCalculator* SuperMolecule::get_DensityOfStatesCalculator(){
    return m_mol1->get_DensityOfStatesCalculator();
  }

  // set composing member of the SuperMolecule, also copy necessary properties
  void SuperMolecule::setMembers(CollidingMolecule* mol1p, ModelledMolecule* mol2p){
    m_mol1 = mol1p; std::vector<double> vfMol1; m_mol1->get_VibFreq(vfMol1);
    m_mol2 = mol2p; std::vector<double> vfMol2; m_mol2->get_VibFreq(vfMol2);
  }

  CollidingMolecule* SuperMolecule::getMember1(){
    return m_mol1;
  }

  ModelledMolecule* SuperMolecule::getMember2(){
    return m_mol2;
  }

}//namespace

