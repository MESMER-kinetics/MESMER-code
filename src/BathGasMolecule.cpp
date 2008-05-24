//
// BathGasMolecule.cpp
//
// Author: Struan Robertson
//-------------------------------------------------------------------------------------------
#include "Molecule.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  BathGasMolecule::BathGasMolecule(const MesmerEnv& Env) 
    :Molecule(Env),
    m_Mass(0.0),
    m_Sigma(sigmaDefault),
    m_Epsilon(epsilonDefault),
    m_Mass_chk(-1),
    m_Sigma_chk(-1),
    m_Epsilon_chk(-1)
  {}

  bool BathGasMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      cerr << "InitializeMolecule failed for " << getName() << " before constructing BathGasMolecule.";
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      cerr << "Cannot find argument me:MW in " << getName();
      setFlag(true); // later put a function to calculate the molecular weight if the user forgot to provide it.
    }
    else { istringstream idata(txt); double mass(0.); idata >> mass; setMass(mass);}

    txt = ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      cerr << "BathGasMolecule::Cannot find argument me:sigma.";
      setFlag(true);
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt = ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      cerr << "BathGasMolecule::Cannot find argument me:epsilon.";
      setFlag(true);
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    if (getFlag()){
      cerr << "Error(s) while initializing: " << getName();
      return false;
    }

    return true;
  }

  BathGasMolecule::~BathGasMolecule()
  {
    if (m_Mass_chk == 0){
      cinfo << "m_Mass is provided but not used in " << getName() << endl;
    }
    if (m_Sigma_chk == 0){
      cinfo << "m_Sigma is provided but not used in " << getName() << endl;
    }
    if (m_Epsilon_chk == 0){
      cinfo << "m_Epsilon is provided but not used in " << getName() << endl;
    }
  };
  
  void   BathGasMolecule::setMass(double value)           {
    m_Mass = value;
    m_Mass_chk = 0;
  } ;

  double BathGasMolecule::getMass()                       {
    if (m_Mass_chk >= 0){
      ++m_Mass_chk;
      return m_Mass ;
    }
    else{
      cerr << "m_Mass was not defined but requested in " << getName();
      exit(1);
    }
  } ;

  void   BathGasMolecule::setSigma(double value)          {
    m_Sigma = value;
    m_Sigma_chk = 0;
  } ;

  double BathGasMolecule::getSigma()                      {
    if (m_Sigma_chk >= 0){
      ++m_Sigma_chk;
      return m_Sigma ;
    }
    else{
      cerr << "m_Sigma was not defined but requested in " << getName()
               << ". Default value " << sigmaDefault << " is used.\n";
      //exit(1);
      return m_Sigma ;
    }
  } ;

  void   BathGasMolecule::setEpsilon(double value)        {
    m_Epsilon = value;
    m_Epsilon_chk = 0;
  } ;

  double BathGasMolecule::getEpsilon()                    {
    if (m_Epsilon_chk >= 0){
      ++m_Epsilon_chk;
      return m_Epsilon ;
    }
    else{
      cerr << "m_Epsilon was not defined but requested in " << getName()
               << ". Default value " << epsilonDefault << " is used.\n";
      //exit(1);
      return m_Epsilon ;
    }
  } ;

}//namespace


