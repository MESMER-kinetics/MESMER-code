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


}//namespace

