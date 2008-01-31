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
    stringstream errorMsg;

    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for " << getName() << " before constructing BathGasMolecule.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt = ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:sigma.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true);
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt = ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:epsilon.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true);
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    if (getFlag()){
      stringstream errorMsg;
      errorMsg << "Error(s) while initializing: " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    return true;
  }


}//namespace

