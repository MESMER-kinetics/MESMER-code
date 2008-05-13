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

  bool SinkMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      cerr << "InitializeMolecule failed for " << getName() << " before constructing SinkMolecule";
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt = ppPropList->XmlReadProperty("me:ZPE");
    if(!txt){
      cinfo << "Cannot find argument me:ZPE" << endl;
      m_ZPE_chk = -1;
    }
    else {
      istringstream idata(txt);
      double tempzpe = 0.0;
      idata >> tempzpe;
      txt= ppPropList->XmlReadPropertyAttribute("me:ZPE", "units");
      string unitsInput;
      if (txt){
        unitsInput = txt;
      }
      else{
        unitsInput = "kJ/mol";
      }
      set_zpe(getConvertedEnergy(unitsInput, tempzpe));
      m_ZPE_chk = 0;
    }

    if (getFlag()){
      cerr << "Error(s) while initializing: " << getName();
      return false;
    }

    return true;
  }

  SinkMolecule::~SinkMolecule()
  {
    if (m_ZPE_chk == 0) cinfo << "m_ZPE is provided but not used in " << getName() << endl;
  }

  double SinkMolecule::get_zpe() {
    if (m_ZPE_chk == -1){
      cinfo << "m_ZPE was not defined but requested in " << getName() << ". Default value " << m_ZPE << " is given." << endl;
      --m_ZPE_chk;
      return m_ZPE;
    }
    else if (m_ZPE_chk < -1){
      --m_ZPE_chk;
      return m_ZPE;
    }
    ++m_ZPE_chk;
    return m_ZPE ;
  }

  void SinkMolecule::set_zpe(double value) {
    m_ZPE = value;
    m_ZPE_chk = 0;
  }

}//namespace

