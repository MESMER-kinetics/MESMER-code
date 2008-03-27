//
// TransitionState.cpp
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
  TransitionState::TransitionState(const MesmerEnv& Env) : ModelledMolecule(Env),
    m_ImFreq(0.0),
    m_ImFreq_chk(-1)
  {}

  TransitionState::~TransitionState()
  {
    if (m_ImFreq_chk == 0) cinfo << "m_ImFreq is provided but not used in " << getName() << endl;
  }

  bool TransitionState::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    PersistPtr oldpp = pp;

    if(!ModelledMolecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule for " << getName() << " before constructing TransitionState with errors.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    bool hasImFreq = false;
    const char* txt;
    txt= ppPropList->XmlReadProperty("me:imFreqs");
    if(!txt){
      cinfo << "No imaginary vibrational frequency for " << getName();
      m_ImFreq_chk = -1;
    }
    else {
      istringstream idata(txt); double x;
      while (idata >> x) m_ImFreq = x;
      hasImFreq = true;
      m_ImFreq_chk = 0;
    }

    return true;
  }

  double TransitionState::get_ImFreq(){
    if (m_ImFreq_chk == -1){
      cerr << "m_ImFreq was not defined but requested in " << getName() << ". Default value " << m_ImFreq << " is given.";
      --m_ImFreq_chk;
      return m_ImFreq;
    }
    else if (m_ImFreq_chk < -1){
      --m_ImFreq_chk;
      return m_ImFreq;
    }
    ++m_ImFreq_chk;
    return m_ImFreq;
  }


}//namespace

