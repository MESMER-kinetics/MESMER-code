//-------------------------------------------------------------------------------------------
//
// System.cpp
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This file contains the implementation of the System class.
//
//-------------------------------------------------------------------------------------------
#include "System.h"


using namespace std ;
using namespace Constants ;

namespace mesmer
{
  System::System(): m_pMoleculeManager(0), m_pReactionManager(0){
      m_pMoleculeManager = new MoleculeManager() ;
      m_pReactionManager = new ReactionManager(m_pMoleculeManager) ;
  }

  System::~System() {
    delete m_pReactionManager;
    delete m_pMoleculeManager;
  }

  //
  // Parse an input data file.
  //
  bool System::parse(PersistPtr ppIOPtr)
  {
    m_ppIOPtr = ppIOPtr;

    //-------------
    //Molecule List // parse this part inside Reaction
    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    m_pMoleculeManager->set_PersistPtr(ppMolList);

    //-------------
    // remove all previous source terms first
//    PersistPtr ppMol = ppMolList->XmlMoveTo("molecule");
//    while (ppMol){
//      string myType = ppMol->XmlReadValue("me:type", false);
//      if (myType == "source"){
//        ppMolList->XmlRemoveChild(ppMol);
//        ppMol = ppMolList;
//      }
//      ppMol = ppMol->XmlMoveTo("molecule");
//    }

    //-------------
    //Model Parameters
    PersistPtr ppParams = ppIOPtr->XmlMoveTo("me:modelParameters");
    if(ppParams)
    {
      const char* txt = ppParams->XmlReadValue("me:grainSize",false);
      if(txt) { istringstream ss(txt); ss >> m_Env.GrainSize; }

      txt = ppParams->XmlReadValue("me:numberOfGrains",false);
      if(txt) { istringstream ss(txt); ss >> m_Env.MaxGrn; }

      txt = ppParams->XmlReadValue("me:maxTemperature",false);
      if(txt) { istringstream ss(txt); ss >> m_Env.MaxT; }

      txt = ppParams->XmlReadValue("me:energyAboveTheTopWell",false);
      if(txt) { istringstream ss(txt); ss >> m_Env.EAboveWell; }

      if(m_Env.GrainSize!=0.0 && m_Env.MaxGrn!=0)
      {
        stringstream errorMsg;
        errorMsg << "No method is provided to specify me:grainSize and me:numberOfGrains.\n"
                 << "me:numberOfGrains has been ignored";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        m_Env.MaxGrn=0;
      }
    }

    //-------------
    //Reaction List
    PersistPtr ppReacList = ppIOPtr->XmlMoveTo("reactionList");
    if(!ppReacList)
    {
      stringstream errorMsg;
      errorMsg << "No reactions have been specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    if(!m_pReactionManager->addreactions(ppReacList, m_Env)) return false;

    //-------------
    //Reaction Conditions
    PersistPtr ppConditions = ppIOPtr->XmlMoveTo("me:conditions");
    if(!ppConditions)
    {
      stringstream errorMsg; errorMsg << "No conditions specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning); return false;
    }
    string Bgtxt = ppConditions->XmlReadValue("me:bathGas");
    if (Bgtxt.empty()){
      stringstream errorMsg; errorMsg << "No bath gas specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning); return false;
    }
    else{
      string molType = "bathGas";
      m_pMoleculeManager->addmol(Bgtxt, molType, ppMolList, m_Env);
      m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;
    }

    if(!ReadRange("me:temperature", Temperatures, ppConditions))   return false;
    if(!ReadRange("me:conc", Concentrations, ppConditions, false)) return false;
    if(!ReadRange("me:pressure", Pressures, ppConditions, false))  return false;

    if(Concentrations.size()==0 && Pressures.size()==0)
    {
      stringstream errorMsg;
      errorMsg << "It is necessary to specify at least one bath gas concentration or pressure";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    else{
      stringstream errorMsg;
      errorMsg << "Number of concentration points: " << Concentrations.size()
               << ", number of pressure points: " << Pressures.size();
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    PersistPtr ppControl = ppIOPtr->XmlMoveTo("me:control");
    if(ppControl)
    {
      m_Env.microRateEnabled = ppControl->XmlReadBoolean("me:testMicroRates");
    }

    return true;
  }

  //
  // Begin calculation.
  //
  void System::calculate()
  {
    TimeCount events; unsigned int timeElapsed =0;

    if(!SetGrainParams())
      return;

    //std::string id;
    //ModelledMolecule* pmol = NULL;

    WriteMetadata();

    /* Needs rethink on how to do looping over T and P
    //Outer loop is temperature
    vector<double>::iterator Titer;
    for(Titer=Temperatures.begin();Titer!=Temperatures.end();++Titer)
    {
      beta = 1.0/(boltzmann_RCpK*(*Titer)) ;

      size_t imax = !Pressures.empty() ? Pressures.size() : Concentrations.size();
      //Inner loop is concentrations, possibly calculated from pressures if these were specified
      // TODO: Get pressure units right. Pressures currently non-functional!
      for(size_t i=0;i<imax;++i)
      {
        m_Env.conc = !Pressures.empty() ? Pressures[i]*beta : Concentrations[i];

        // Build collison matrix for system.
        m_pReactionManager->BuildSystemCollisionOperator(beta, m_Env) ;

        m_pReactionManager->diagCollisionOperator() ;
      }
    }
    */

    m_Env.beta = 1.0 / (boltzmann_RCpK * Temperatures[0]) ; //temporary statements
    double beta = m_Env.beta;
    m_Env.conc = Concentrations[0];

    //---------------
    //About precision
    string precisionMethod;
    switch(precisionTag)
    {
      case varUseDouble:         precisionMethod = "Double";              break;
      case varUseDoubleDouble:   precisionMethod = "Double-double";       break;
      case varUseQuadDouble:     precisionMethod = "Quad-double";         break;
    }

    {
      stringstream errorMsg;
      errorMsg << "Precision: " << precisionMethod;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    //---------------

    // Build collison matrix for system.
    {
      stringstream errorMsg;
      string thisEvent = "Build Collison Matrix";
      errorMsg << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    if (!m_pReactionManager->BuildSystemCollisionOperator(beta, m_Env)){
        stringstream errorMsg;
        errorMsg << "Failed building system collison operator.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    }

    {
      stringstream errorMsg;
      string thisEvent = "Diagonlize Collision Operator";
      errorMsg << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed)  << " -- Time elapsed: " << timeElapsed << " seconds.\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    m_pReactionManager->diagCollisionOperator() ;

    {
      stringstream errorMsg;
      string thisEvent = "Finish Calculation";
      errorMsg << endl << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed)  << " -- Time elapsed: " << timeElapsed << " seconds.\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    /*
    for (size_t i=0; i < m_pReactionManager->size() ; i++) {

    Reaction *reaction = (*m_pReactionManager)[i] ;

    // Work space to hold microcanonical rates.

    vector<double> kmc(m_Env.MaxGrn,0.0) ;

    reaction->get_MicroRateCoeffs(kmc, m_Env) ;

    CollidingMolecule *pmolecule = reaction->m_rct1 ;

    // for (int i(0) ; i < MAXCELL ; i++)
    //     kmc[i] /= omega ;

    pmolecule->diagCollisionOperator() ;

    // Calculate matrix elements

    double kinf = pmolecule->matrixElement(m_Env.MaxGrn-1,m_Env.MaxGrn-1,kmc,m_Env.MaxGrn) ;

    cout << endl ;
    formatFloat(cout, kinf, 6, 15) ;
    }*/

    {
      stringstream errorMsg;
      errorMsg << events;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
  }

  bool System::SetGrainParams()
  {
    /*
    Grain size and number of grain:

      - Either grain size or number of grains can be specified, but not both.

      - Uses the value of grain size in the datafile, if specified.

      - If grain size is not specified but number of grains is, use a grain size to fit the energy range.
        If neither is specified, the grain size is set to 100cm-1 and the number of grains set so that
        the energy range is sufficient.

    Energy Range:

      - The required total energy domain extends from the lowest zero point energy of the lowest molecule
        to 10kT above the highest. <me:maxTemperature> is used, if specified in the data file.
        But, if the input energy of the system is higher, should this be used?
    */

    //Calculate the energy range covering all modelled molecules
    const double BIG = 1e100;
    m_Env.EMax = 0.0; m_Env.EMin = BIG;

    std::string id; ModelledMolecule* pmol = NULL;
    while(m_pMoleculeManager->GetNextMolecule(id, pmol))
    {
      double zpe = pmol->get_zpe() * KcalPerMolToRC ;
      m_Env.EMax = std::max(m_Env.EMax, zpe);
      m_Env.EMin = std::min(m_Env.EMin, zpe);
    }

    if(m_Env.EMin==BIG || m_Env.EMax==0.0)
    {
      stringstream errorMsg;
      errorMsg << "The modelled molecules do not cover an appropriate energy range";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    //m_Env.MaxT gives the option of widening the energy range
    if(m_Env.MaxT <= 0.0){
      m_Env.MaxT = *max_element(Temperatures.begin(), Temperatures.end());
//      stringstream errorMsg;
//      errorMsg << "Maximum Temperature was not set. Reset Maximum Temperature, me:maxTemperature to remove this message.";
//      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    //EAboveWell: Max energy above the highest well.
    m_Env.EMax = m_Env.EMax + m_Env.EAboveWell * boltzmann_RCpK * m_Env.MaxT;
    if(m_Env.GrainSize <= 0.0){
      m_Env.GrainSize = 100.; //default 100cm-1
      stringstream errorMsg;
      errorMsg << "Grain size was invalid. Reset grain size to default: 100";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }

    if(m_Env.MaxGrn <= 0)
      m_Env.MaxGrn = (int)((m_Env.EMax-m_Env.EMin)/m_Env.GrainSize + 0.5);
    else
      m_Env.GrainSize = (m_Env.EMax-m_Env.EMin)/m_Env.MaxGrn;

    m_Env.MaxCell = (int)(m_Env.GrainSize * m_Env.MaxGrn + 0.5);
    {
      stringstream errorMsg;
      errorMsg << "Cell number = " << m_Env.MaxCell << ", Grain number = " << m_Env.MaxGrn;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    return true;
    /*
     //Hardwired
     m_Env.MaxCell=50000;
     m_Env.MaxGrn =500;
     m_Env.GrainSize = 100;
     return true;
    */
  }

  bool System::ReadRange(const string& name, vector<double>& vals, PersistPtr ppbase, bool MustBeThere)
  {
    PersistPtr pp=ppbase;
    for(;;)
    {
      const char* txt;
      pp = pp->XmlMoveTo(name);
      if(pp)
        txt = pp->XmlRead(); //element may have a value
      else //no more elements
        break;
      if(!txt)
        txt = pp->XmlReadValue("initial"); //or use value of "initial" attribute
      if(!txt)
        return false;
      vals.push_back(atof(txt));

      if((txt=pp->XmlReadValue("increment",false)))//optional attribute
      {
        double incr = atof(txt);
        txt = pp->XmlReadValue("final"); //if have "increment" must have "final"
        if(!txt)
          return false;
        for(double val=vals.back()+incr; val<=atof(txt); val+=incr)
          vals.push_back(val);
      }
    }
    if(MustBeThere && vals.size()==0)
    {
      stringstream errorMsg;
      errorMsg << "Must specify at least one value of " << name;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }
    return true;
  }

  void System::WriteMetadata()
  {
    PersistPtr ppList = m_ppIOPtr->XmlWriteMainElement("metadataList", "");
    PersistPtr ppItem = ppList->XmlWriteElement("metadata");
    ppItem->XmlWriteAttribute("name", "dc:creator");
    ppItem->XmlWriteAttribute("content", "Mesmer v0.1");

    ppItem = ppList->XmlWriteElement("metadata");
    ppItem->XmlWriteAttribute("name", "dc:description");
    ppItem->XmlWriteAttribute("content",
    "Calculation of the interaction between collisional energy transfer and chemical reaction"
    " for dissociation, isomerization and association processes");

    ppItem = ppList->XmlWriteElement("metadata");
    ppItem->XmlWriteAttribute("name", "dc:date");

    //----------------------------------------
    TimeCount events;
    string thisEvent, timeString;
    {
      stringstream errorMsg;
      thisEvent = "Write XML attribute";
      timeString = events.setTimeStamp(thisEvent);
      errorMsg << thisEvent << " at " << timeString;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    ppItem->XmlWriteAttribute("content", timeString);
    //----------------------------------------

    //The user's name should be in an environment variable attached to his account (not a System variable)
    const char* author = getenv("MESMER_AUTHOR");
    if(!author)
      author = "unknown";
    ppItem = ppList->XmlWriteElement("metadata");
    ppItem->XmlWriteAttribute("name", "dc:contributor");
    ppItem->XmlWriteAttribute("content", author);
  }

}//namespace
