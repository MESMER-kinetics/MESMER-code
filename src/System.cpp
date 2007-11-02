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

  System::~System() { }

  //
  // Parse an input data file.
  //
  bool System::parse(PersistPtr ppIOPtr)
  {
    //-------------
    //Molecule List
    m_ppIOPtr = ppIOPtr;
    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    if(!ppMolList)
    {
      stringstream errorMsg;
      errorMsg << "No molecules have been specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    if(!m_pMoleculeManager->addmols(ppMolList))
      return false;
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
    if(!m_pReactionManager->addreactions(ppReacList))
      return false;

    //-------------
    //Reaction Conditions
    PersistPtr ppConditions = ppIOPtr->XmlMoveTo("me:conditions");
    if(!ppConditions)
    {
      stringstream errorMsg;
      errorMsg << "No conditions have been specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    const string Bgtxt = ppConditions->XmlReadValue("me:bathGas");
    if(!Bgtxt.size() || !(m_pMoleculeManager->find(Bgtxt)) )
    {
      stringstream errorMsg;
      errorMsg << "No bath gas specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    } else {
      m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;
    }

    if(!ReadRange("me:temperature", Temperatures, ppConditions))
      return false;

    if(!ReadRange("me:conc", Concentrations, ppConditions, false))
      return false;

    if(!ReadRange("me:pressure", Pressures, ppConditions, false))
      return false;

    if(Concentrations.size()==0 && Pressures.size()==0)
    {
      stringstream errorMsg;
      errorMsg << "It is necessary to specify at least one bath gas concentration or pressure";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    //-------------
    //Model Parameters
    PersistPtr ppParams = ppIOPtr->XmlMoveTo("me:modelParameters");

    if(ppParams)
    {
      const char* txt = ppParams->XmlReadValue("me:grainSize",false);
      if(txt)
      {
        istringstream ss(txt);
        ss >> mEnv.GrainSize;
      }

      txt = ppParams->XmlReadValue("me:numberOfGrains",false);
      if(txt)
      {
        istringstream ss(txt);
        ss >> mEnv.MaxGrn;
      }

      txt = ppParams->XmlReadValue("me:maxTemperature",false);
      if(txt)
      {
        istringstream ss(txt);
        ss >> mEnv.MaxT;
      }

      if(mEnv.GrainSize!=0.0 && mEnv.MaxGrn!=0)
      {
        stringstream errorMsg;
        errorMsg << "No method is provided to specify me:grainSize and me:numberOfGrains.\n"
                 << "me:numberOfGrains has been ignored";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        mEnv.MaxGrn=0;
      }
    }

    PersistPtr ppControl = ppIOPtr->XmlMoveTo("me:control");
    if(ppControl)
    {
      mEnv.microRateEnabled = ppControl->XmlReadBoolean("me:testMicroRates");
    }

    return true;
  }

  //
  // Begin calculation.
  //
  void System::calculate()
  {
    double beta = .0;
    TimeCount events; string thisEvent; unsigned int timeElapsed =0;

    if(!SetGrainParams())
      return;

    std::string id;
    ModelledMolecule* pmol = NULL;

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
        mEnv.conc = !Pressures.empty() ? Pressures[i]*beta : Concentrations[i];

        // Build collison matrix for system.
        m_pReactionManager->BuildSystemCollisionOperator(beta, mEnv) ;

        m_pReactionManager->diagCollisionOperator() ;
      }
    }
    */

    mEnv.temp = Temperatures[0]; //temporary statements
    mEnv.conc = Concentrations[0];

    beta = 1.0/(boltzmann_RCpK* mEnv.temp) ;

    //---------------
    //About precision
    string precisionMethod;
    switch(precisionTag)
    {
      case varUseDouble:         precisionMethod = "Double";              break;
      case varUseDoubleDouble:   precisionMethod = "Double-double";       break;
      case varUseQuadDouble:     precisionMethod = "Quad-double";         break;
    }
    cout << "Precision: " << precisionMethod << endl;
    //---------------

    // Build collison matrix for system.
    thisEvent = "Build Collison Matrix";
    cout << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;
    m_pReactionManager->BuildSystemCollisionOperator(beta, mEnv) ;

    thisEvent = "Diagonlize Collision Operator";
    cout << endl << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed)  << " -- Time elapsed: " << timeElapsed << " seconds.\n";
    m_pReactionManager->diagCollisionOperator() ;

    thisEvent = "Finish Calculation";
    cout << endl << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed)  << " -- Time elapsed: " << timeElapsed << " seconds.\n";


    /*
    for (size_t i=0; i < m_pReactionManager->size() ; i++) {

    Reaction *reaction = (*m_pReactionManager)[i] ;

    // Work space to hold microcanonical rates.

    vector<double> kmc(mEnv.MaxGrn,0.0) ;

    reaction->get_MicroRateCoeffs(kmc, mEnv) ;

    CollidingMolecule *pmolecule = reaction->m_rct1 ;

    // for (int i(0) ; i < MAXCELL ; i++)
    //     kmc[i] /= omega ;

    pmolecule->diagCollisionOperator() ;

    // Calculate matrix elements

    double kinf = pmolecule->matrixElement(mEnv.MaxGrn-1,mEnv.MaxGrn-1,kmc,mEnv.MaxGrn) ;

    cout << endl ;
    formatFloat(cout, kinf, 6, 15) ;
    }*/

    cout << events << endl;
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
    mEnv.EMax = 0.0; mEnv.EMin = BIG;

    std::string id; ModelledMolecule* pmol = NULL;
    while(m_pMoleculeManager->GetNextMolecule(id, pmol))
    {
      double zpe = pmol->get_zpe() * KcalPerMolToRC ;
      mEnv.EMax = std::max(mEnv.EMax, zpe);
      mEnv.EMin = std::min(mEnv.EMin, zpe);
    }

    if(mEnv.EMin==BIG || mEnv.EMax==0.0)
    {
      stringstream errorMsg;
      errorMsg << "The modelled molecules do not cover an appropriate energy range";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    //mEnv.MaxT gives the option of widening the energy range
    if(mEnv.MaxT <= 0.0){
      mEnv.MaxT = *max_element(Temperatures.begin(), Temperatures.end());
      stringstream errorMsg;
      errorMsg << "Maximum Temperature was invalid. Reset Maximum Temperature, me:maxTemperature, to remove this warning.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }

    //Max energy is ** 20kT  ** above the highest well [was 10kT]
    mEnv.EMax = mEnv.EMax + 20 * boltzmann_RCpK * mEnv.MaxT;
    if(mEnv.GrainSize <= 0.0){
      mEnv.GrainSize = 100.; //default 100cm-1
      stringstream errorMsg;
      errorMsg << "Grain size was invalid. Reset grain size to default: 100";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }

    if(mEnv.MaxGrn <= 0)
      mEnv.MaxGrn = (int)((mEnv.EMax-mEnv.EMin)/mEnv.GrainSize + 0.5);
    else
      mEnv.GrainSize = (mEnv.EMax-mEnv.EMin)/mEnv.MaxGrn;

    mEnv.MaxCell = (int)(mEnv.GrainSize * mEnv.MaxGrn + 0.5);
    return true;
    /*
     //Hardwired
     mEnv.MaxCell=50000;
     mEnv.MaxGrn =500;
     mEnv.GrainSize = 100;
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
    std::string thisEvent, timeString;
    thisEvent = "Write XML attribute";
    timeString = events.setTimeStamp(thisEvent);
    cout << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;
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
