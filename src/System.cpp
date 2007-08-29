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
#include <algorithm>
#include "System.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  System::System(): m_pMoleculeManager(0), m_pReactionManager(0) {
      m_pMoleculeManager = new MoleculeManager(this) ;
      m_pReactionManager = new ReactionManager(this, m_pMoleculeManager) ;
  }

  System::~System() { }

  //
  // Parse an input data file.
  //
  bool System::parse(PersistPtr ppIOPtr)
  {
    GrainSize=0.0;
    MaxGrn=0;
    MaxT=0.0;

    m_ppIOPtr = ppIOPtr;
    PersistPtr ppMolList = ppIOPtr->MoveTo("moleculeList");
    if(!ppMolList)
    {
      stringstream errorMsg;
      errorMsg << "No molecules have been specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    if(!m_pMoleculeManager->addmols(ppMolList))
      return false;;

    PersistPtr ppReacList = ppIOPtr->MoveTo("reactionList");
    if(!ppReacList)
    {
      stringstream errorMsg;
      errorMsg << "No reactions have been specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }if(!m_pReactionManager->addreactions(ppReacList))
        return false;

    PersistPtr ppConditions = ppIOPtr->MoveTo("me:conditions");
    if(!ppConditions)
    {
      stringstream errorMsg;
      errorMsg << "No conditions have been specified";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    const string Bgtxt = ppConditions->ReadValue("me:bathGas");
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
    
    PersistPtr ppParams = ppIOPtr->MoveTo("me:modelParameters");
    
    if(ppParams)
    {
      const char* txt = ppParams->ReadValue("me:grainSize",false);
      if(txt)
      {
        istringstream ss(txt);
        ss >> GrainSize;
      }
      
      txt = ppParams->ReadValue("me:numberOfGrains",false);
      if(txt)
      {
        istringstream ss(txt);
        ss >> MaxGrn;
      }
      
      
      txt = ppParams->ReadValue("me:maxTemperature",false);
      if(txt)
      {
        istringstream ss(txt);
        ss >> MaxT;
      }
      
      if(GrainSize!=0.0 && MaxGrn!=0)
      {
        stringstream errorMsg;
        errorMsg << "No method is provided to specify me:grainSize and me:numberOfGrains.\n"
                 << "me:numberOfGrains has been ignored";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        MaxGrn=0;
      }
    }

    PersistPtr ppControl = ppIOPtr->MoveTo("me:control");
    if(ppControl)
    {
      bTestMicroRates = ppControl->ReadBoolean("me:testMicroRates");
    }

    return true;
  }

  //
  // Begin calculation.
  //
  void System::calculate() 
  { 
    if(!SetGrainParams())
      return;
    
    std::string id;
    ModelledMolecule* pmol=NULL;

    WriteMetadata();

    /* Needs rethink on how to do looping over T and P
    //Outer loop is temperature
    vector<double>::iterator Titer;
    for(Titer=Temperatures.begin();Titer!=Temperatures.end();++Titer)
    {
      double beta = 1.0/(boltzmann*(*Titer)) ;
      
      size_t imax = !Pressures.empty() ? Pressures.size() : Concentrations.size();
      //Inner loop is concentrations, possibly calculated from pressures if these were specified
      //***TODO Get pressure units right. Pressures currently non-functional!
      for(size_t i=0;i<imax;++i)
      {
        double conc = !Pressures.empty() ? Pressures[i]*beta : Concentrations[i];

        // Build collison matrix for system.
        m_pReactionManager->BuildSystemCollisionOperator(beta, conc) ;

        m_pReactionManager->diagCollisionOperator() ;
      }
    }
    */
    temp = Temperatures[0]; //temporary statements
    conc = Concentrations[0];

    double beta = 1.0/(boltzmann*temp) ;

    // Build collison matrix for system.

    m_pReactionManager->BuildSystemCollisionOperator(beta, conc) ;

    m_pReactionManager->diagCollisionOperator() ;


    /*      for (size_t i=0; i < m_pReactionManager->size() ; i++) {

            Reaction *reaction = (*m_pReactionManager)[i] ;

            // Work space to hold microcanonical rates.

            vector<double> kmc(MaxGrn,0.0) ;

            reaction->get_MicroRateCoeffs(kmc) ;

            CollidingMolecule *pmolecule = reaction->m_Reactant ;

            // for (int i(0) ; i < MAXCELL ; i++) 
            //     kmc[i] /= omega ;

            pmolecule->diagCollisionOperator() ;

            // Calculate matrix elements

            double kinf = pmolecule->matrixElement(MaxGrn-1,MaxGrn-1,kmc,MaxGrn) ;

            cout << endl ;
            formatFloat(cout, kinf, 6, 15) ; 

    }*/

  }

  bool System::SetGrainParams()
  {
    /* 
    Either grain size or number of grains can be specified, but not both.

    Uses the value of grain size in the datafile, if specified .
    If grain size is not specified but number of grains is,
    use a grain size to fit the energy range.
    If neither is specified, the grain size is set to 100cm-1
    and the number of grains set so that the energy range is sufficient.

    The required total energy domain extends from the lowest zero point energy
    of the lowest molecule to 10kT above the highest. <me:maxTemperature>
    is used, if specified in the data file.
    But if the input energy of the system is higher, should this be used?
    */

    //Calculate the energy range covering all modelled molecules
    const double BIG =1e100;
    EMax=0.0;
    EMin=BIG;
    std::string id;
    ModelledMolecule* pmol=NULL;
    while(m_pMoleculeManager->GetNextMolecule(id, pmol))
    {
      double zpe = pmol->get_zpe() * KCMLTOPCM ;
      EMax = std::max(EMax, zpe);
      EMin = std::min(EMin, zpe);
    }

    if(EMin==BIG || EMax==0.0)
    {
      stringstream errorMsg;
      errorMsg << "The modelled molecules do not cover an appropriate energy range";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    
    //MaxT gives the option of widening the energy range
    if(MaxT==0.0) 
      MaxT = *max_element(Temperatures.begin(), Temperatures.end());

    //Max energy is ** 20kT  ** above the highest well [was 10kT]
    EMax = EMax + 20 * boltzmann * MaxT;
    if(GrainSize==0.0)
      GrainSize = 100; //default 100cm-1

    if(MaxGrn==0)
      MaxGrn = (int)((EMax-EMin)/GrainSize + 0.5);
    else
      GrainSize = (EMax-EMin)/MaxGrn;

    MaxCell = (int)(GrainSize * MaxGrn + 0.5);
    return true;
    /*
     //Hardwired
     MaxCell=50000;
     MaxGrn =500;
     GrainSize = 100;
     return true;
    */
  }

  bool System::ReadRange(const string& name, vector<double>& vals, PersistPtr ppbase, bool MustBeThere)
  {
    PersistPtr pp=ppbase;
    for(;;)
    {
      const char* txt;
      pp = pp->MoveTo(name);
      if(pp)
        txt = pp->Read(); //element may have a value
      else //no more elements
        break;
      if(!txt)
        txt = pp->ReadValue("initial"); //or use value of "initial" attribute
      if(!txt)
        return false;
      vals.push_back(atof(txt));

      if(txt=pp->ReadValue("increment",false))//optional attribute
      {
        double incr = atof(txt);
        txt = pp->ReadValue("final"); //if have "increment" must have "final"
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
    PersistPtr ppList = m_ppIOPtr->WriteMainElement("metadataList", "");
    PersistPtr ppItem = ppList->WriteElement("metadata");
    ppItem->WriteAttribute("name", "dc:creator");
    ppItem->WriteAttribute("content", "Mesmer v0.1");

    ppItem = ppList->WriteElement("metadata");
    ppItem->WriteAttribute("name", "dc:description");
    ppItem->WriteAttribute("content", 
    "Calculation of the interaction between collisional energy transfer and chemical reaction"
    " for dissociation, isomerization and association processes");

    ppItem = ppList->WriteElement("metadata");
    ppItem->WriteAttribute("name", "dc:date");
    ppItem->WriteAttribute("content", IPersist::TimeString());

    //The user's name should be in an environment variable attached to his account (not a System variable)
    const char* author = getenv("MESMER_AUTHOR");
    if(!author)
      author = "unknown";
    ppItem = ppList->WriteElement("metadata");
    ppItem->WriteAttribute("name", "dc:contributor");
    ppItem->WriteAttribute("content", author);

  }

}//namespace
