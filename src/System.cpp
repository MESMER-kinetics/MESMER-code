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
#include "Persistence.h"
#include "Constants.h"
#include "System.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  //namespace global variable
  System* pSys;

    System::System(): m_pMoleculeManager(0), m_pReactionManager(0) {

        m_pMoleculeManager = new MoleculeManager() ;

        m_pReactionManager = new ReactionManager(m_pMoleculeManager) ;

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

        PersistPtr ppMolList = ppIOPtr->MoveTo("moleculeList");
        if(!ppMolList)
        {
            cerr << "No molecules have been specified" << endl;
            return false;
        }
        if(!m_pMoleculeManager->addmols(ppMolList))
            return false;;

        PersistPtr ppReacList = ppIOPtr->MoveTo("reactionList");
        if(!ppReacList)
        {
            cerr << "No reactions have been specified" << endl;
            return false;
        }if(!m_pReactionManager->addreactions(ppReacList))
            return false;

        PersistPtr ppConditions = ppIOPtr->MoveTo("me:conditions");
        if(!ppConditions)
        {
            cerr << "No conditions have been specified" << endl;
            return false;
        }
        const string Bgtxt = ppConditions->ReadValue("me:bathGas");
        if(!Bgtxt.size() || !(m_pMoleculeManager->find(Bgtxt)) )
        {
            cerr << "No bath gas specified" << endl;
            return false;
		} else {
			m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;
		}

        const char* Ttxt = ppConditions->ReadValue("me:temperature");
        const char* Ctxt = ppConditions->ReadValue("me:conc");
        if(!Ttxt || !Ctxt)
            return false;

        istringstream Tss(Ttxt);
        Tss >> temp;
        istringstream Css(Ctxt);
        Css >> conc;

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
            cerr << "You cannot specify me:grainSize and me:numberOfGrains\n"
                 << "me:numberOfGrains has been ignored" << endl;
            MaxGrn=0;
          }
        }

        return true;
    }

    //
    // Begin calculation.
    //
    void System::calculate() 
    { 
        double beta = 1.0/(boltzmann*temp) ;
        
        if(!SetGrainParams())
          return;
        
        std::string id;
        ModelledMolecule* pmol=NULL;
        while(m_pMoleculeManager->GetNextMolecule(id, pmol))
          pmol->calcDensityOfStates();

        // Build collison matrix for system.

        m_pReactionManager->BuildSystemCollisionOperator(beta, conc) ;

        m_pReactionManager->diagCollisionOperator() ;


        for (size_t i=0; i < m_pReactionManager->size() ; i++) {

            Reaction *reaction = (*m_pReactionManager)[i] ;

            reaction->CalcMicroRateCoeffs() ;

            // Work space to hold microcanonical rates.

            vector<double> kmc(MaxGrn,0.0) ;

            reaction->get_MicroRateCoeffs(kmc, MaxGrn) ;

            CollidingMolecule *pmolecule = reaction->m_Reactant ;

            // for (int i(0) ; i < MAXCELL ; i++) 
            //     kmc[i] /= omega ;

            pmolecule->diagCollisionOperator() ;

            // Calculate matrix elements

            double kinf = pmolecule->matrixElement(MaxGrn-1,MaxGrn-1,kmc,MaxGrn) ;

            cout << endl ;
            formatFloat(cout, kinf, 6, 15) ;

        }

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
        cerr << "The modelled molecules do not cover an appropriate energy range" << endl;
        return false;
      }
      
      //MaxT gives the option of widening the energy range
      if(MaxT==0.0) 
        MaxT = temp;

      //Max energy is ** 10kT  ** above the highest well
      EMax = EMax + 10 * boltzmann * MaxT;
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

}//namespace
