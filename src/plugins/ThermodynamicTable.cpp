//-------------------------------------------------------------------------------------------
//
// ThermodynamicTable.cpp
//
// Author: Struan Robertson
// Date:   06/Mar/2011
//
// This file contains the declaration and implementation of the plug-in class that calculates
// the thermodynamics tables for all the molecules specified in the input file.
//
//-------------------------------------------------------------------------------------------

#include "../System.h"
#include "../calcmethod.h"

namespace mesmer
{
  class ThermodynamicTable : public CalcMethod
  {
  public:
    ThermodynamicTable(const std::string& id) : CalcMethod(id) {}
    virtual ~ThermodynamicTable() {}
    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:
    void ThermodynamicTable::underlineText(const string& text) ;
  } ;

  ////////////////////////////////////////////////
  //Global instance
  ThermodynamicTable theThermodynamicTable("ThermodynamicTable");
  ///////////////////////////////////////////////

  bool ThermodynamicTable::DoCalculation(System* pSys)
  {

    //Read in fitting parameters, or use values from defaults.xml.
    PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");
    // unsigned maxIterations= ppControl->XmlReadInteger("me:fittingIterations");
    double unitFctr(1.0/kJPerMol_in_RC) ;
    double interval(50.0) ;
    int nintervals(60) ;
    
    MesmerEnv& Env = pSys->getEnv() ;
    Env.MaxCell = 100000 ;
    
    underlineText(string("Thermodynamic Tables")) ; 

    const MoleculeManager* pMoleculeManager = pSys->getMoleculeManager() ;

    // Loop over all molecules producing a table for each molecule.

    MoleculeManager::constMolIter molItr = pMoleculeManager->begin() ;
    MoleculeManager::constMolIter molItrEnd = pMoleculeManager->end() ;
    for (; molItr != molItrEnd ; molItr++) {
      Molecule *pmol = molItr->second;

      underlineText(pmol->getName()) ; 
      underlineText(string("     Temp          H(T)          S(T)        G(T)")) ; 
      
      for (int i(1); i < nintervals ; i++) {
        double temp(interval*double(i)) ;
        double enthalpy, entropy, gibbsFreeEnergy ;
        pmol->getDOS().thermodynamicsFunctions(temp, unitFctr, enthalpy, entropy, gibbsFreeEnergy) ;
        ctest << formatFloat(temp, 6, 10) << formatFloat(enthalpy, 6, 15) 
          << formatFloat(entropy, 6, 15) << formatFloat(gibbsFreeEnergy, 6, 15) << endl ;
        if (!(i % 5)) 
          ctest << endl ;
      }

    }

    return true ;

  }

  void ThermodynamicTable::underlineText(const string& text) {

    ctest << endl ;
    ctest << " " << text << endl ;
    ctest << " " ;
    for (size_t i(0) ; i < text.size() ; i++ ) 
      ctest << "-" ;
    ctest << endl ;

  }

}//namespace

