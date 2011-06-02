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

    ThermodynamicTable(const std::string& id) : CalcMethod(id),
      m_nTemp(20),
      m_TempInterval(50.0),
      m_Unit("kJ/mol") {}

    virtual ~ThermodynamicTable() {}

    // Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:

    // Read any data from XML and store in this instance. 
    bool ReadParameters(PersistPtr ppControl) ;

    void underlineText(const string& text) const ;

    void writeTableHeader(const string& text, const string& unit) const ;

    int m_nTemp ;
    double m_TempInterval ;
    string m_Unit ;

  } ;

  ////////////////////////////////////////////////
  //Global instance
  ThermodynamicTable theThermodynamicTable("ThermodynamicTable");
  ///////////////////////////////////////////////

  bool ThermodynamicTable::DoCalculation(System* pSys)
  {

    //Read in fitting parameters, or use values from defaults.xml.
    PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");

    ReadParameters(ppControl) ;

    double unitFctr(1.0/kJPerMol_in_RC) ;
    if (m_Unit == "kcal/mol") 
       unitFctr = 1.0/kCalPerMol_in_RC ;

    MesmerEnv& Env = pSys->getEnv() ;
    Env.MaxCell = 100000 ;

    ctest << endl ;
    underlineText(string("Thermodynamic Tables")) ; 

    const MoleculeManager* pMoleculeManager = pSys->getMoleculeManager() ;

    // Loop over all molecules producing a table for each molecule.

    MoleculeManager::constMolIter molItr = pMoleculeManager->begin() ;
    MoleculeManager::constMolIter molItrEnd = pMoleculeManager->end() ;
    for (; molItr != molItrEnd ; molItr++) {

      Molecule *pmol = molItr->second;
      writeTableHeader(pmol->getName(), m_Unit) ;       
      for (int i(1); i < m_nTemp ; i++) {
        double temp(m_TempInterval*double(i)) ;
        double enthalpy, entropy, gibbsFreeEnergy ;
        pmol->getDOS().thermodynamicsFunctions(temp, unitFctr, enthalpy, entropy, gibbsFreeEnergy) ;
        ctest << formatFloat(temp, 6, 11) << formatFloat(enthalpy, 6, 15) 
          << formatFloat(entropy, 6, 15) << formatFloat(gibbsFreeEnergy, 6, 15) << endl ;
        if (!(i % 5)) 
          ctest << endl ;
      }

    }

    return true ;

  }

  void ThermodynamicTable::underlineText(const string& text) const {

    ctest << " " << text << endl ;
    ctest << " " ;
    for (size_t i(0) ; i < text.size() ; i++ ) 
      ctest << "-" ;
    ctest << endl ;

  }

  void ThermodynamicTable::writeTableHeader(const string& text, const string& unit) const {

    ostringstream sstrdatum ;

    sstrdatum.setf(ios::right, ios::adjustfield) ;

    ctest << "\n" << text << "\n " ;
    sstrdatum << setw(10) << "Temp" ;
    sstrdatum << setw(15) << "H(T)" ;
    sstrdatum << setw(15) << "S(T)" ;
    sstrdatum << setw(15) << "G(T)" ;
    ctest << sstrdatum.str() << endl ;

    ostringstream sstrunit ;
    sstrunit << "(" << unit << ")" ;
    ostringstream sstrunitk ;
    sstrunitk << "(" << unit << "/K)" ;

    sstrdatum.str("") ;
    sstrdatum << setw(10) << "(K)" ;
    sstrdatum << setw(15) << sstrunit.str()  ;
    sstrdatum << setw(15) << sstrunitk.str() ;
    sstrdatum << setw(15) << sstrunit.str()  ;

    underlineText(sstrdatum.str());
  }

  bool ThermodynamicTable::ReadParameters(PersistPtr ppControl) {

    PersistPtr ppProp = ppControl->XmlMoveTo("me:calcMethod") ;

    const char* utxt= ppProp->XmlReadValue("me:Units", false);
    if (utxt) {
      string unit(utxt) ;
      for (size_t i(0) ; i < unit.size(); i++)
        unit[i] = tolower(unit[i]) ;

      if (unit == "kcal/mol"){ 
        m_Unit = unit ;
      } else if (unit == "kj/mol") {
        m_Unit = "kJ/mol" ;
      } else {
        clog << "Un-supported unit, the default units of kJ/mol will be used for thermodynamics tables." ;
      } 

    } else {
      clog << "The default units of kJ/mol will be used for thermodynamics tables." ;
    }

    int nTemp = ppProp->XmlReadInteger("me:NumberOfTemp", false) ;
    if (nTemp > 0) m_nTemp = nTemp ;

    double TempInterval = ppProp->XmlReadDouble("me:TempInterval", false) ;
    if (TempInterval > 0.0) m_TempInterval = TempInterval ;

    return true ;

  }

}//namespace

