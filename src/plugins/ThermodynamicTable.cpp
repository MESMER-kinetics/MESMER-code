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

	ThermodynamicTable(const char* id) : m_id(id),
	  m_nTemp(20),
	  m_TempInterval(50.0),
	  m_Unit("kJ/mol") 
  { Register(); }

	virtual ~ThermodynamicTable() {}
  virtual const char* getID()  { return m_id; }

	// Function to do the work
	virtual bool DoCalculation(System* pSys);

  private:

	// Read any data from XML and store in this instance. 
	bool ReadParameters(PersistPtr ppControl) ;

	string underlineText(const string& text) const ;

	string writeTableHeader(const string& unit) const ;

	void writeTableEntry(Molecule *pmol, double temp, double unitFctr, string & header) const ;

	int m_nTemp ;
	double m_TempInterval ;
	string m_Unit ;

  private:
    const char* m_id;
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

	// Make provision for the special case of T = 298.15.

	bool   tempLessThan298(true) ;
	double temp289(298.15) ;

	// Determine if DOS test information is to appear.

	pSys->m_Flags.testDOSEnabled = ppControl->XmlReadBoolean("me:testDOS");
	if (pSys->m_Flags.testDOSEnabled) {
	  pSys->getEnv().beta = 1.0/(boltzmann_RCpK*double(m_nTemp)*m_TempInterval) ;
	}

	// Begin table.

	ctest << endl ;
	ctest << underlineText(string("Thermodynamic Tables")) ; 

	// Parse molecule data. 

	MoleculeManager* pMoleculeManager = pSys->getMoleculeManager() ;

	PersistPtr ppMolList = pMoleculeManager->get_PersistPtr();
	if (!ppMolList) {
	  cerr << "No molecules have been specified." << endl;
	  return false;
	}

	PersistPtr ppmol = ppMolList ;
	while(ppmol = ppmol->XmlMoveTo("molecule")) {

	  // Get the name of the molcule.
	  const char* reftxt = ppmol->XmlReadValue("id");
	  if (reftxt) {
      //Try molType="" rather than "modelled". Avoid deltaEdown being required.
		pMoleculeManager->addmol(string(reftxt), string(""), pSys->getEnv(), pSys->m_Flags);
	  }
	}

	// Loop over all molecules producing a table for each molecule.

	MoleculeManager::constMolIter molItr = pMoleculeManager->begin() ;
	MoleculeManager::constMolIter molItrEnd = pMoleculeManager->end() ;
	for (; molItr != molItrEnd ; molItr++) {

	  Molecule *pmol = molItr->second;

	  ctest << endl ;
	  ctest << underlineText(pmol->getName()) ;

	  string header = writeTableHeader(m_Unit) ;
	  tempLessThan298 = true ;
	  for (int i(1); i < m_nTemp ; i++) {
		double temp(m_TempInterval*double(i)) ;
		if (tempLessThan298 && temp > temp289) {

		  // Special case of T = 289.15

		  tempLessThan298 = false ;
		  ctest << endl ;
		  writeTableEntry(pmol, temp289, unitFctr, header) ;
		  ctest << endl ;
		}
		writeTableEntry(pmol, temp, unitFctr, header) ;
		if (!(i % 5)) 
		  ctest << endl ;
	  }
	  if (tempLessThan298) {

		// Special case of T = 289.15

		tempLessThan298 = false ;
		ctest << endl ;
		writeTableEntry(pmol, temp289, unitFctr, header) ;
		ctest << endl ;
	  }

	}

	return true ;

  }

  string ThermodynamicTable::underlineText(const string& text) const {

	ostringstream sstrdatum ;
	sstrdatum << " " << text << endl ;
	sstrdatum << " " ;
	for (size_t i(0) ; i < text.size() ; i++ ) 
	  sstrdatum << "-" ;
	sstrdatum << endl ;

	return sstrdatum.str() ;

  }

  string ThermodynamicTable::writeTableHeader(const string& unit) const {

	ostringstream sstrdatum ;

	sstrdatum.setf(ios::right, ios::adjustfield) ;

	sstrdatum << " " ;
	sstrdatum << setw(10) << "Temp" ;
	sstrdatum << setw(15) << "H(T)" ;
	sstrdatum << setw(15) << "S(T)" ;
	sstrdatum << setw(15) << "G(T)" ;
	sstrdatum << endl ;

	ostringstream sstrunit ;
	sstrunit << "(" << unit << ")" ;
	ostringstream sstrunitk ;
	sstrunitk << "(" << unit << "/K)" ;

	ostringstream sstrdatum2 ;
	sstrdatum2 << setw(10) << "(K)" ;
	sstrdatum2 << setw(15) << sstrunit.str()  ;
	sstrdatum2 << setw(15) << sstrunitk.str() ;
	sstrdatum2 << setw(15) << sstrunit.str()  ;

	sstrdatum << underlineText(sstrdatum2.str()) ;

	return sstrdatum.str() ;
  }

  void ThermodynamicTable::writeTableEntry(Molecule *pmol, double temp, double unitFctr, string& header) const {

	double enthalpy(0.0), entropy(0.0), gibbsFreeEnergy(0.0) ;
	pmol->getDOS().thermodynamicsFunctions(temp, unitFctr, enthalpy, entropy, gibbsFreeEnergy) ;

	if (header.length() > 0) {
	  ctest << endl ;
	  ctest << header ;
	  header.clear() ;
	}
	ctest << formatFloat(temp,    6, 11) << formatFloat(enthalpy,        6, 15) 
	  << formatFloat(entropy, 6, 15) << formatFloat(gibbsFreeEnergy, 6, 15) << endl ;
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
		cwarn << "Un-supported unit, the default units of kJ/mol will be used for thermodynamics tables." << endl;
	  } 

	} else {
	  cinfo << "The default units of kJ/mol will be used for thermodynamics tables." << endl;
	}

	int nTemp = ppProp->XmlReadInteger("me:NumberOfTemp", false) ;
	if (nTemp > 0) m_nTemp = nTemp ;

	double TempInterval = ppProp->XmlReadDouble("me:TempInterval", false) ;
	if (TempInterval > 0.0) m_TempInterval = TempInterval ;

	return true ;

  }

}//namespace

