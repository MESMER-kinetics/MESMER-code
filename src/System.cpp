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
#include <fstream>

#include "AssociationReaction.h"
#include "IrreversibleUnimolecularReaction.h"
#include "IsomerizationReaction.h"
#include "IrreversibleExchangeReaction.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{

  System::System(const std::string& libraryfilename): 
m_pMoleculeManager(0), 
m_pReactionManager(0), 
m_pTitle(NULL),
m_pDescription(NULL),
m_isomers(),
m_sources(),
m_sinkRxns(),
m_SinkSequence(),
m_meanOmega(0.0),
m_reactionOperator(0),
m_eigenvectors(0),
m_eigenvalues(),
m_SpeciesSequence(),
m_eqVector()
{
  m_pMoleculeManager = new MoleculeManager(libraryfilename) ;
  m_pReactionManager = new ReactionManager(m_pMoleculeManager) ;
}

System::~System() {
  if (m_reactionOperator) delete m_reactionOperator;
  delete m_pReactionManager;
  delete m_pMoleculeManager;
}

//
// Parse an input data file.
//
bool System::parse(PersistPtr ppIOPtr)
{
  initializeConversionMaps();
  m_ppIOPtr = ppIOPtr;

  m_pTitle       = ppIOPtr->XmlReadValue("title", false);
  m_pDescription = ppIOPtr->XmlReadValue("description", false);

  //-------------
  //Molecule List (parse this part inside Reaction)
  PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
  m_pMoleculeManager->set_PersistPtr(ppMolList);

  //-------------
  //Model Parameters
  PersistPtr ppParams = ppIOPtr->XmlMoveTo("me:modelParameters");
  if(ppParams)
  {
	m_Env.GrainSize          = ppParams->XmlReadInteger("me:grainSize");

	m_Env.MaximumTemperature = ppParams->XmlReadDouble("me:maxTemperature",optional);
	m_Env.EAboveHill         = ppParams->XmlReadDouble("me:energyAboveTheTopHill");
	m_Env.useBasisSetMethod  = ppParams->XmlReadBoolean("me:runBasisSetMethodroutines");
	if (m_Env.useBasisSetMethod) {
	  PersistPtr ppBasisSet = ppParams->XmlMoveTo("me:runBasisSetMethodroutines");
	  if(ppBasisSet) {
		m_Env.nBasisSet = ppBasisSet->XmlReadInteger("me:numberBasisFunctions");
	  } else {
		cerr << "Basis set method requested but number of basis functions unspecified.";
		return false;
	  }
	}
  }
  cinfo.flush();

  //-------------
  //Reaction List
  PersistPtr ppReacList = ppIOPtr->XmlMoveTo("reactionList");
  if(!ppReacList)
  {
	cerr << "No reactions have been specified";
	return false;
  }
  if(!m_pReactionManager->addreactions(ppReacList, m_Env, m_Flags)) return false;

  //Check that the energy baseline is the same for all the modelled molecules
  string energyConvention = m_pMoleculeManager->checkEnergyConventions();
  if(energyConvention.empty())
  {
	cerr << "Not all the molecule energies use the same baseline and need to be.\n";
	return false;
  }
  else
	cinfo << "All molecules are on the same energy basis: " << energyConvention << endl;
  cinfo.flush();
  //-------------

  //Reaction Conditions
  PersistPtr ppConditions = ppIOPtr->XmlMoveTo("me:conditions");
  if(!ppConditions)
  {
	cerr << "No conditions specified";
	return false;
  }
  string Bgtxt = ppConditions->XmlReadValue("me:bathGas");
  if(!m_pMoleculeManager->addmol(Bgtxt, "bathGas", ppMolList, m_Env, m_Flags))
	return false;
  m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;

  //--------------
  //  The concentration/pressure units are of following formats:
  //  units:
  //   PPCC: particles per cubic centimeter
  //   Torr: Torr
  //
  //  Allowed input formats are shown below (example units in particles per cubic centimeter).
  //
  //  <me:PTs>
  //    <me:PTset me:units="Torr">
  //      <me:Prange initial="1e8" increment="2e7" final="2e8" />
  //      <me:Trange initial="100" increment="20" final="200" />
  //    </me:PTset>
  //  </me:PTs>
  //
  //  The above example creates a matrix of concentration/temperature points of the size:
  //      (number of P points) x (number of T points)
  //
  //  Another example of specifying small numbers of PT points (example units in Torr):
  //
  //  <me:PTs>
  //    <me:PTpair me:units="Torr" me:P="100" me:T="200" />
  //    <me:PTpair me:units="PPCC" me:P="1e18" me:T="298" />
  //  </me:PTs>
  //
  //  The looping of the PT points are easy, they are first parsed and all the points are stored in pairs in
  //  vector PandTs, and Mesmer simply loop through all its members.
  //
  //  Or, if someone wants a higher precision on some condition they are interested, they can specify additional
  //  precision flag with small numbers of PT points.
  //
  //  <me:PTs>
  //    <me:PTpair me:units="Torr" me:P="100" me:T="200" me:precision="double-double" />
  //    <me:PTpair me:units="PPCC" me:P="1e18" me:T="298" me:precision="quad-double" />
  //  </me:PTs>
  //
  //  The description above will do first a double-double precision calculation and a quad-double calculation.
  //--------------

  PersistPtr ppPTs = ppConditions->XmlMoveTo("me:PTs");
  if(ppPTs)
	readPTs(ppPTs);
  if (!PandTs.size())
	cerr << "No pressure and temperature specified.";

  // read initial population (needs to be normalized later if their sum not equals to 1.0)
  PersistPtr ppInitalPopulation = ppConditions->XmlMoveTo("me:InitalPopulation");
  if (ppInitalPopulation)
	m_pReactionManager->setInitialPopulation(ppInitalPopulation);

  PersistPtr ppControl = ppIOPtr->XmlMoveTo("me:control");
  if(ppControl)
  {
	m_Flags.testDOSEnabled              = ppControl->XmlReadBoolean("me:testDOS");
	m_Flags.testRateConstantEnabled     = ppControl->XmlReadBoolean("me:testRateConstants");
	m_Flags.microRateEnabled            = ppControl->XmlReadBoolean("me:testMicroRates");
	m_Flags.grainBoltzmannEnabled       = ppControl->XmlReadBoolean("me:printGrainBoltzmann");
	m_Flags.grainDOSEnabled             = ppControl->XmlReadBoolean("me:printGrainDOS");
	m_Flags.grainTSsosEnabled           = ppControl->XmlReadBoolean("me:printTSsos");
	m_Flags.cellDOSEnabled              = ppControl->XmlReadBoolean("me:printCellDOS");
	m_Flags.reactionOCSEnabled          = ppControl->XmlReadBoolean("me:printReactionOperatorColumnSums");
	m_Flags.kfEGrainsEnabled            = ppControl->XmlReadBoolean("me:printGrainkfE");
	m_Flags.kbEGrainsEnabled            = ppControl->XmlReadBoolean("me:printGrainkbE");
	// Both Tunnelling and Tunneling will work
	m_Flags.TunnellingCoeffEnabled      = ppControl->XmlReadBoolean("me:printTunnellingCoefficients");
	if (!m_Flags.TunnellingCoeffEnabled)
	  m_Flags.TunnellingCoeffEnabled    = ppControl->XmlReadBoolean("me:printTunnelingCoefficients");
	m_Flags.CrossingCoeffEnabled        = ppControl->XmlReadBoolean("me:printCrossingCoefficients");
	m_Flags.cellFluxEnabled             = ppControl->XmlReadBoolean("me:printCellTransitionStateFlux");
	m_Flags.grainFluxEnabled            = ppControl->XmlReadBoolean("me:printGrainTransitionStateFlux");
	m_Flags.rateCoefficientsOnly        = ppControl->XmlReadBoolean("me:calculateRateCoefficientsOnly");
	m_Flags.useTheSameCellNumber        = ppControl->XmlReadBoolean("me:useTheSameCellNumberForAllConditions");
	m_Flags.grainedProfileEnabled       = ppControl->XmlReadBoolean("me:printGrainedSpeciesProfile");
	m_Flags.speciesProfileEnabled       = ppControl->XmlReadBoolean("me:printSpeciesProfile");
	m_Flags.viewEvents                  = ppControl->XmlReadBoolean("me:printEventsTimeStamps");
	m_Flags.allowSmallerDEDown          = ppControl->XmlReadBoolean("me:allowSmallerDeltaEDown");
	m_Flags.print_TabbedMatrices        = ppControl->XmlReadBoolean("me:printTabbedMatrices");
	m_Flags.useDOSweightedDT             = ppControl->XmlReadBoolean("me:useDOSweighedDownWardTransition");
	if (!m_Flags.useTheSameCellNumber && m_Env.MaximumTemperature != 0.0){
	  m_Flags.useTheSameCellNumber = true;
	}

	// System configuration information
	if (ppControl->XmlReadBoolean("me:runPlatformDependentPrecisionCheck")) configuration();

	m_CalcMethod = CalcMethod::GetCalcMethod(ppControl);

	//if (m_Flags.grainedProfileEnabled && (m_Flags.speciesProfileEnabled)){
	//  cinfo << "Turn off grained species profile to prevent disk flooding." << endl;
	//  m_Flags.grainedProfileEnabled = false;
	//}

	const char* txtPCOP = ppControl->XmlReadValue("me:printCollisionOperatorLevel",false);
	if(txtPCOP) {
	  istringstream ss(txtPCOP);
	  ss >> m_Flags.showCollisionOperator;
	}

	const char* txtEV = ppControl->XmlReadValue("me:eigenvalues",false);
	if(txtEV) {
	  istringstream ss(txtEV);
	  ss >> m_Flags.printEigenValuesNum;
	}

	const char* txtROS = ppControl->XmlReadValue("me:printReactionOperatorSize",optional);
	if(txtROS) {
	  istringstream sROS(txtROS);
	  sROS >> m_Flags.printReactionOperatorNum;
	}

	const char* txtSTI = ppControl->XmlReadValue("me:shortestTimeOfInterest", optional);
	if (txtSTI){
	  istringstream ss(txtSTI);
	  ss >> m_Flags.shortestTimeOfInterest;
	}

	const char* txtMET = ppControl->XmlReadValue("me:MaximumEvolutionTime", optional);
	if (txtMET){
	  istringstream ss(txtMET);
	  ss >> m_Flags.maxEvolutionTime;
	}
  }

  return true;
}

//
// Main calculation method.
//
void System::executeCalculation()
{
  assert(m_CalcMethod);
  m_CalcMethod->DoCalculation(this);
}

// pop the P and T points into PandTs
// This is a function for reading concentration/pressure and temperature conditions.
void System::readPTs(PersistPtr anchor)
{
  PersistPtr pp=anchor;
  const char* txt;

  //default unit, pressure and temperature
  const string default_unit = "PPCC";
  const double default_P = 1e17;
  const double default_T = 298.;

  //
  // check for grid values of temperatures and concentrations
  //
  PersistPtr ppPTset = pp->XmlMoveTo("me:PTset");
  while (ppPTset){
	txt = ppPTset->XmlReadValue("me:units",optional);
	string this_units = txt;
	if (!txt){
	  cerr << "No units provided. Default units " << default_unit << " are used.";
	  this_units = default_unit;
	}

	//
	// If user does not input any value for temperature and concentration,
	// give a Default set of concentration and pressurefor simulation.
	// 
	std::vector<double> Pvals, Tvals;
	if(!ReadRange("me:Prange", Pvals, ppPTset)) Pvals.push_back(default_P);
	if(!ReadRange("me:Trange", Tvals, ppPTset)) Tvals.push_back(default_T);

	for (unsigned int i = 0; i < Pvals.size(); ++i){
	  for (unsigned int j = 0; j < Tvals.size(); ++j){
		CandTpair thisPair(getConvertedP(this_units, Pvals[i], Tvals[j]), Tvals[j]);
		PandTs.push_back(thisPair);
	  }
	}
	ppPTset = ppPTset->XmlMoveTo("me:PTset");
  }

  //
  // check for indivually specified concentration/temperature points
  //
  PersistPtr ppPTpair = pp->XmlMoveTo("me:PTpair");
  while (ppPTpair){
	string this_units;
	txt = ppPTpair->XmlReadValue("me:units", optional);
	if (txt)
	  this_units = txt;
	double this_P = default_P;
	double this_T = default_T;
	int this_precision = 0;

	txt = ppPTpair->XmlReadValue("me:P", optional);
	if (txt){ stringstream s1(txt); s1 >> this_P; }
	txt = ppPTpair->XmlReadValue("me:T", optional);
	if (txt){ stringstream s1(txt); s1 >> this_T; }
	txt = ppPTpair->XmlReadValue("me:precision");
	// Can specify abbreviation
	if (txt){
	  if (!strcmp(txt,"1d")) this_precision = 0;
	  if (!strcmp(txt,"double")) this_precision = 0;
	  if (!strcmp(txt,"2d")) this_precision = 1;
	  if (!strcmp(txt,"dd")) this_precision = 1;
	  if (!strcmp(txt,"double-double")) this_precision = 1;
	  if (!strcmp(txt,"4d")) this_precision = 2;
	  if (!strcmp(txt,"qd")) this_precision = 2;
	  if (!strcmp(txt,"quad-double")) this_precision = 2;
	}
	CandTpair thisPair(getConvertedP(this_units, this_P, this_T), this_T, this_precision);
	cinfo << this_P << this_units << ", " << this_T << "K at " << txt << " precision" <<endl; 

	// Set experimental conditions for chiSquare calculation
	txt = ppPTpair->XmlReadValue("me:experimentalRate", false);
	PersistPtr ppExpRate = ppPTpair->XmlMoveTo("me:experimentalRate");
	while (ppExpRate){
	  double rateValue(0.0), errorValue(0.0); 
	  string ref1, ref2, refReaction;
	  stringstream s1(txt); s1 >> rateValue;
	  txt = ppExpRate->XmlReadValue("ref1");
	  stringstream s2(txt); s2 >> ref1;
	  txt = ppExpRate->XmlReadValue("ref2");
	  stringstream s3(txt); s3 >> ref2;
	  txt = ppExpRate->XmlReadValue("refReaction", false);
	  if (txt) {
		stringstream s3(txt); s3 >> refReaction ;
	  }
	  txt = ppExpRate->XmlReadValue("error");
	  stringstream s4(txt); s4 >> errorValue;
	  thisPair.set_experimentalRate(ref1, ref2, refReaction, rateValue, errorValue);
	  ppExpRate = ppExpRate->XmlMoveTo("me:experimentalRate");
	}

	PandTs.push_back(thisPair);
	ppPTpair = ppPTpair->XmlMoveTo("me:PTpair");
  }
}

//
// Begin calculation.
// over all PT values, constant parameters
bool System::calculate(double& chiSquare, bool writeReport)
{
  // Controls the print-out of grain/cell DOS in each cycle (This is only for source term)
  if (m_Flags.cellDOSEnabled) m_Flags.cyclePrintCellDOS = true;
  if (m_Flags.grainDOSEnabled) m_Flags.cyclePrintGrainDOS = true;

  chiSquare = 0.0; // reset the value to zero

  TimeCount events; unsigned int timeElapsed =0;

  // Find the highest temperature
  for (unsigned int i = 0; i < PandTs.size(); ++i){
	m_Env.MaximumTemperature = max(m_Env.MaximumTemperature, PandTs[i].get_temperature());
  }

  //---------------------------------------------
  // looping over temperatures and concentrations
  unsigned int calPoint(0);
  //XML output
  //Considered putting this output under each PT pair in me:conditions.
  //But doesn't work with a range of Ps or Ts. So has to have its own section.

  //There will usually be an <analysis> section for every calculate()
  //When fitting set m_Flags.overwriteXmlAnalysis true
  string comment = m_Flags.overwriteXmlAnalysis ?
	"Only selected calculations shown here" : "All calculations shown";
  PersistPtr ppAnalysis = m_ppIOPtr->XmlWriteMainElement("me:analysis", comment, m_Flags.overwriteXmlAnalysis);
  if(Rdouble::withRange().size()!=0)
  {
	PersistPtr ppParams = ppAnalysis->XmlWriteElement("me:parameters");
	for(size_t i=0;i!=Rdouble::withRange().size();++i)
	{
	  stringstream ss;
	  ss << *Rdouble::withRange()[i];
	  ppParams->XmlWriteAttribute(Rdouble::withRange()[i]->get_varname(), ss.str());
	}
  }

  stringstream rateCoeffTable ;

  rateCoeffTable << endl ;
  rateCoeffTable << "    Temperature  Concentration    Exp. Coeff.    Cal. Coeff." << endl ;
  rateCoeffTable << endl ;

  for (calPoint = 0; calPoint < PandTs.size(); ++calPoint) {

	m_Env.beta = 1.0 / (boltzmann_RCpK * PandTs[calPoint].get_temperature()) ; //temporary statements
	m_Env.conc = PandTs[calPoint].get_concentration();
	// unit of conc: particles per cubic centimeter
	//clog << "PT Grid " << calPoint << endl;
	cinfo << "PT Grid " << calPoint << endl;
	int precision = PandTs[calPoint].get_precision();
	ctest << "PT Grid " << calPoint << " Condition: conc = " << m_Env.conc << ", temp = " << PandTs[calPoint].get_temperature();

	switch (precision) {
	  case 1: ctest << ", diagonalization precision: double-double\n{\n"; break;
	  case 2: ctest << ", diagonalization precision: quad-double\n{\n"; break;
	  default: ctest << ", diagonalization precision: double\n{\n";
	}

	// Build collison matrix for system.
	if (!BuildReactionOperator(m_Env, m_Flags))
	  throw (std::runtime_error("Failed building system collison operator.")); 

	{string thisEvent = "Build Collison Operator" ;
	events.setTimeStamp(thisEvent, timeElapsed) ;
	cinfo << thisEvent << " -- Time elapsed: " << timeElapsed << " seconds." << endl ;}

	if (!m_Flags.rateCoefficientsOnly){

	  if (!calculateEquilibriumFractions(m_Env.beta))
		throw (std::runtime_error("Failed calculating equilibrium fractions.")); 

	  // Diagonalise the collision operator.
	  diagReactionOperator(m_Flags, m_Env, precision, ppAnalysis) ;

	  {string thisEvent = "Diagonalize the Reaction Operator" ;
	  events.setTimeStamp(thisEvent, timeElapsed) ;
	  cinfo << thisEvent << " -- Time elapsed: " << timeElapsed << " seconds." << endl ;}

	  // Locate all sink terms.
	  locateSinks() ;

	  if (!m_Env.useBasisSetMethod) {

		PersistPtr ppPopList;
		if(m_Flags.speciesProfileEnabled)
		{
		  ppPopList  = ppAnalysis->XmlWriteElement("me:populationList");
		  ppPopList->XmlWriteAttribute("T", toString(PandTs[calPoint].get_temperature()));
		  ppPopList->XmlWriteAttribute("conc", toString(m_Env.conc));
		}
		// Time steps loop
		timeEvolution(m_Flags, ppPopList); 
		PersistPtr ppList = ppAnalysis->XmlWriteElement("me:rateList");
		ppList->XmlWriteAttribute("T", toString(PandTs[calPoint].get_temperature()));
		ppList->XmlWriteAttribute("conc", toString(m_Env.conc));
		ppList->XmlWriteAttribute("me:units", "s-1");
		qdMatrix mesmerRates(1);
		BartisWidomPhenomenologicalRates(mesmerRates, m_Flags, ppList);

		vector<conditionSet> expRates;
		PandTs[calPoint].get_experimentalRates(expRates);

		rateCoeffTable << formatFloat(PandTs[calPoint].get_temperature(), 6, 15) ;
		rateCoeffTable << formatFloat(PandTs[calPoint].get_concentration(), 6, 15) ;

		chiSquare += calcChiSquare(mesmerRates, expRates, rateCoeffTable);

		ctest << "}\n";

	  } else {

		qdMatrix mesmerRates(1);
		BartisWidomBasisSetRates(mesmerRates, m_Flags);

	  }
	}

  } // End of temperature and concentration loop. 

  rateCoeffTable << endl ;

  if (writeReport) cinfo << rateCoeffTable.str() ;

  {string thisEvent = "Finish Calculation";
  events.setTimeStamp(thisEvent, timeElapsed);
  cinfo << endl << thisEvent << " -- Time elapsed: " << timeElapsed << " seconds.\n";
  cwarn << calPoint << " temperature/concentration-pressure points calculated." << endl;}

  if (m_Flags.viewEvents) cinfo << events;

  return true;
}

double System::calcChiSquare(const qdMatrix& mesmerRates, vector<conditionSet>& expRates, stringstream &rateCoeffTable){

  double chiSquare(0.0) ;

  for (size_t i(0); i < expRates.size(); ++i){
	string ref1, ref2, refReaction; 
	double expRate(0.0), expErr(0.0); 
	int seqMatrixLoc1(-1), seqMatrixLoc2(-1);

	expRates[i].get_conditionSet(ref1, ref2, refReaction, expRate, expErr);

	// Get the position of ref1 and ref2 inside m_SpeciesSequence vector
	seqMatrixLoc1 = getSpeciesSequenceIndex(ref1);
	seqMatrixLoc2 = getSpeciesSequenceIndex(ref2);

	if(seqMatrixLoc1<0 || seqMatrixLoc2<0)
	  throw(std::runtime_error("Failed to locate species in rate coefficient matrix.")) ;

	// 
	// In the following it is assumed that experimental rate coefficients will always 
	// be quoted as a absolute values. Since the diagonal values of the BW matrix are
	// negative, their absolute value is required for comparision with experimental
	// values hence the fabs invocation.
	//
	double rateCoeff = fabs(to_double(mesmerRates[seqMatrixLoc2][seqMatrixLoc1])) ;

	// Is a bimolecular rate coefficient required?

	Reaction *reaction = m_pReactionManager->find(refReaction) ;
	if (reaction && reaction->getReactionType() == ASSOCIATION ) {
	  double concExcessReactant = reaction->get_concExcessReactant() ;

	  // Test concentration and reaction sense.

	  if (concExcessReactant > 0.0 && (reaction->get_reactant()->getName() == ref1)) {
		rateCoeff /= concExcessReactant ;
	  }
	} else {
	  // No reference reaction. Assume reaction is unimolecular.
	}

	double diff = rateCoeff  - expRate ;
	chiSquare += (diff * diff) / (expErr * expErr);

	rateCoeffTable << formatFloat(expRate, 6, 15) << formatFloat(rateCoeff, 6, 15) << endl ;

  }
  return chiSquare;
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
	cerr << "Must specify at least one value of " << name;
	return false;
  }
  return true;
}

void System::WriteMetadata(const string& infilename)
{
  PersistPtr ppList = m_ppIOPtr->XmlWriteMainElement("metadataList", "");
  ppList->XmlWriteAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/");
  if(m_pTitle)
    ppList->XmlWriteValueElement("dc:title", m_pTitle);
  if(m_pDescription)
    ppList->XmlWriteValueElement("dc:description", m_pDescription);
  ppList->XmlWriteValueElement("dc:source", infilename);

  ppList->XmlWriteValueElement("dc:creator","Mesmer v"+string(MESMER_VERSION));
  TimeCount events;
  string timeString = events.setTimeStamp("");
  cinfo << "Write metadata " << timeString << endl;
  ppList->XmlWriteValueElement("dc:date", timeString);

  //The user's name should be in an environment variable attached to his account (not a System variable)
  const char* author = getenv("MESMER_AUTHOR");
  if(!author)
    author = "unknown";
  ppList->XmlWriteValueElement("dc:contributor", author);
}

  void System::configuration(void){
  clog << "\nPrinting system precision configuration:" << endl;
  clog << "Size of float = " << sizeof(float) << endl;
  clog << "Size of double = " << sizeof(double) << endl;
  clog << "Size of long double = " << sizeof(long double) << endl;
  clog << "Size of double-double = " << sizeof(dd_real) << endl;
  clog << "Size of quad-double = " << sizeof(qd_real) << endl;

  clog << "\nEpsilon is the difference between 1 and the smallest value greater than 1 that is representable for the data type." << endl;
  clog << "float epsilon == " << numeric_limits<float>::epsilon() << endl;
  clog << "double epsilon == " << numeric_limits<double>::epsilon() << endl;
  clog << "long double epsilon == " << numeric_limits<long double>::epsilon() << endl;
  clog << "dd_real epsilon == " << numeric_limits<dd_real>::epsilon() << endl;
  clog << "qd_real epsilon == " << numeric_limits<qd_real>::epsilon() << endl;

  clog << "\nfloat max == " << numeric_limits<float>::max() << endl;
  clog << "double max == " << numeric_limits<double>::max() << endl;
  clog << "long double max == " << numeric_limits<long double>::max() << endl;
  clog << "dd_real max == " << numeric_limits<dd_real>::max() << endl;
  clog << "qd_real max == " << numeric_limits<qd_real>::max() << endl;

  clog << "\nfloat min == " << numeric_limits<float>::min() << endl;
  clog << "double min == " << numeric_limits<double>::min() << endl;
  clog << "long double min == " << numeric_limits<long double>::min() << endl;
  clog << "dd_real min == " << numeric_limits<dd_real>::min() << endl;
  clog << "qd_real min == " << numeric_limits<qd_real>::min() << endl;
}


bool System::BuildReactionOperator(MesmerEnv &mEnv, MesmerFlags& mFlags)
{
  const double SUPREMUM =  9e23 ;
  const double INFIMUM  = -SUPREMUM ;
  //
  // Find all the unique wells and lowest zero point energy.
  //
  m_isomers.clear();

  double minEnergy(SUPREMUM) ; // The minimum & maximum ZPE amongst all wells, set artificially large and small
  double maxEnergy(INFIMUM) ;  // to guarantee that each is overwritten in setting minEnergy and maxEnergy.
  Molecule *pBathGasMolecule = m_pMoleculeManager->get_BathGasMolecule();

  // populate molMapType with unimolecular species and determine minimum/maximum energy on the PES
  for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
	double TS_ZPE(INFIMUM);

	Reaction *pReaction = (*m_pReactionManager)[i] ;

	// Reset the DOS calculation flags before building the reaction operator.
	pReaction->resetCalcFlag();

	// Transition State
	// third check for the transition state in this reaction
	Molecule *pTransitionState = pReaction->get_TransitionState();
	if (pTransitionState){
	  TS_ZPE = pTransitionState->getDOS().get_zpe();
	  maxEnergy = max(maxEnergy, TS_ZPE) ;
	}

	// unimolecular species
	vector<Molecule *> unimolecules ;
	pReaction->get_unimolecularspecies(unimolecules) ;
	// populate molMapType with unimolecular species
	for (size_t j(0) ; j < unimolecules.size() ; ++j) {
	  // wells
	  Molecule *pCollidingMolecule = unimolecules[j] ;
	  const double collidingMolZPE(pCollidingMolecule->getDOS().get_zpe());
	  if(pCollidingMolecule && m_isomers.find(pCollidingMolecule) == m_isomers.end()){ // New isomer
		m_isomers[pCollidingMolecule] = 0 ; //initialize to a trivial location

		minEnergy = min(minEnergy, collidingMolZPE) ;
		maxEnergy = max(maxEnergy, collidingMolZPE) ;
	  }

	  //calculate the lowest barrier associated with this well(species)
	  if (TS_ZPE != INFIMUM){
		const double barrierHeight = TS_ZPE - collidingMolZPE;
		if (barrierHeight < pCollidingMolecule->getColl().getLowestBarrier()){
		  pCollidingMolecule->getColl().setLowestBarrier(barrierHeight);
		}
	  }
	}

	//
	// For Association reactions determine zero point energy location of the
	// associating pair.
	//
	AssociationReaction *pAReaction = dynamic_cast<AssociationReaction*>(pReaction) ;
	if (pAReaction) {
	  double pseudoIsomerZPE = pAReaction->get_pseudoIsomer()->getDOS().get_zpe();
	  double excessReactantZPE = pAReaction->get_excessReactant()->getDOS().get_zpe();
	  double sourceTermZPE = pseudoIsomerZPE + excessReactantZPE;
	  minEnergy = min(minEnergy, sourceTermZPE) ;
	  maxEnergy = max(maxEnergy, sourceTermZPE) ;

	  // Calculate the lowest barrier associated with this well(species)
	  // For association reaction, it is assumed that the barrier height is close to the source term energy
	  // and in a sense, it is preferable to set this variable to the source term energy even there is an explicit
	  // transition state.
	  double adductZPE = unimolecules[0]->getDOS().get_zpe();
	  double barrierHeight = sourceTermZPE - adductZPE;
	  if (barrierHeight < unimolecules[0]->getColl().getLowestBarrier()){
		unimolecules[0]->getColl().setLowestBarrier(barrierHeight);
	  }
	}

	//
	// For irreversible exchange reactions determine zero point energy location of the
	// associating pair.
	//
	IrreversibleExchangeReaction *pIEReaction = dynamic_cast<IrreversibleExchangeReaction*>(pReaction) ;
	if (pIEReaction) {
	  double pseudoIsomerZPE = pIEReaction->get_pseudoIsomer()->getDOS().get_zpe();
	  double excessReactantZPE = pIEReaction->get_excessReactant()->getDOS().get_zpe();
	  double sourceTermZPE = pseudoIsomerZPE + excessReactantZPE;
	  minEnergy = min(minEnergy, sourceTermZPE) ;
	  maxEnergy = max(maxEnergy, sourceTermZPE) ;

	  // There is no well for this reaction
	}

	//
	// For dissociation reactions determine zero point energy location of the barrier
	//
	IrreversibleUnimolecularReaction *pDissnRtn = dynamic_cast<IrreversibleUnimolecularReaction*>(pReaction) ;
	if (pDissnRtn) {
	  const double rctZPE = pDissnRtn->get_reactant()->getDOS().get_zpe();
	  double barrierZPE = rctZPE + pDissnRtn->get_ThresholdEnergy();
	  minEnergy = min(minEnergy, barrierZPE) ;
	  maxEnergy = max(maxEnergy, barrierZPE) ;

	  // Calculate the lowest barrier associated with this well(species).
	  if (barrierZPE < unimolecules[0]->getColl().getLowestBarrier()){
		unimolecules[0]->getColl().setLowestBarrier(barrierZPE);
	  }
	}

  }

  // set grain parameters for the current Temperature/pressure condition
  if(!SetGrainParams(mEnv, mFlags, minEnergy, maxEnergy))
	return false;

  // Calculate flux and k(E)s
  for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
	if(!(*m_pReactionManager)[i]->calcGrnAvrgMicroRateCoeffs())
	  return false;
  }

  if (!mFlags.rateCoefficientsOnly){
	//
	// Shift all wells to the same origin, calculate the size of the reaction operator,
	// calculate the mean collision frequency and initialize all collision operators.
	//
	int msize(0) ; // size of the collision matrix
	m_meanOmega = 0.0;

	Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
	for (; isomeritr != m_isomers.end() ; ++isomeritr) {

	  Molecule *isomer = isomeritr->first ;
	  isomeritr->second = msize ; //set location

	  int grnZpe = isomer->getColl().get_grnZPE() ; //set grain ZPE (with respect to the minimum of all wells)

	  int colloptrsize = mEnv.MaxGrn - grnZpe ;
	  isomer->getColl().set_colloptrsize(colloptrsize) ;
	  msize += colloptrsize ;

	  if(!isomer->getColl().initCollisionOperator(mEnv.beta, pBathGasMolecule)){
		cerr << "Failed initializing collision operator for " << isomer->getName();
		return false;
	  }

	  // update the size of the collision operator if it is different.
	  int nGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
	  if (nGroupedGrains != 0){
		msize -= (nGroupedGrains - 1);
		if (isomer->getColl().isCemetery()) msize -= 1;
	  }

	  m_meanOmega += isomer->getColl().get_collisionFrequency() ;
	}
	m_meanOmega /= double(m_isomers.size());

	//
	// Find all source terms. Note: a source term contains the deficient reactant.
	// It is possible for there to be more than one source term.
	//
	m_sources.clear(); // Maps the location of source in the system matrix.
	for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
	  (*m_pReactionManager)[i]->updateSourceMap(m_sources) ;
	}

	// Build reaction operator.
	//
	// One of two methods for building the reaction operator are available:
	// the conventional energy grained master equation method which is based
	// on energy grains and a contracted basis set method in which a basis
	// set is generated from the individual collision operators and a
	// representation of the reaction operator build upon this basis.

	if (!mEnv.useBasisSetMethod) {

	  // Full energy grained reaction operator.

	  constructGrainMatrix(msize);

	} else {

	  // Contracted basis set reaction operator.

	  constructBasisMatrix();

	}
  }

  return true;
}

bool System::SetGrainParams(MesmerEnv &mEnv, const MesmerFlags& mFlags, const double minEne, const double maxEne)
{
  //  Grain size and number of grain:
  //
  //  - Either grain size or number of grains can be specified, but not both.
  //
  //  - Uses the value of grain size in the datafile, if specified.
  //
  //  - If grain size is not specified but number of grains is, use a grain size to fit the energy range.
  //  If neither is specified, the grain size is set to 100cm-1 and the number of grains set so that
  //  the energy range is sufficient.
  //
  //  Energy Range:
  //
  //  - The required total energy domain extends from the lowest zero point energy of the lowest molecule
  //  to 10 k_B T above the highest.

  mEnv.EMin = minEne;
  mEnv.EMax = maxEne;

  /*For testing purposes, set the maxGrn based on the highest temperature we use in all calculations.*/
  const double MaximumTemperature = mEnv.MaximumTemperature;

  /*EAboveHill: Max energy above the highest hill. The temperature refers to the current condition.*/
  if (mFlags.useTheSameCellNumber){
	mEnv.EMax += mEnv.EAboveHill * MaximumTemperature * boltzmann_RCpK;
  }
  else{
	mEnv.EMax += mEnv.EAboveHill / mEnv.beta;
  }

  if(mEnv.GrainSize <= 0.0){
	mEnv.GrainSize = 100; //default 100cm-1
	cerr << "Grain size was invalid. Reset grain size to default: 100";
  }

  mEnv.MaxGrn = (int)((mEnv.EMax-mEnv.EMin)/mEnv.GrainSize + 0.5);
  mEnv.MaxCell = mEnv.GrainSize * mEnv.MaxGrn;

  //clog << "Cell number = " << mEnv.MaxCell << ", Grain number = " << mEnv.MaxGrn << endl;
  cinfo << "Cell number = " << mEnv.MaxCell << ", Grain number = " << mEnv.MaxGrn << endl;

  return true;
}

// This method constructs a transition matrix based on energy grains.
//
void System::constructGrainMatrix(int msize){

  // Determine the size and location of various blocks.

  // 1. Isomers.

  //size_t msize(0) ;
  //Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
  //for (; isomeritr != m_isomers.end() ; ++isomeritr) {
  //  Molecule *isomer = isomeritr->first ;
  //  isomeritr->second = static_cast<int>(msize) ; //set location
  //  msize += isomer->getColl().get_nbasis() ;
  //}

  // 2. Pseudoisomers.

  Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin() ;
  for (; pseudoIsomeritr != m_sources.end() ; ++pseudoIsomeritr) {
	pseudoIsomeritr->second = static_cast<int>(msize) ; //set location
	msize++ ;
  }

  // Allocate space for the full system collision operator.
  if (m_reactionOperator) delete m_reactionOperator;
  m_reactionOperator = new qdMatrix(msize, 0.0) ;

  // Insert collision operators to reaction operator from individual wells.
  Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
  for (isomeritr = m_isomers.begin() ; isomeritr != m_isomers.end() ; ++isomeritr) {

	Molecule *isomer = isomeritr->first ;
	int colloptrsize = isomer->getColl().getNumberOfGroupedGrains() != 0
	  ? ( isomer->getColl().isCemetery()
	  ? isomer->getColl().get_colloptrsize() - isomer->getColl().getNumberOfGroupedGrains()
	  : isomer->getColl().get_colloptrsize() - isomer->getColl().getNumberOfGroupedGrains() + 1)
	  : isomer->getColl().get_colloptrsize();
	double omega = isomer->getColl().get_collisionFrequency();
	int idx = isomeritr->second ;

	isomer->getColl().copyCollisionOperator(m_reactionOperator, colloptrsize, idx, omega/m_meanOmega) ;

  }

  // Add connecting rate coefficients.
  for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
	(*m_pReactionManager)[i]->AddReactionTerms(m_reactionOperator,m_isomers,1.0/m_meanOmega) ;
  }

}

// This is a routine to construct the big basis matrix based on the alternative basis set method.
// The full reaction operator is subject to a similarity transformation process by a set of eigenvectors.
// If there are three wells and two sources in the system, and the eigenvectors of each well under the assumption
// of the conservation of the wells are U_0, U_1 and U_2, respectively. The transformer matrix should look like
//
//        [  U_0   0    0   0   0 ]
//        [   0   U_1   0   0   0 ]
//    U = [   0    0   U_2  0   0 ]
//        [   0    0    0   1   0 ]
//        [   0    0    0   0   1 ]
//
// This transformer matrix operates on the reaction operator to give the basis matrix by doing
//
//     M'' = U^-1 M U
//
// One then needs to decide how many members of this basis matrix to include in the reduced basis matrix for
// diagonalization.
//
void System::constructBasisMatrix(void){

  // Determine the size and location of various blocks.

  // 1. Isomers.

  size_t msize(0) ;
  Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
  for (; isomeritr != m_isomers.end() ; ++isomeritr) {
	Molecule *isomer = isomeritr->first ;
	isomeritr->second = static_cast<int>(msize) ; //set location
	msize += isomer->getColl().get_nbasis() ;
  }

  // 2. Pseudoisomers.

  Reaction::molMapType::iterator pseudoIsomeritr = m_sources.begin() ;
  for (; pseudoIsomeritr != m_sources.end() ; ++pseudoIsomeritr) {
	pseudoIsomeritr->second = static_cast<int>(msize) ; //set location
	msize++ ;
  }

  // Allocate space for the reaction operator.

  if (m_reactionOperator) delete m_reactionOperator;
  m_reactionOperator = new qdMatrix(msize, 0.0) ;

  // Insert collision operators: in the contracted basis these are the eignvalues
  // of the isomer collision operators.
  for (isomeritr = m_isomers.begin() ; isomeritr != m_isomers.end() ; ++isomeritr) {

	Molecule *isomer = isomeritr->first ;
	double omega = isomer->getColl().get_collisionFrequency() ;
	int idx = isomeritr->second ;

	isomer->getColl().copyCollisionOperatorEigenValues(m_reactionOperator, idx, omega) ;
  }

  // Add rate coefficients.
  for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
	(*m_pReactionManager)[i]->AddContractedBasisReactionTerms(m_reactionOperator,m_isomers) ;
  }

  // Print out system matrix.

  //ctest << endl << "System matrix:" << endl << endl ;
  //for (size_t i(0) ; i < msize ; ++i) {
  //  for (size_t j(0) ; j < msize ; ++j) {
  //    formatFloat(ctest, (*m_reactionOperator)[i][j],  6,  15) ;
  //  }
  //  ctest << endl ;
  //}

}

bool System::calculateEquilibriumFractions(const double beta)
{ /* Consider a three well system: e.g., A <-> B <-> C where A <-> B has Keq = K1 & B <-> C has Keq = K2.
  This routine uses the fact that the normalized equilibrated system may be described
  by a 3x3 matrix and a vector which satisfy the following:
  |-K1  1   0| |A|   |0|
  | 0  -K2  1| |B| = |0|
  | 1   1   1| |C|   |1|
  The equilibrium fraction of each isomer (or pseudo isomer, in the case of a source term) may be
  obtained by inverting the matrix shown above, and taking the elements in the final column of the inverse.
  Any system, with an arbitrary number of wells and connections, may be described by such a Matrix */

  // determine the total number of isomers + sources from the m_isomers and m_sources maps
  int eqMatrixSize = int(m_isomers.size() + m_sources.size());

  // intialize the matrix which holds the system of equations that describe the equilibrium distribution
  dMatrix  eqMatrix(eqMatrixSize);

  // initialize a map of equilibrium fractions
  m_SpeciesSequence.clear();

  // loop over the number of reactions in order to assign elements to the m_SpeciesSequence map
  // and then update the corresponding matrix elements in eqMatrix

  int counter(0);   //counter keeps track of how may elements are in the m_SpeciesSequence map
  for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {  //iterate through m_reactions

	Molecule* rct;
	Molecule* pdt;
	double Keq(0.0);

	//only need eq fracs for species in isom & assoc rxns
	if ((*m_pReactionManager)[i]->isEquilibratingReaction(Keq, &rct, &pdt)){

	  int ploc, rloc ;

	  Reaction::molMapType::iterator rctitr = m_SpeciesSequence.find(rct);   //check if the reactant is in the map
	  bool rval = (rctitr != m_SpeciesSequence.end()) ;       //if the reactant isnt in the map
	  if (rval)
		rloc = rctitr->second ;        //if the reactant is in the map, get the location

	  Reaction::molMapType::iterator pdtitr = m_SpeciesSequence.find(pdt);   //check if the product is in the map
	  bool pval = (pdtitr != m_SpeciesSequence.end()) ;       //if the product isnt in the map
	  if (pval)
		ploc = pdtitr->second;        //if the product is in the map, get the location

	  if(!rval && !pval){             // if neither reactant nor product are in the m_SpeciesSequence map
		m_SpeciesSequence[rct] = counter;            // update the eqMatrix elements
		counter++ ;
		m_SpeciesSequence[pdt] = counter;
		eqMatrix[counter-1][counter-1] -= Keq;
		eqMatrix[counter-1][counter] += 1.0;
		counter++ ;
	  }
	  else if(!rval && pval){        // if reactant isnt in m_SpeciesSequence map & product is
		m_SpeciesSequence[rct] = counter;            // update the eqMatrix matrix elements
		eqMatrix[counter-1][ploc] += 1.0;
		eqMatrix[counter-1][counter] -= Keq;
		counter++ ;
	  }
	  else if(rval && !pval){        // if reactant is in m_SpeciesSequence map & product isnt
		m_SpeciesSequence[pdt] = counter;            // update the eqMatrix matrix elements
		eqMatrix[counter-1][rloc] -= Keq;
		eqMatrix[counter-1][counter] += 1.0 ;
		counter++ ;
	  }
	  else if(rval && pval){        // if both reactant & product are in m_SpeciesSequence map

		double pdtRowSum(0.0), rctRowSum(0.0);

		for(int j(0);j<counter;++j){           // calculate pdt & rct rowSums of EqMatrix to see if the rxn is redundant
		  pdtRowSum += eqMatrix[ploc][j];
		  rctRowSum += eqMatrix[rloc][j];
		}

		if(pdtRowSum!=0.0 && rctRowSum!=0.0){ // connection is redundant
		  eqMatrix[counter-1][ploc] += 1.0 ;
		  eqMatrix[counter-1][rloc] -= Keq ;
		}
		else if(rctRowSum==0.0){              // connection is not redundant, pdts lack specification
		  eqMatrix[rloc][ploc] += 1.0 ;
		  eqMatrix[rloc][rloc] -= Keq ;
		}
		else if(pdtRowSum==0.0){
		  eqMatrix[ploc][ploc] += 1.0 ;        // connection is not redundant, rcts lack specification
		  eqMatrix[ploc][rloc] -= Keq ;
		}
	  }
	}
  }

  // if counter==0 after the for loop above, then there are no equilibrating reactions (i.e., all the reactions
  // are irreversible).  In that case, the lone isomer has an equilibrium fraction of 1.  Thus, we increment
  // counter so that the 1 is added to the eqMatrix in the for loop immediately following
  if (counter==0){
	if (m_isomers.size()){
	  Molecule* rct=(m_isomers.begin())->first;
	  m_SpeciesSequence[rct] = counter;
	}
	else if (m_sources.size()){
	  Molecule* rct=(m_sources.begin())->first;
	  m_SpeciesSequence[rct] = counter;
	}
	else{
	  return false;
	}
	++counter;
  }

  for(int i=0; i < counter; ++i){         // add ones to the final row of the matrix
	eqMatrix[counter-1][i]= 1.0;
  }

  //    ctest << "matrix elements for calculating isomer equilibrium fractions:" << endl;
  //    eqMatrix.showFinalBits(counter);

  dMatrix backup(eqMatrix);  //backup EqMatrix for error reporting

  ctest << endl << "Eq fraction matrix:" << endl;
  backup.showFinalBits(counter);

  if(eqMatrix.invertGaussianJordan()){
	cerr << "Inversion of matrix for calculating Eq fractions failed.  Matrix before inversion is: ";
	backup.showFinalBits(counter);
  }

  ctest << "inverse of Eq fraction matrix:" << endl;
  eqMatrix.showFinalBits(counter);

  Reaction::molMapType::iterator itr1;

  for(itr1= m_SpeciesSequence.begin(); itr1!=m_SpeciesSequence.end(); ++itr1){  //assign Eq fraction to appropriate Molecule
	int seqMatrixLoc = itr1->second;                          //in the Eq frac map
	Molecule* key = itr1->first;
	key->getPop().setEqFraction(eqMatrix[seqMatrixLoc][counter-1]);    //set Eq fraction to last column in eqMatrix
	string speciesName = key->getName();
	ctest << "Equilibrium Fraction for " << speciesName << " = " << key->getPop().getEqFraction() << endl;
  }
  return true;
}

void System::diagReactionOperator(const MesmerFlags &mFlags, const MesmerEnv &mEnv, const int precision, PersistPtr ppAnalysis)
{
  // Allocate space for eigenvalues.
  const size_t smsize = m_reactionOperator->size() ;
  m_eigenvalues.clear();
  m_eigenvalues.resize(smsize, 0.0);
  if (m_eigenvectors) delete m_eigenvectors;
  m_eigenvectors = new qdMatrix(smsize, 0.0) ;

  // This block prints Reaction Operator before diagonalization
  if (mFlags.printReactionOperatorNum){
	ctest << "Reaction operator --- ";
	printReactionOperator(mFlags);
  }

  //-------------------------------------------------------------
  // diagonalize the whole matrix
  switch (precision){
	case 0: // diagonalize in double
	  {
		dMatrix dDiagM(smsize);
		for ( size_t i = 0 ; i < smsize ; ++i )
		  for ( size_t j = 0 ; j < smsize ; ++j )
			dDiagM[i][j] = to_double((*m_reactionOperator)[i][j]) ;
		vector<double>  dEigenValue(smsize, 0.0);
		dDiagM.diagonalize(&dEigenValue[0]) ;
		for ( size_t i = 0 ; i < smsize ; ++i )
		  m_eigenvalues[i] = dEigenValue[i];
		for ( size_t i = 0 ; i < smsize ; ++i )
		  for ( size_t j = 0 ; j < smsize ; ++j )
			(*m_eigenvectors)[i][j] = dDiagM[i][j] ;
		break;
	  }
	case 1: // diagonalize in double double
	  {
		ddMatrix ddDiagM(smsize);
		for ( size_t i = 0 ; i < smsize ; ++i )
		  for ( size_t j = 0 ; j < smsize ; ++j )
			ddDiagM[i][j] = to_dd_real((*m_reactionOperator)[i][j]) ;
		vector<dd_real> ddEigenValue(smsize, 0.0);
		ddDiagM.diagonalize(&ddEigenValue[0]) ;
		for ( size_t i = 0 ; i < smsize ; ++i )
		  m_eigenvalues[i] = ddEigenValue[i];
		for ( size_t i = 0 ; i < smsize ; ++i )
		  for ( size_t j = 0 ; j < smsize ; ++j )
			(*m_eigenvectors)[i][j] = ddDiagM[i][j] ;
		break;
	  }
	default: // diagonalize in quad double
	  {
		(*m_eigenvectors) = (*m_reactionOperator) ;
		m_eigenvectors->diagonalize(&m_eigenvalues[0]) ;
	  }

  }
  // diagonalize the whole matrix
  //-------------------------------------------------------------

  // This block prints Eigenvectors
  if (mFlags.printReactionOperatorNum){
	ctest << "Eigenvectors --- ";
	stringstream os;
	printEigenvectors(mFlags, os);
	ctest << os.str();
  }

  size_t numberStarted = 0;
  size_t numberPrinted = smsize; // Default prints all of the eigenvalues
  if (mFlags.printEigenValuesNum > 0 && mFlags.printEigenValuesNum <= int(smsize)){ //at least prints 1 eigenvalue
	numberPrinted = mFlags.printEigenValuesNum;
	numberStarted = smsize - mFlags.printEigenValuesNum;
  }

  PersistPtr ppEigenList = ppAnalysis->XmlWriteElement("me:eigenvalueList");
  ppEigenList->XmlWriteAttribute("number",toString(smsize));
  ppEigenList->XmlWriteAttribute("selection",toString(mFlags.printEigenValuesNum));//TODO improve this
  ctest << "\nTotal number of eigenvalues = " << smsize << endl;
  ctest << "Eigenvalues\n{\n";
  for (size_t i = numberStarted ; i < smsize; ++i) {
	qd_real tmp = (mEnv.useBasisSetMethod)? m_eigenvalues[i] : m_eigenvalues[i] * m_meanOmega ;
	formatFloat(ctest, tmp, 6, 15) ;
	ctest << endl ;
	ppEigenList->XmlWriteValueElement("me:eigenvalue", to_double(tmp), 6);
  }
  ctest << "}\n";


}

void System::printReactionOperator(const MesmerFlags& mFlags)
{
  const int smsize = int(m_reactionOperator->size()) ;

  switch (mFlags.printReactionOperatorNum)
  {
  case -1:
	ctest << "Printing all (" << smsize << ") columns/rows of the Reaction Operator:\n";
	(*m_reactionOperator).showFinalBits(smsize, mFlags.print_TabbedMatrices);
	break;
  case -2:
	ctest << "Printing final 1/2 (" << smsize/2 << ") columns/rows of the Reaction Operator:\n";
	(*m_reactionOperator).showFinalBits(smsize/2, mFlags.print_TabbedMatrices);
	break;
  case -3:
	ctest << "Printing final 1/3 (" << smsize/3 << ") columns/rows of the Reaction Operator:\n";
	(*m_reactionOperator).showFinalBits(smsize/3, mFlags.print_TabbedMatrices);
	break;
  default: // the number is either smaller than -3 or positive
	if (abs(mFlags.printReactionOperatorNum) > smsize){
	  ctest << "Printing all (" << smsize << ") columns/rows of the Reaction Operator:\n";
	  (*m_reactionOperator).showFinalBits(smsize, mFlags.print_TabbedMatrices);
	}
	else{
	  ctest << "Printing final " << abs(mFlags.printReactionOperatorNum) << " columns/rows of the Reaction Operator:\n";
	  (*m_reactionOperator).showFinalBits(abs(mFlags.printReactionOperatorNum), mFlags.print_TabbedMatrices);
	}
  }
}

void System::printEigenvectors(const MesmerFlags& mFlags, std::ostream& os)
{
  const size_t smsize = m_eigenvectors->size() ;

  switch (mFlags.printReactionOperatorNum)
  {
  case -1:
	os << "Printing all (" << smsize << ") columns/rows of the eigenvectors:\n";
	(*m_eigenvectors).showFinalBits(smsize, mFlags.print_TabbedMatrices);
	break;
  case -2:
	os << "Printing final 1/2 (" << smsize/2 << ") columns/rows of the eigenvectors:\n";
	(*m_eigenvectors).showFinalBits(smsize/2, mFlags.print_TabbedMatrices);
	break;
  case -3:
	os << "Printing final 1/3 (" << smsize/3 << ") columns/rows of the eigenvectors:\n";
	(*m_eigenvectors).showFinalBits(smsize/3, mFlags.print_TabbedMatrices);
	break;
  default: // the number is either smaller than -3 or positive
	if (abs(mFlags.printReactionOperatorNum) > int(smsize)){
	  os << "Printing all (" << smsize << ") columns/rows of the eigenvectors:\n";
	  (*m_eigenvectors).showFinalBits(smsize, mFlags.print_TabbedMatrices);
	}
	else{
	  os << "Printing final " << abs(mFlags.printReactionOperatorNum) << " columns/rows of the eigenvectors:\n";
	  (*m_eigenvectors).showFinalBits(abs(mFlags.printReactionOperatorNum), mFlags.print_TabbedMatrices);
	}
  }
}  

bool System::timeEvolution(MesmerFlags& mFlags, PersistPtr ppPopList)
{
  ErrorContext c(__FUNCTION__);
  int smsize = int(m_eigenvectors->size());

  if (!produceEquilibriumVector()){
	cerr << "Calculation of equilibrium vector failed.";
	return false;
  }

  vector<double> n_0(smsize, 0.); // initial distribution
  if (!produceInitialPopulationVector(n_0)){
	cerr << "Calculation of initial conditions vector failed.";
	return false;
  }

  //Cut short if species profiles not needed
  if(!mFlags.speciesProfileEnabled)
	return true;

  // |n_0> = F^(-1)*|n_0>
  for (int j = 0; j < smsize; ++j) {
	n_0[j] /= to_double(m_eqVector[j]) ;
  }
  // Converts the initial population vector into Boltzmann weighted population vector.
  // All transitions in the reaction matrix are Boltzmann weighted for symmetry.

  // |n_0> is the initial populations of the grains for all species
  // |n_t> = U exp(Lamda t) U^-1 |n_0>
  // |r_0> = U^-1 |n_0>
  vector<double> r_0(smsize, 0.);

  double shortestTime = 0.;
  // set the default maximum evolution time
  if (mFlags.shortestTimeOfInterest < 1.0e-20 || mFlags.shortestTimeOfInterest > 1.0)
	shortestTime = 1.0e-11;
  else
	shortestTime = mFlags.shortestTimeOfInterest;

  double maxEvoTime = 0.;
  // set the default maximum evolution time
  if (mFlags.maxEvolutionTime <= 0.001 || mFlags.maxEvolutionTime > 1.0e8)
	maxEvoTime = 1.2e5;
  else
	maxEvoTime = mFlags.maxEvolutionTime;

  // Calculates the time points
  vector<double> timePoints;
  for (int i = 0; i <= 300; ++i){
	double thetime = pow(10., static_cast<double>(i) / 10. - 20.);
	if (thetime < shortestTime) continue;
	if (thetime > maxEvoTime) break;
	timePoints.push_back(thetime);
  }

  //Initialises dt vector for calculating product yields
  vector<double> dt(timePoints.size()-1,0.0);
  dt[0] = timePoints[0];
  for (int i = 1; i < int(dt.size()); ++i){
	dt[i] = timePoints[i] - timePoints[i-1];
  }


  dMatrix totalEigenVecs(smsize); // copy full eigenvectors of the system
  for ( int i = 0 ; i < smsize ; ++i )
	for ( int j = 0 ; j < smsize ; ++j )
	  totalEigenVecs[i][j] = to_double((*m_eigenvectors)[i][j]);


  for (int i = 0; i < smsize; ++i) {
	double sum = 0.;
	for (int j = 0; j < smsize; ++j) {
	  sum += n_0[j] * totalEigenVecs[j][i];
	}
	r_0[i] = sum;  // now |r_0> = V^(T)*|init> = U^(-1)*|n_0>
	// Times the initial population with the inverse of the eigenvector
	// which converts the populations into the "decay modes" domain.
  }

  for (int i = 0; i < smsize; ++i) {
	double tmp = to_double(m_eqVector[i]);
	for (int j = 0; j < smsize; ++j) {
	  totalEigenVecs[i][j] *= tmp;
	}
  }

  const int maxTimeStep = int(dt.size());
  db2D grnProfile(smsize, maxTimeStep); // numbers inside the parentheses are dummies
  vector<double> work2(smsize, 0.);

  for (int timestep = 0; timestep < maxTimeStep; ++timestep){
	double numColl = m_meanOmega * timePoints[timestep];
	for (int j = 0; j < smsize; ++j) {
	  work2[j] = r_0[j] * exp(to_double(m_eigenvalues[j]) * numColl);
	} // now |wk2> = exp(Lambda*t)*V^(T)*|init> = exp(Lambda*t)*U^(-1)*|n_0>
	for (int j = 0; j < smsize; ++j) {
	  double sum = 0.;
	  for (int l = 0; l < smsize; ++l) {
		sum += work2[l] * totalEigenVecs[j][l];
	  }
	  grnProfile[j][timestep] = sum;
	} // now |grnProfile(t)> = |grnProfile(i)> = F*V*exp(Lambda*t)*V^(T)*|init> = U*exp(Lambda*t)*U^(-1)*|n_0>
  }

  //------------------------------
  // print grained species profile
  if (mFlags.grainedProfileEnabled) {
	ctest << "\nGrained species profile (the first row is time points in unit of second):\n{\n";
	for (int timestep = 0; timestep < maxTimeStep; ++timestep){
	  formatFloat(ctest, timePoints[timestep], 6,  15);
	}
	ctest << endl;
	for (int j = 0; j < smsize; ++j) {
	  for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		formatFloat(ctest, grnProfile[j][timestep], 6,  15);
	  }
	  ctest << endl;
	}
	ctest << "}\n";
  }
  //------------------------------

  ctest<<"mean collision frequency = " << m_meanOmega << "/s" << endl;

  vector<double> totalIsomerPop(maxTimeStep, 0.);
  vector<double> totalPdtPop(maxTimeStep, 0.);

  for(int timestep(0); timestep<maxTimeStep; ++timestep){
	for(int j(0);j<smsize;++j){
	  totalIsomerPop[timestep] += grnProfile[j][timestep];
	}
	double popTime = totalIsomerPop[timestep];
	if (popTime > 1.0){
	  popTime = 1.0; // correct some numerical error
	  //totalIsomerPop[timestep] = 1.0; // Not very sure if we should cover up this numerical error entirely!!?
	}
	else if (popTime < 0.0){
	  popTime = 0.0;
	  //totalIsomerPop[timestep] = 0.0; // Not very sure if we should cover up this numerical error entirely!!?
	}
	totalPdtPop[timestep] = 1.0 - popTime;
  }

  if (mFlags.speciesProfileEnabled){
	ctest << endl << "Print time dependent species and product profiles" << endl << "{" << endl;
	int numberOfSpecies = static_cast<int>(m_isomers.size() + m_sources.size() + m_sinkRxns.size());

	//---------------------------------------------------------------------------------------------
	// Need to include the cemetery states too, so loop into isomers and see how many have cemetery.
	Reaction::molMapType::iterator iposC;
	for (iposC = m_isomers.begin(); iposC != m_isomers.end(); ++iposC){  // Iterate through the isomer map
	  Molecule* isomer = iposC->first;
	  if (isomer->getColl().isCemetery()) ++numberOfSpecies;
	}
	//---------------------------------------------------------------------------------------------

	db2D speciesProfile(numberOfSpecies, maxTimeStep);
	int speciesProfileidx(0);

	ctest << setw(16) << "Timestep/s";

	vector<string> speciesNames;
	Reaction::molMapType::iterator spos;
	for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map
	  Molecule* source = spos->first ;                        // to get source profile vs t
	  ctest << setw(16) << source->getName();
	  speciesNames.push_back(source->getName());
	  int rxnMatrixLoc = spos->second;
	  for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		double gPf = grnProfile[rxnMatrixLoc][timestep];
		speciesProfile[speciesProfileidx][timestep] = gPf;
	  }
	  ++speciesProfileidx;
	}

	Reaction::molMapType::iterator ipos;
	for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
	  Molecule* isomer = ipos->first;                        // to get isomer profile vs t
	  const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
	  string isomerName = isomer->getName();
	  if (nrg){
		ctest << setw(16) << isomerName;
	  }
	  else{
		isomerName += "(+)"; // active states
		ctest << setw(16) << isomerName;
	  }
	  speciesNames.push_back(isomerName);
	  int rxnMatrixLoc = ipos->second;
	  const int colloptrsize = isomer->getColl().get_colloptrsize();
	  const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
	  if (numberGrouped == 0){
		for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		  for(int i = 0; i < colloptrsize; ++i){
			speciesProfile[speciesProfileidx][timestep] += grnProfile[i+rxnMatrixLoc][timestep];
		  }
		}
	  }
	  else{
		for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		  for(int i = 0; i < colloptrsize - numberGrouped + nrg; ++i){
			speciesProfile[speciesProfileidx][timestep] += grnProfile[i+rxnMatrixLoc][timestep];
		  }
		}
	  }
	  ++speciesProfileidx;
	}

	int pdtProfileStartIdx = speciesProfileidx;

	// Taking account of the cemetery states in all wells.
	for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
	  Molecule* isomer = ipos->first;
	  if (isomer->getColl().isCemetery()){
		vector<double> grainKdmc = isomer->getColl().get_GrainKdmc();
		string cemName = isomer->getName() + "(-)";
		ctest << setw(16) << cemName;
		speciesNames.push_back(cemName);
		int rxnMatrixLoc = ipos->second;                       // get isomer location
		double TimeIntegratedCemeteryPop(isomer->getPop().getInitCemeteryPopulation());
		for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		  for(size_t i = 0; i < grainKdmc.size(); ++i){
			speciesProfile[speciesProfileidx][timestep] += m_meanOmega * grainKdmc[i]*grnProfile[i+rxnMatrixLoc][timestep]*dt[timestep];;
		  }
		  TimeIntegratedCemeteryPop += speciesProfile[speciesProfileidx][timestep];
		  speciesProfile[speciesProfileidx][timestep]= TimeIntegratedCemeteryPop;
		}
		++speciesProfileidx;
	  }
	}

	sinkMap::iterator pos;      // iterate through sink map to get product profile vs t
	for (pos = m_sinkRxns.begin(); pos != m_sinkRxns.end(); ++pos){
	  vector<double> KofEs;                             // vector to hold sink k(E)s
	  Reaction* sinkReaction = pos->first;
	  const int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
	  vector<Molecule*> pdts;                               // in the sink reaction
	  sinkReaction->get_products(pdts);

	  int numberGrouped(0);
	  string pdtName = pdts[0]->getName();
	  if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
		KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
		pdtName += "(bim)";
		ctest << setw(16) << pdtName;
	  }
	  else{   // if the collision operator size is >1, there are k(E)s for the irreversible loss
		KofEs = sinkReaction->get_GrainKfmc();          // assign sink k(E)s, the vector size == maxgrn
		ctest << setw(16) << pdtName;
		numberGrouped = sinkReaction->get_reactant()->getColl().getNumberOfGroupedGrains();
	  }
	  speciesNames.push_back(pdtName);
	  int rxnMatrixLoc = pos->second;                       // get sink location
	  double TimeIntegratedProductPop(0.0);
	  if (numberGrouped == 0){
		for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		  for(int i = 0; i < colloptrsize; ++i){
			speciesProfile[speciesProfileidx][timestep] += KofEs[i]*grnProfile[i+rxnMatrixLoc][timestep]*dt[timestep];
		  }
		  TimeIntegratedProductPop += speciesProfile[speciesProfileidx][timestep];
		  speciesProfile[speciesProfileidx][timestep]= TimeIntegratedProductPop;
		}
		++speciesProfileidx;
	  }
	  else{
		Molecule* rctMol = pos->first->get_reactant();
		int nrg = rctMol->getColl().isCemetery() ? 0 : 1;
		for (int timestep = 0; timestep < maxTimeStep; ++timestep){
		  for(int i = 0; i < colloptrsize - numberGrouped + nrg; ++i){
			speciesProfile[speciesProfileidx][timestep] += KofEs[i + numberGrouped - nrg]*grnProfile[i+rxnMatrixLoc][timestep]*dt[timestep];
		  }
		  TimeIntegratedProductPop += speciesProfile[speciesProfileidx][timestep];
		  speciesProfile[speciesProfileidx][timestep]= TimeIntegratedProductPop;
		}
		++speciesProfileidx;
	  }
	  KofEs.clear();
	}

	if (pdtProfileStartIdx < speciesProfileidx){
	  for(int timestep = 0; timestep < maxTimeStep; ++timestep){    // normalize product profile to account for small
		double normConst(0.0);                          // numerical errors in TimeIntegratedProductPop
		double pdtYield(0.0);
		for(int i(pdtProfileStartIdx); i<speciesProfileidx; ++i){   // calculate normalization constant
		  pdtYield += speciesProfile[i][timestep];
		}
		normConst = totalPdtPop[timestep] / pdtYield;
		for(int i(pdtProfileStartIdx); i<speciesProfileidx; ++i){   // apply normalization constant
		  speciesProfile[i][timestep] *= normConst;
		}
	  }
	}

	//Write to ctest and XML
	ctest << setw(16)<< "totalIsomerPop" << setw(16)<< "totalPdtPop"  << endl;
	for(int timestep = 0; timestep < maxTimeStep; ++timestep){
	  ctest << setw(16) << timePoints[timestep];
	  PersistPtr ppPop =  ppPopList->XmlWriteElement("me:population");
	  ppPop->XmlWriteAttribute("time", toString(timePoints[timestep]));
	  ppPop->XmlWriteAttribute("logTime", toString(log10(timePoints[timestep])));
	  for(int i(0); i<speciesProfileidx; ++i){
		ctest << setw(16) << speciesProfile[i][timestep];
		PersistPtr ppVal = ppPop->XmlWriteValueElement("me:pop", speciesProfile[i][timestep]);
		ppVal->XmlWriteAttribute("ref", speciesNames[i]);
	  }
	  ctest << setw(16) << totalIsomerPop[timestep] << setw(16) << totalPdtPop[timestep] << endl;
	}
	ctest << "}" << endl;
  }
  return true;
}

//
// The vector produced by this function contains the sqrt of
// the normalized equilibrium distribution.
//
bool System::produceEquilibriumVector()
{

  m_eqVector.clear();
  m_eqVector.resize(m_reactionOperator->size());

  Reaction::molMapType::iterator spos;
  for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // Iterate through the source map to get
	Molecule* source = spos->first;                                 // the equilibrum fractions.
	int rxnMatrixLoc = spos->second;
	qd_real eqFrac = source->getPop().getEqFraction();
	m_eqVector[rxnMatrixLoc] = sqrt(eqFrac) ;
  }

  Reaction::molMapType::iterator ipos;
  for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // Iterate through the isomer map
	Molecule* isomer = ipos->first;                                 // to get the equilibrium fractions.
	int rxnMatrixLoc = ipos->second;
	qd_real eqFrac = isomer->getPop().getEqFraction();
	const int colloptrsize = isomer->getColl().get_colloptrsize();
	const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
	const int nrg = isomer->getColl().isCemetery() ? 0 : 1 ;
	vector<double> boltzFrac;
	isomer->getColl().normalizedGrnBoltzmannDistribution(boltzFrac, colloptrsize, numberGrouped);
	if (numberGrouped == 0) {
	  for(int i(0);i<colloptrsize;++i){
		m_eqVector[rxnMatrixLoc + i] = sqrt(eqFrac * qd_real(boltzFrac[i]) ) ;
	  }
	}
	else{
	  for(int i(1-nrg), j(0);i<colloptrsize - numberGrouped + 1;++i, ++j){
		m_eqVector[rxnMatrixLoc + j] = sqrt(eqFrac * qd_real(boltzFrac[i]) ) ;
	  }
	}
  }
  return true;
}

bool System::produceInitialPopulationVector(vector<double>& n_0){

  double populationSum = 0.0;

  Reaction::molMapType::iterator ipos;
  for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){  // iterate through isomer map
	Molecule* isomer = ipos->first;                        // to get isomer initial populations
	populationSum += isomer->getPop().getInitPopulation();
  }

  Reaction::molMapType::iterator spos;
  for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  // iterate through source map to get
	Molecule* source = spos->first;                         // source initial populations
	populationSum += source->getPop().getInitPopulation();
  }

  for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
	Molecule* isomer = ipos->first;                        // get initial population of each isomer
	double initFrac = isomer->getPop().getInitPopulation();
	if (initFrac != 0.0){                                           // if isomer initial populations are nonzero
	  initFrac /= populationSum;                                    // normalize initial pop fraction
	  int rxnMatrixLoc = ipos->second;
	  const int colloptrsize = isomer->getColl().get_colloptrsize();
	  const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
	  vector<double> boltzFrac;
	  isomer->getColl().normalizedGrnBoltzmannDistribution(boltzFrac, colloptrsize, numberGrouped);
	  const int nrg = isomer->getColl().isCemetery() ? 0 : 1 ;
	  if (numberGrouped == 0){
		for (int i = 0; i < colloptrsize; ++i){
		  n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
		}
	  }
	  else{
		if (!nrg){
		  isomer->getPop().setInitCemeteryPopulation(initFrac * boltzFrac[0]);
		}
		for (int i(1-nrg), j(0); i < colloptrsize - numberGrouped + 1; ++i, ++j){
		  n_0[j + rxnMatrixLoc] = initFrac * boltzFrac[i];
		}
	  }
	}
  }

  // if there is no source term and the populationSum is still zero, set population = 1.0 for the first isomer
  int sizeSource = static_cast<int>(m_sources.size());
  if (populationSum == 0. && sizeSource == 0){
	ipos = m_isomers.begin();
	Molecule* isomer = ipos->first;
	isomer->getPop().setInitPopulation(1.0); // set initial population for the first isomer
	double initFrac = isomer->getPop().getInitPopulation();
	cinfo << "No population was assigned, and there is no source term."  << endl
	  << "Initialize a Boltzmann distribution in the first isomer." << endl;
	int rxnMatrixLoc = ipos->second;
	const int colloptrsize = isomer->getColl().get_colloptrsize();
	const int numberGrouped = isomer->getColl().getNumberOfGroupedGrains();
	vector<double> boltzFrac;
	isomer->getColl().normalizedInitialDistribution(boltzFrac, colloptrsize, numberGrouped);
	const int nrg = isomer->getColl().isCemetery() ? 0 : 1 ;
	if (numberGrouped == 0){
	  for (int i = 0; i < colloptrsize; ++i){
		n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
	  }
	}
	else{
	  for (int i = 0; i < colloptrsize - numberGrouped + nrg; ++i){
		n_0[i + rxnMatrixLoc] = initFrac * boltzFrac[i];
	  }
	}
  }

  for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){
	Molecule* source = spos->first;
	int rxnMatrixLoc = spos->second;
	if (populationSum == 0. && spos == m_sources.begin()){
	  cinfo << "No population was assigned. Initialize the first source term to 1.0." << endl;
	  n_0[rxnMatrixLoc] = 1.0;
	}else{
	  double initFrac = source->getPop().getInitPopulation() / populationSum;
	  n_0[rxnMatrixLoc] = initFrac;
	}
  }

  return true;
}

bool System::BartisWidomPhenomenologicalRates(qdMatrix& mesmerRates, MesmerFlags& mFlags, PersistPtr ppList)
{
  // Constants.
  const size_t smsize   = m_eigenvectors->size() ;
  const size_t nchem    = m_isomers.size() + m_sources.size() ;  // number of isomers+pseudoisomers
  const size_t nchemIdx = smsize - nchem ;       // Location of chemically significant eigenvalues & vectors
  const size_t nsinks   = m_sinkRxns.size() ;    // Number of Sinks.

  ctest << "\nBartis Widom eigenvalue/eigenvector analysis\n" << endl ;
  ctest << "Number of sinks in this system: " << nsinks << endl;

  if(nsinks > 0){
	ctest << "\nThere should be " << nchem << " chemically significant eigenvalues (CSEs)" << endl;
  } else {
	ctest << "\nThere should be 1 zero eigenvalue (zero within numerical precision) and " << nchem-1
	      << " chemically significant eigenvalues (CSEs)" << endl;
  }

  //
  // If there are no sinks, replace the equilibrium vector with the eigenvector whose
  // associated eigenvalue is zero, as this is a consistent estimate of the equilibrium 
  // with respect to the other eigenvalues. Also, as the system is conservative, set the 
  // smallest eigenvalue explicitly to zero.
  //
  if (nsinks < 1) {
	m_eigenvalues[smsize-1] = 0.0 ;
	for(size_t i(0) ; i<smsize ; ++i){
	  m_eqVector[i] = (*m_eigenvectors)[i][smsize-1];
	}
  }

  //
  // Construct assymmetric eigenvectors reuquied for the z matrix.
  //
  qdMatrix assymInvEigenVec(smsize);   // U^(-1)
  qdMatrix assymEigenVec(smsize);      // U
  for(size_t i(0) ; i<smsize ; ++i){
	qd_real tmp = m_eqVector[i];
	qd_real sm(0) ;
	for(size_t j(0) ; j<smsize ; ++j){
	  assymInvEigenVec[j][i] = (*m_eigenvectors)[i][j]/tmp ;          //calculation of U^(-1) = (FV)^-1 = V^T * F^-1
	  assymEigenVec[j][i] = m_eqVector[j] * (*m_eigenvectors)[j][i] ; //calculation of U = FV
	  sm += assymEigenVec[j][i] ;
	}
  }

  //------------------------- TEST block ----------------------------------------
  for(size_t i(nchemIdx) ; i<smsize ; ++i){         // multiply U*U^(-1) for testing
	qd_real test = 0.0;
	for(size_t j(nchemIdx) ; j<smsize ; ++j){
	  qd_real sm = 0.0;
	  for(size_t k(0) ; k<smsize ; ++k){
		sm += assymEigenVec[i][k] * assymInvEigenVec[k][j];
	  }
	  test += sm;
	}
	if( test < 0.999 || test > 1.001)      // test that U*U^(-1) = 1
	  ctest << "row " << i << " of the U*U^(-1) matrix does not equal unity. It sums to " << test << endl;
  }
  //------------------------- TEST block ----------------------------------------
  if (!m_Flags.rateCoefficientsOnly){
	qdMatrix Z_matrix(nchem);  // definitions of Y_matrix and Z_matrix taken from PCCP 2007(9), p.4085
	qdb2D Y_matrix;
	Reaction::molMapType::iterator ipos;  // set up an iterator through the isomer map
	Reaction::molMapType::iterator spos;  // set up an iterator through the source map
	sinkMap::iterator sinkpos;           // set up an iterator through the irreversible rxn map

	// check the separation between chemically significant eigenvalues (CSEs)
	// and internal energy relaxation eigenvalues (IEREs); if it's not good, print a warning

	const double last_CSE   = (to_double(m_eigenvalues[nchemIdx]))* m_meanOmega;
	const double first_IERE = (to_double(m_eigenvalues[nchemIdx-1]))* m_meanOmega;
	const double CSE_IERE_separation = to_double(m_eigenvalues[nchemIdx]/m_eigenvalues[nchemIdx-1]);
	if(CSE_IERE_separation > 0.1){
	  stringstream ss1 ;
	  ss1 << "\nWarning: CSEs not well separated from internal energy relaxation eigenvals (IEREs)" << endl;
	  ss1 << "\nThe last CSE = " << last_CSE << " and the first IERE = " << first_IERE << endl;
	  ss1 << "(last CSE)/(first IERE) ratio = " << CSE_IERE_separation << ", which is less than an order of magnitude" << endl;
	  ss1 << "\nResults obtained from Bartis Widom eigenvalue-vector analysis may be unreliable" << endl;
	  ctest << ss1.str() ;
	  ppList->XmlWriteValueElement("me:warning", ss1.str());
	}

	int numberOfCemeteries(0); // initialize the number of cemeteries to zero
	for(size_t i(0); i<nchem; ++i){

	  numberOfCemeteries = 0; // re-initialize for every nchem calculation.

	  // Calculate Z matrix elements for all the isomers in the system.

	  for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
		qd_real sm(0.0) ; 
		Molecule* isomer = ipos->first;
		const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
		const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
		int colloptrsize = isomer->getColl().get_colloptrsize() ;       // get colloptrsize for isomer
		colloptrsize += (numberGroupedGrains) ? nrg - numberGroupedGrains : 0 ;
		int rxnMatrixLoc = ipos->second + colloptrsize - 1 ;            // get location for isomer in the rxn matrix
		int seqMatrixLoc = m_SpeciesSequence[isomer];                   // get sequence position for isomer
		for(int j(0) ; j<colloptrsize ; ++j){
		  sm += assymEigenVec[rxnMatrixLoc-j][nchemIdx+i];
		}
		Z_matrix[seqMatrixLoc][i] = sm;
	  }

	  // Calculate Z_matrix matrix elements for all sources in the system.

	  for (spos = m_sources.begin(); spos != m_sources.end(); ++spos){  
		Molecule* pPseudoIsomer = spos->first ;
		const int rxnMatrixLoc = spos->second;
		const int seqMatrixLoc = m_SpeciesSequence[pPseudoIsomer];
		Z_matrix[seqMatrixLoc][i] = assymEigenVec[rxnMatrixLoc][nchemIdx+i];
	  }

	  // Calculate Y_matrix elements for sinks.

	  if(nsinks) {
		for(sinkpos=m_sinkRxns.begin(); sinkpos!=m_sinkRxns.end(); ++sinkpos){
		  qd_real sm = 0.0;
		  vector<double> KofEs;                                         // vector to hold sink k(E)s
		  vector<double> KofEsTemp;                                     // vector to hold sink k(E)s
		  Reaction* sinkReaction = sinkpos->first;
		  int colloptrsize = sinkReaction->getRctColloptrsize();  // get collisionoptrsize of reactant
		  if(colloptrsize == 1){  // if the collision operator size is 1, there is one canonical loss rate coefficient
			KofEs.push_back(sinkReaction->get_fwdGrnCanonicalRate());
			KofEsTemp.push_back(KofEs[0]);
		  } else {                   // if the collision operator size is >1, there are k(E)s for the irreversible loss
			KofEs = sinkReaction->get_GrainKfmc();                      // assign sink k(E)s, the vector size == maxgrn
			Molecule* isomer = sinkReaction->get_reactant();
			const int nrg = isomer->getColl().isCemetery() ? 0 : 1;
			const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();

			// DO NOT MOVE THIS SECTION --- INDEX SENSITIVE
			int ll = (numberGroupedGrains != 0) ? numberGroupedGrains - nrg : 0 ;
			for (int i(ll) ; i < colloptrsize ; ++i)
			  KofEsTemp.push_back(KofEs[i]);
			// DO NOT MOVE THIS SECTION --- INDEX SENSITIVE

			colloptrsize -= ll ;
		  }
		  int rxnMatrixLoc = sinkpos->second;                               // get sink location
		  int seqMatrixLoc = m_SinkSequence[sinkReaction];                  // get sink sequence position
		  for(int j(0);j<colloptrsize;++j){
			sm += assymEigenVec[rxnMatrixLoc+j][nchemIdx+i] * KofEsTemp[j];
		  }
		  Y_matrix[seqMatrixLoc][i] = sm;
		  KofEs.clear();
		}
	  }

	  // calculate Y_matrix matrix elements for cemetery states
	  for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
		Molecule* isomer = ipos->first;
		if (isomer->getColl().isCemetery()){ // if it is in cemetery state
		  qd_real sm = 0.0;
		  vector<double> KofEs;                                         // vector to hold sink k(E)s
		  KofEs = isomer->getColl().get_GrainKdmc();
		  const int numberGroupedGrains = isomer->getColl().getNumberOfGroupedGrains();
		  const int colloptrsize =  isomer->getColl().get_colloptrsize() - numberGroupedGrains;
		  // get colloptrsize for isomer
		  int rxnMatrixLoc = ipos->second;                                // get location for isomer in the rxn matrix
		  int seqMatrixLoc = int(nsinks) + numberOfCemeteries;      // get sequence position for isomer
		  for(int j(0);j<colloptrsize;++j){
			sm += assymEigenVec[rxnMatrixLoc+j][nchemIdx+i] * KofEs[j] * m_meanOmega;
		  }
		  Y_matrix[seqMatrixLoc][i] = sm;
		  KofEs.clear();
		  ++numberOfCemeteries;
		}
	  }
	}

	// Print out Y_matrix for testing.
	if (nsinks + numberOfCemeteries){
	  ctest << "Y_matrix:" << endl;
	  Y_matrix.print((int)(nsinks) + numberOfCemeteries, (int)(m_SpeciesSequence.size())); 
	}

	qdMatrix Zinv(Z_matrix) ;
	if (nsinks + numberOfCemeteries) {

	  // Apply standard inversion method.

	  if(Zinv.invertGaussianJordan()){
		cerr << "Inversion of Z_matrix failed.  Matrix before inversion is: ";
		Z_matrix.showFinalBits(nchem);
	  }

	} else {

	  // Apply Gram-Schmit orthogonalization in order to invert the matrix.
	  // This imposes detailed balance at the macroscopic level.
	  //
	  // SHR 25/Apr/2010 : It remains unclear that this is correct at the time
	  // of writting, however for some systems it is difficult to realize mass
	  // conservation without it.

	  // Decompose the reduced eigenvector matrix.

	  qdMatrix Fr(nchem), Fr_inv(nchem) ;
	  for(size_t i(0) ; i<nchem ; ++i){
		Fr[i][i]     = sqrt(Z_matrix[i][nchem-1]) ;
		Fr_inv[i][i] = 1.0/Fr[i][i] ;
	  }

	  qdMatrix Er = Fr_inv * Z_matrix ;

	  // Orthogonalize the reduced symmetric eigenvectro matrix.

	  Er.GramSchimdt(nchem - 1) ;

	  Z_matrix = Fr * Er ;

	  // Transpose the orthonormal matrix and form inverse.

	  Er.Transpose() ;

	  Zinv = Er * Fr_inv ;

	}

	ctest << "\nZ_matrix: ";
	Z_matrix.showFinalBits(nchem, true);

	ctest << endl << "Z_matrix^(-1):" << endl;
	Zinv.showFinalBits(nchem, true);

	qdMatrix Zidentity = Z_matrix * Zinv ;

	ctest << "\nZ_matrix * Z_matrix^(-1) [Identity matrix]:" << endl;
	Zidentity.showFinalBits(nchem, true);

	// Construct phenomenological rate coefficient matrix.

	qdMatrix Egv(nchem) ;
	for (size_t i(0) ; i<nchem ; ++i){
	  Egv[i][i] = m_eigenvalues[nchemIdx+i] * m_meanOmega ; 
	} 
	qdMatrix Kr = Z_matrix * Egv * Zinv ;

	ctest << "\nKr matrix:" << endl;
	Kr.showFinalBits(nchem, true);       // print out Kr_matrix

	// Construct loss matrix.

	qdb2D Kp;
	if (nsinks > 0) {
	  for(size_t i(0); i != nsinks + numberOfCemeteries; ++i){    // calculate Kp (definition taken from PCCP 2007(9), p.4085)
		for(size_t j(0);j<nchem;++j){
		  qd_real sm = 0.0;
		  for(size_t k(0);k<nchem;++k){
			sm += Y_matrix[i][k] * Zinv[k][j];
		  }
		  Kp[i][j] = sm;
		}
	  }
	  ctest << "\nKp matrix:" << endl;    // print out Kp_matrix
	  Kp.print(nsinks + numberOfCemeteries, m_SpeciesSequence.size());
	}

	// Write out phenomenological rate coefficients.
	PrintPhenomenologicalRates(Kr, Kp, numberOfCemeteries, mFlags, ppList) ;

	mesmerRates = Kr;
  }
  return true;    

}

// Write out phenomenological rate coefficients.
bool System::PrintPhenomenologicalRates(qdMatrix& Kr, qdb2D& Kp, int numberOfCemeteries, MesmerFlags& mFlags, PersistPtr ppList) {

  Reaction::molMapType::iterator ipos;  // set up an iterator through the isomer map

  ctest << "\nFirst order & pseudo first order rate coefficients for loss rxns:\n{\n";
  Reaction::molMapType::iterator lossitr, rctitr, pdtitr;

  stringstream puSymbols;
  stringstream puNumbers;
  // print pseudo 1st order k loss for isomers
  for(lossitr=m_SpeciesSequence.begin(); lossitr!=m_SpeciesSequence.end(); ++lossitr){
	Molecule* iso = lossitr->first;
	int losspos = lossitr->second;
	string isomerName = iso->isCemetery() ? iso->getName() + "(+)" : iso->getName();
	ctest << isomerName << " loss = " << Kr[losspos][losspos] << endl;
	PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderLoss", to_double(Kr[losspos][losspos]));
	ppItem->XmlWriteAttribute("ref", isomerName);

	puNumbers << Kr[losspos][losspos] << "\t";
	if (m_punchSymbolGathered == false){
	  puSymbols << isomerName << " loss\t";
	}
  }
  ctest << "}\n";

  if(m_SpeciesSequence.size()>1){
	ctest << "\nFirst order & pseudo first order rate coefficients for isomerization rxns:\n{\n";

	// print pseudo first order connecting ks
	for (rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
	  Molecule* rct = rctitr->first;
	  string rctName = rct->isCemetery() ? rct->getName() + "(+)" : rct->getName();
	  int rctpos = rctitr->second;
	  for (pdtitr=m_SpeciesSequence.begin(); pdtitr!=m_SpeciesSequence.end(); ++pdtitr){
		Molecule* pdt = pdtitr->first;
		string pdtName = pdt->isCemetery() ? pdt->getName() + "(+)" : pdt->getName();
		int pdtpos = pdtitr->second;
		if(rctpos != pdtpos){
		  ctest << rctName << " -> " << pdtName << " = " << Kr[pdtpos][rctpos] << endl;

		  PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", to_double(Kr[pdtpos][rctpos]));
		  ppItem->XmlWriteAttribute("fromRef", rctName);
		  ppItem->XmlWriteAttribute("toRef",   pdtName);
		  ppItem->XmlWriteAttribute("reactionType", "isomerization");
		}

		puNumbers << Kr[pdtpos][rctpos] << "\t";
		if (m_punchSymbolGathered == false){
		  puSymbols << rctName << " -> " << pdtName << "\t";
		}
	  }
	}
	ctest << "}\n";
  }

  if(m_sinkRxns.size()!=0){
	ctest << "\nFirst order & pseudo first order rate coefficients for irreversible rxns:\n{\n";
	sinkMap::iterator sinkitr;

	for(sinkitr=m_SinkSequence.begin(); sinkitr!=m_SinkSequence.end(); ++sinkitr){
	  Reaction* sinkReaction = sinkitr->first;          // get Irreversible Rxn
	  int sinkpos = m_SinkSequence[sinkReaction];                   // get products & their position
	  vector<Molecule*> pdts;
	  sinkReaction->get_products(pdts);
	  string pdtsName = pdts[0]->getName();
	  if (pdts.size() == 2) {pdtsName += + "+"; pdtsName += pdts[1]->getName();}
	  for(rctitr=m_SpeciesSequence.begin(); rctitr!=m_SpeciesSequence.end(); ++rctitr){
		Molecule* rcts = rctitr->first;     // get reactants & their position
		int rctpos = rctitr->second;
		if(sinkReaction->getRctColloptrsize()==1){
		  ctest << rcts->getName() << " -> "  << pdtsName << "(bim) = " << Kp[sinkpos][rctpos] << endl;
		  puNumbers << Kp[sinkpos][rctpos] << "\t";
		  if (m_punchSymbolGathered == false){
			puSymbols << rcts->getName() << " -> " << pdtsName << "(bim)\t";
		  }
		}
		else{
		  string rctName = rcts->isCemetery() ? rcts->getName() + "(+)" : rcts->getName();

		  ctest << rctName << " -> "  << pdtsName << " = " << Kp[sinkpos][rctpos] << endl;

		  PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", to_double(Kp[sinkpos][rctpos]));
		  ppItem->XmlWriteAttribute("fromRef", rctName);
		  ppItem->XmlWriteAttribute("toRef",   pdtsName);
		  ppItem->XmlWriteAttribute("reactionType", "irreversible");
		  puNumbers << Kp[sinkpos][rctpos] << "\t";
		  if (m_punchSymbolGathered == false){
			puSymbols << rctName << " -> " << pdtsName << "\t";
		  }
		}
	  }
	}
	ctest << "}\n\n";
  }

  Reaction::molMapType::iterator speciesitr;
  int tempNumCemetery(0);
  if(numberOfCemeteries){
	ctest << "\nFirst order & pseudo first order rate coefficients for cemetery states:\n{\n";
	for (ipos = m_isomers.begin(); ipos != m_isomers.end(); ++ipos){
	  Molecule* isomer = ipos->first;
	  if (isomer->isCemetery()){ // if it is in cemetery state
		string cemName = isomer->getName() + "(-)";
		for (speciesitr=m_SpeciesSequence.begin(); speciesitr!=m_SpeciesSequence.end(); ++speciesitr){
		  Molecule* rct = speciesitr->first;
		  string rctName = rct->isCemetery() ? rct->getName() + "(+)" : rct->getName();
		  int rctpos = speciesitr->second;
		  ctest << rctName << " -> " << cemName << " = " << Kp[m_sinkRxns.size()+tempNumCemetery][rctpos] << endl;

		  PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", to_double(Kp[m_sinkRxns.size()+tempNumCemetery][rctpos]));
		  ppItem->XmlWriteAttribute("fromRef", rctName);
		  ppItem->XmlWriteAttribute("toRef",   cemName);
		  ppItem->XmlWriteAttribute("reactionType", "deactivation");

		  puNumbers << Kp[m_sinkRxns.size()+tempNumCemetery][rctpos] << "\t";
		  if (m_punchSymbolGathered == false){
			puSymbols << rctName << " -> " << cemName << "\t";
		  }
		}
		++tempNumCemetery;
	  }
	}
	ctest << "}\n\n";
  }

  if (puSymbols.str().size()) {
	puSymbols << "\n";
	mFlags.punchSymbols = puSymbols.str();
	m_punchSymbolGathered = true;
  }

  if (puNumbers.str().size()) {
	puNumbers << "\n";
	mFlags.punchNumbers = puNumbers.str();
  }

  return true;
}

//
// Calculates the Bartis-Widom macroscopic rate coefficients, using the contracted basis set eigenvectors.
//
bool System::BartisWidomBasisSetRates(qdMatrix& mesmerRates, MesmerFlags& mFlags) {

  // Constants.
  const size_t smsize   = m_eigenvectors->size() ;
  const size_t nchem    = m_isomers.size() + m_sources.size() ;  // number of isomers+pseudoisomers
  const size_t nchemIdx = smsize - nchem ;                       // idx for chemically significant eigenvalues & vectors

  // Print out eigenvector matrix.

  //ctest << endl << "Eigenvector matrix:" << endl << endl ;
  //for (size_t i(0) ; i < smsize ; ++i) {
  //  for (size_t j(0) ; j < smsize ; ++j) {
  //    formatFloat(ctest, (*m_eigenvectors)[i][j],  6,  15) ;
  //  }
  //  ctest << endl ;
  //}

  qdMatrix Z(nchem), Zinv(nchem), Kr(nchem);

  if (m_sinkRxns.size()==0){

	//
	// Conservative system.
	//

	// 1. Isomers.

	size_t location(0) ;
	Reaction::molMapType::iterator isomeritr = m_isomers.begin() ;
	for (size_t i(0); isomeritr != m_isomers.end() ; ++i, ++isomeritr) {
	  location = isomeritr->second ;
	  for (size_t j(1); j<=nchem; ++j){
		Z[i][nchem - j] = (*m_eigenvectors)[location][smsize - j] ;
	  }
	}

	// Invert Z matrix. 

	ctest << endl << "BW coefficient matrix:" << endl << endl ;
	for (size_t i(0) ; i < nchem ; ++i) {
	  for (size_t j(0) ; j < nchem ; ++j) {
		formatFloat(ctest, Z[i][j],  6,  15) ;
		Zinv[j][i] = Z[i][j] ;
	  }
	  ctest << endl ;
	}

	// Calculate symmetric rate matrix.

	m_eigenvalues[smsize - 1] = 0.0 ;

	for (size_t i(0) ; i < nchem ; ++i) {
	  for (size_t j(0) ; j < nchem ; ++j) {
		qd_real sm = 0.0;
		for (size_t k(0) ; k < nchem ; ++k) {
		  // sm += Z[i][k] * to_double(m_eigenvalues[nchemIdx+k]) * Zinv[k][j];
		  sm += Zinv[i][k]*Z[k][j] ;
		}
		Kr[i][j] = sm ; // * m_meanOmega;
	  }
	}

	// Apply similarity transform. 

	//for (size_t i(0) ; i < nchem ; ++i) {
	//  for (size_t j(0) ; j < nchem ; ++j) {
	//    Kr[i][j] *= Z[i][nchem]/Z[j][nchem];
	//  }
	//}

	string rcm(string("Rate coefficient matrix:"));
	Kr.print(rcm, ctest) ;

  } else {

	//
	// Non-conservative system.
	//

  }

  mesmerRates = Kr;

  return true;

}

int System::getSpeciesSequenceIndex(const std::string ref)
{
  Reaction::molMapType::iterator spcitr;
  for (spcitr = m_SpeciesSequence.begin(); spcitr != m_SpeciesSequence.end(); ++spcitr)
  {
	if (ref == (spcitr->first)->getName())
	  return spcitr->second;
  }
  cerr << "No molecule named " << ref << " is available in the reaction species.";
  return -1;
}

// This method locates all sinks and determines their location in the relevant
// isomer or source map. 
void System::locateSinks()
{
  int sinkpos(0);
  m_sinkRxns.clear();
  m_SinkSequence.clear();                      
  for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {

	Reaction* pReaction = (*m_pReactionManager)[i];
	ReactionType reactionType = pReaction->getReactionType() ;

	bool Irreversible = (reactionType == IRREVERSIBLE_ISOMERIZATION || reactionType == IRREVERSIBLE_EXCHANGE || reactionType == DISSOCIATION );
	if (Irreversible && m_sinkRxns.find(pReaction) == m_sinkRxns.end()) {   
	  // Add an irreversible rxn to the map.
	  Molecule* rctnt = pReaction->get_reactant();
	  int rxnMatrixLoc ;
	  if(reactionType == IRREVERSIBLE_ISOMERIZATION || reactionType == DISSOCIATION ){
		rxnMatrixLoc = m_isomers[rctnt];
	  } else { // Irreversible exchange reaction.
		rxnMatrixLoc = m_sources[rctnt];
	  }
	  m_sinkRxns[pReaction] = rxnMatrixLoc;
	  m_SinkSequence[pReaction] = sinkpos;
	  ++sinkpos;
	}
  }

}

} //namespace

