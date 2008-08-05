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
    initializeConversionMaps();
    m_ppIOPtr = ppIOPtr;

    //-------------
    //Molecule List (parse this part inside Reaction)
    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    m_pMoleculeManager->set_PersistPtr(ppMolList);

    //-------------
    //Model Parameters
    PersistPtr ppParams = ppIOPtr->XmlMoveTo("me:modelParameters");
    if(ppParams)
    {
      const char* txt = ppParams->XmlReadValue("me:grainSize",false);
      if(txt){
        istringstream ss(txt); ss >> m_Env.GrainSize;
      }
      else{
        cinfo << "Grain size is not provided. Default grain size = " << m_Env.GrainSize << " is used";
      }

      txt = ppParams->XmlReadValue("me:maxTemperature",false);
      if(txt) {
        istringstream ss(txt);
        ss >> m_Env.MaximumTemperature;
      }

      txt = ppParams->XmlReadValue("me:energyAboveTheTopHill",false);
      if(txt) { istringstream ss(txt); ss >> m_Env.EAboveHill; }

    }

    //-------------
    //Reaction List
    PersistPtr ppReacList = ppIOPtr->XmlMoveTo("reactionList");
    if(!ppReacList)
    {
      cerr << "No reactions have been specified";
      return false;
    }
    if(!m_pReactionManager->addreactions(ppReacList, m_Env)) return false;

    //-------------
    //Reaction Conditions
    PersistPtr ppConditions = ppIOPtr->XmlMoveTo("me:conditions");
    if(!ppConditions)
    {
      cerr << "No conditions specified";
      return false;
    }
    string Bgtxt = ppConditions->XmlReadValue("me:bathGas");
    if (Bgtxt.empty()){
      cerr << "No bath gas specified";
      return false;
    }
    else{
      string molType = "bathGas";
      if(!m_pMoleculeManager->addmol(Bgtxt, molType, ppMolList, m_Env))
        return false;
      m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;
    }

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
      m_Env.testDOSEnabled              = ppControl->XmlReadBoolean("me:testDOS");
      m_Env.testRateConstantEnabled     = ppControl->XmlReadBoolean("me:testRateConstant");
      m_Env.microRateEnabled            = ppControl->XmlReadBoolean("me:testMicroRates");
      m_Env.grainBoltzmannEnabled       = ppControl->XmlReadBoolean("me:printGrainBoltzmann");
      m_Env.grainDOSEnabled             = ppControl->XmlReadBoolean("me:printGrainDOS");
      m_Env.cellDOSEnabled              = ppControl->XmlReadBoolean("me:printCellDOS");
      m_Env.collisionOCSEnabled         = ppControl->XmlReadBoolean("me:printCollisionOperatorColumnSums");
      m_Env.kfEGrainsEnabled            = ppControl->XmlReadBoolean("me:printGrainkfE");
      m_Env.kbEGrainsEnabled            = ppControl->XmlReadBoolean("me:printGrainkbE");
      // Both Tunnelling and Tunneling will work
      m_Env.TunnellingCoeffEnabled      = ppControl->XmlReadBoolean("me:printTunnellingCoefficients");
      if (!m_Env.TunnellingCoeffEnabled)
        m_Env.TunnellingCoeffEnabled    = ppControl->XmlReadBoolean("me:printTunnelingCoefficients");
      m_Env.cellTSFluxEnabled           = ppControl->XmlReadBoolean("me:printCellTransitionStateFlux");
      m_Env.grainTSFluxEnabled          = ppControl->XmlReadBoolean("me:printGrainTransitionStateFlux");
      m_Env.rateCoefficientsOnly        = ppControl->XmlReadBoolean("me:calculateRateCoefficinetsOnly");
      m_Env.useTheSameCellNumber        = ppControl->XmlReadBoolean("me:useTheSameCellNumberForAllConditions");
      m_Env.grainedProfileEnabled       = ppControl->XmlReadBoolean("me:printGrainedSpeciesProfile");
      m_Env.speciesProfileEnabled       = ppControl->XmlReadBoolean("me:printSpeciesProfile");
      if (!m_Env.useTheSameCellNumber && m_Env.MaximumTemperature != 0.0){
        m_Env.useTheSameCellNumber = true;
      }
      
      // System configuration information
      if (ppControl->XmlReadBoolean("me:runPlatformDependentPrecisionCheck")) configuration();

      if (ppControl->XmlReadBoolean("me:gridSearch")) m_Env.searchMethod = 1;
      else if (ppControl->XmlReadBoolean("me:fitting")) m_Env.searchMethod = 2;
      else m_Env.searchMethod = 0;

      const char* txtEV = ppControl->XmlReadValue("me:eigenvalues",false);
      if(txtEV) {
        istringstream ss(txtEV);
        ss >> m_Env.printEigenValuesNum;
      }
      const char* txtMET = ppControl->XmlReadValue("me:MaximumEvolutionTime");
      if (txtMET){
        istringstream ss(txtMET);
        ss >> m_Env.maxEvolutionTime;
      }
    }

    return true;
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
      txt = ppPTset->XmlReadValue("me:units");
      string this_units = txt;
      if (!txt){
        cerr << "No units provided. Default units " << default_unit << " are used.";
        this_units = default_unit;
      }

      // if user does not input any value for temperature and concentration, give a Default set of concentration and pressure
      // for simulation
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
      txt = ppPTpair->XmlReadValue("me:units");
      if (txt)
        this_units = txt;
      double this_P = default_P;
      double this_T = default_T;
      int this_precision = 0;

      txt = ppPTpair->XmlReadValue("me:P");
      if (txt){ stringstream s1(txt); s1 >> this_P; }
      txt = ppPTpair->XmlReadValue("me:T");
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

      // Set experimental conditions for chiSquare calculation
      txt = ppPTpair->XmlReadValue("me:experimentalRate", false);
      PersistPtr ppExpRate = ppPTpair->XmlMoveTo("me:experimentalRate");
      while (ppExpRate){
        double rateValue(0.0), errorValue(0.0); string ref1, ref2;
        stringstream s1(txt); s1 >> rateValue;
        txt = ppExpRate->XmlReadValue("ref1");
        stringstream s2(txt); s2 >> ref1;
        txt = ppExpRate->XmlReadValue("ref2");
        stringstream s3(txt); s3 >> ref2;
        txt = ppExpRate->XmlReadValue("error");
        stringstream s4(txt); s4 >> errorValue;
        thisPair.set_experimentalRate(ref1, ref2, rateValue, errorValue);
        ppExpRate = ppExpRate->XmlMoveTo("me:experimentalRate");
      }

      PandTs.push_back(thisPair);
      ppPTpair = ppPTpair->XmlMoveTo("me:PTpair");
    }
  }

  void System::gridSearch(void){

    // produce a grid for search
    db2D gridArray; // this array grows up freely without needing dimensions

    int totalSteps = 1;
    for (size_t varID(0); varID < fitDP.size(); ++varID){
      // for every dimension create a series and duplicate the serie while the indices going forward
      double lower(0.0), upper(0.0), stepsize(0.0);
      fitDP[varID].get_range(lower, upper, stepsize);
      int numSteps = int((upper - lower) / stepsize) + 1;
      totalSteps *= numSteps;
    }

    int spanSteps = 1;
    for (int varID(0); varID < int(fitDP.size()); ++varID){
      double lower(0.0), upper(0.0), stepsize(0.0);
      fitDP[varID].get_range(lower, upper, stepsize);
      int numSteps = int((upper - lower) / stepsize) + 1;
      int stack(0), block(0);
      while (stack < totalSteps){
        int step(0);
        while(step < spanSteps){
          double stepValue = lower + block * stepsize;
          gridArray[stack][varID] = stepValue;
          ++step;
          ++stack;
        }
        ++block;
        if (block == numSteps) block = 0;
      }
      spanSteps *= numSteps;
    }

    // TimeCount events; unsigned int timeElapsed; 
    int calPoint(0);

    for (int i(0); i < totalSteps; ++i){
      double chiSquare(1000.0);

      // assign values
      for (int varID(0); varID < int(fitDP.size()); ++varID) fitDP[varID] = gridArray[i][varID];

      // calculate
      cerr << "Parameter Grid " << calPoint;
      ctest << "Parameter Grid " << calPoint << "\n{\n";
      calculate(chiSquare);

      if (fitDP.size()){
        ctest << "Parameters: ( ";
        for (int varID(0); varID < int(fitDP.size()); ++varID) ctest << gridArray[i][varID] << " ";
        ctest << "chiSquare = " << chiSquare << " )\n}\n";
      }
      ++calPoint;
    }

  }

  // This function calls calculate() to obtain a serie of
  void System::fitting(void){

    double chiSquare = 1000.0;

    int steps(0);
    while (1){

      for (size_t obj(0); obj < fitDP.size(); ++obj){

      }

      ++steps;
      ctest << "Step " << steps << " of fitting. chiSquare = " << chiSquare << endl;
      if (chiSquare < 1 || steps > 100) break;
    }

  }

  //
  // Begin calculation.
  //
  void System::calculate(double& chiSquare)
  {
    chiSquare = 0.0; // reset the value to zero

    TimeCount events; unsigned int timeElapsed =0;

    WriteMetadata();

    // Find the highest temperature
    for (unsigned int i = 0; i < PandTs.size(); ++i){
      m_Env.MaximumTemperature = max(m_Env.MaximumTemperature, PandTs[i].get_temperature());
    }

    //---------------------------------------------
    // looping over temperatures and concentrations
    unsigned int calPoint(0);
    for (calPoint = 0; calPoint < PandTs.size(); ++calPoint){
      m_Env.beta = 1.0 / (boltzmann_RCpK * PandTs[calPoint].get_temperature()) ; //temporary statements
      m_Env.conc = PandTs[calPoint].get_concentration();
      // unit of conc: particles per cubic centimeter
      cerr << "PT Grid " << calPoint;
      int precision = PandTs[calPoint].get_precision();
      ctest << "PT Grid " << calPoint << " Condition: conc = " << m_Env.conc << ", temp = " << PandTs[calPoint].get_temperature();
      switch (precision){
        case 1: ctest << ", diagonalization precision: double-double\n{\n"; break;
        case 2: ctest << ", diagonalization precision: quad-double\n{\n"; break;
        default: ctest << ", diagonalization precision: double\n{\n";
      }
      // Build collison matrix for system.
      {string thisEvent = "Build Collison Operator";
      cinfo << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;}
      if (!m_pReactionManager->BuildSystemCollisionOperator(m_Env)){
        cerr << "Failed building system collison operator.";
      }

      if (!m_pReactionManager->calculateEquilibriumFractions(m_Env.beta)){
        cerr << "Failed calculating equilibrium fractions.";
      }

      // Calculate eigenvectors and eigenvalues.
      {string thisEvent = "Diagonlize Collision Operator";
      cinfo << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed)  << " -- Time elapsed: " << timeElapsed << " seconds.\n";}

      m_pReactionManager->diagCollisionOperator(m_Env, precision) ;

      // Time steps loop
      int timestep = 160;
      m_pReactionManager->timeEvolution(timestep, m_Env);

      dMatrix mesmerRates(1);
      m_pReactionManager->BartisWidomPhenomenologicalRates(mesmerRates);

      if (m_Env.searchMethod){
        vector<conditionSet> expRates;
        PandTs[calPoint].get_experimentalRates(expRates);
        chiSquare += m_pReactionManager->calcChiSquare(mesmerRates, expRates);
      }

      ctest << "}\n";

    }
    //---------------------------------------------

    {string thisEvent = "Finish Calculation";
    cinfo << endl << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed)  << " -- Time elapsed: " << timeElapsed << " seconds.\n";
    cinfo << "In total, " << calPoint << " temperature/concentration-pressure points calculated." << endl;}

    cinfo << events;
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
      thisEvent = "Write XML attribute";
      timeString = events.setTimeStamp(thisEvent);
      cinfo << thisEvent << " at " << timeString << endl;
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

  void System::configuration(void){
    cerr << "\nPrinting system precision configuration:";
    cerr << "Size of float = " << sizeof(float);
    cerr << "Size of double = " << sizeof(double);
    cerr << "Size of long double = " << sizeof(long double);
    cerr << "Size of double-double = " << sizeof(dd_real);
    cerr << "Size of quad-double = " << sizeof(qd_real);

    cerr << "\nEpsilon is the difference between 1 and the smallest value greater than 1 that is representable for the data type.";
    cerr << "float epsilon == " << numeric_limits<float>::epsilon() << endl;
    cerr << "double epsilon == " << numeric_limits<double>::epsilon() << endl;
    cerr << "long double epsilon == " << numeric_limits<long double>::epsilon() << endl;
    cerr << "dd_real epsilon == " << numeric_limits<dd_real>::epsilon() << endl;
    cerr << "qd_real epsilon == " << numeric_limits<qd_real>::epsilon() << endl;

    cerr << "\nfloat max == " << numeric_limits<float>::max() << endl;
    cerr << "double max == " << numeric_limits<double>::max() << endl;
    cerr << "long double max == " << numeric_limits<long double>::max() << endl;
    cerr << "dd_real max == " << numeric_limits<dd_real>::max() << endl;
    cerr << "qd_real max == " << numeric_limits<qd_real>::max() << endl;

    cerr << "\nfloat min == " << numeric_limits<float>::min() << endl;
    cerr << "double min == " << numeric_limits<double>::min() << endl;
    cerr << "long double min == " << numeric_limits<long double>::min() << endl;
    cerr << "dd_real min == " << numeric_limits<dd_real>::min() << endl;
    cerr << "qd_real min == " << numeric_limits<qd_real>::min() << endl;
  }

}//namespace
