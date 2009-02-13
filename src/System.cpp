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

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  System::System(const std::string& libraryfilename): m_pMoleculeManager(0), m_pReactionManager(0), m_pTitle(NULL){
    m_pMoleculeManager = new MoleculeManager(libraryfilename) ;
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

    m_pTitle = ppIOPtr->XmlReadValue("title", false);

    //-------------
    //Molecule List (parse this part inside Reaction)
    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    m_pMoleculeManager->set_PersistPtr(ppMolList);

    //-------------
    //Model Parameters
    PersistPtr ppParams = ppIOPtr->XmlMoveTo("me:modelParameters");
    if(ppParams)
    {
      m_Env.GrainSize = ppParams->XmlReadDouble("me:grainSize");

      m_Env.MaximumTemperature = ppParams->XmlReadDouble("me:maxTemperature",optional);
      m_Env.EAboveHill = ppParams->XmlReadDouble("me:energyAboveTheTopHill");
    }

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
      m_Flags.testRateConstantEnabled     = ppControl->XmlReadBoolean("me:testRateConstant");
      m_Flags.microRateEnabled            = ppControl->XmlReadBoolean("me:testMicroRates");
      m_Flags.grainBoltzmannEnabled       = ppControl->XmlReadBoolean("me:printGrainBoltzmann");
      m_Flags.grainDOSEnabled             = ppControl->XmlReadBoolean("me:printGrainDOS");
      m_Flags.cellDOSEnabled              = ppControl->XmlReadBoolean("me:printCellDOS");
      m_Flags.reactionOCSEnabled          = ppControl->XmlReadBoolean("me:printReactionOperatorColumnSums");
      m_Flags.kfEGrainsEnabled            = ppControl->XmlReadBoolean("me:printGrainkfE");
      m_Flags.kbEGrainsEnabled            = ppControl->XmlReadBoolean("me:printGrainkbE");
      // Both Tunnelling and Tunneling will work
      m_Flags.TunnellingCoeffEnabled      = ppControl->XmlReadBoolean("me:printTunnellingCoefficients");
      if (!m_Flags.TunnellingCoeffEnabled)
        m_Flags.TunnellingCoeffEnabled    = ppControl->XmlReadBoolean("me:printTunnellingCoefficients");
      m_Flags.cellFluxEnabled           = ppControl->XmlReadBoolean("me:printCellTransitionStateFlux");
      m_Flags.grainFluxEnabled          = ppControl->XmlReadBoolean("me:printGrainTransitionStateFlux");
      m_Flags.rateCoefficientsOnly        = ppControl->XmlReadBoolean("me:calculateRateCoefficientsOnly");
      m_Flags.useTheSameCellNumber        = ppControl->XmlReadBoolean("me:useTheSameCellNumberForAllConditions");
      m_Flags.grainedProfileEnabled       = ppControl->XmlReadBoolean("me:printGrainedSpeciesProfile");
      m_Flags.speciesProfileEnabled       = ppControl->XmlReadBoolean("me:printSpeciesProfile");
      m_Flags.viewEvents                  = ppControl->XmlReadBoolean("me:printEventsTimeStamps");
      m_Flags.allowSmallerDEDown          = ppControl->XmlReadBoolean("me:allowSmallerDeltaEDown");
      m_Flags.print_TabbedMatrices        = ppControl->XmlReadBoolean("me:printTabbedMatrices");
      m_Flags.useDOSweighedDT             = ppControl->XmlReadBoolean("me:useDOSweighedDownWardTransition");
      m_Flags.doBasisSetMethod            = ppControl->XmlReadBoolean("me:runBasisSetMethodroutines");
      if (!m_Flags.useTheSameCellNumber && m_Env.MaximumTemperature != 0.0){
        m_Flags.useTheSameCellNumber = true;
      }

      // System configuration information
      if (ppControl->XmlReadBoolean("me:runPlatformDependentPrecisionCheck")) configuration();

      if (ppControl->XmlReadBoolean("me:gridSearch")) m_Flags.searchMethod = 1;
      else if (ppControl->XmlReadBoolean("me:fitting")) m_Flags.searchMethod = 2;
      else if (ppControl->XmlReadBoolean("me:gridSearchWithPunch")) m_Flags.searchMethod = 3;
      else m_Flags.searchMethod = 0;

      if (m_Flags.grainedProfileEnabled && (m_Flags.speciesProfileEnabled || m_Flags.searchMethod)){
        cinfo << "Turn off grained species profile to prevent disk flooding." << endl;
        m_Flags.grainedProfileEnabled = false;
      }

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

      const char* txtMET = ppControl->XmlReadValue("me:MaximumEvolutionTime", optional);
      if (txtMET){
        istringstream ss(txtMET);
        ss >> m_Flags.maxEvolutionTime;
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
      txt = ppPTset->XmlReadValue("me:units",optional);
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

    int totalSteps = 1, dataPointSize = int(fitDP.size());
    for (int varID(0); varID < dataPointSize; ++varID){
      // for every dimension create a series and duplicate the serie while the indices going forward
      double lower(0.0), upper(0.0), stepsize(0.0);
      fitDP[varID].get_range(lower, upper, stepsize);
      int numSteps = int((upper - lower) / stepsize) + 1;
      totalSteps *= numSteps;
    }

    int spanSteps = 1;
    for (int varID(0); varID < dataPointSize; ++varID){
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

    ofstream punchStream(m_Flags.punchFileName.c_str());

    for (int i(0); i < totalSteps; ++i){
      double chiSquare(1000.0);

      // assign values
      for (int varID(0); varID < dataPointSize; ++varID) fitDP[varID] = gridArray[i][varID];

      // calculate
      cerr << "Parameter Grid " << calPoint;
      ctest << "Parameter Grid " << calPoint << "\n{\n";
      calculate(chiSquare);

      if (dataPointSize){
        ctest << "Parameters: ( ";
        if (m_Flags.punchSymbols.size()){
          for (int varID(0); varID < dataPointSize; ++varID){
            punchStream << "Para" << varID << "\t";
          }
          punchStream << "Temperature (K)\tNumber density\t";
          punchStream << m_Flags.punchSymbols;
          m_Flags.punchSymbols.clear();
        }
        for (int varID(0); varID < dataPointSize; ++varID){
          ctest << gridArray[i][varID] << " ";
          punchStream << gridArray[i][varID] << "\t";
        }
        punchStream << m_Env.beta << "\t" << m_Env.conc << "\t";
        punchStream << m_Flags.punchNumbers;
        m_Flags.punchNumbers.clear();

        punchStream.flush();

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
    // Controls the print-out of grain/cell DOS in each cycle (This is only for source term)
    if (m_Flags.cellDOSEnabled) m_Flags.cyclePrintCellDOS = true;
    if (m_Flags.grainDOSEnabled) m_Flags.cyclePrintGrainDOS = true;

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
    //XML output
    //Considered putting this output under each PT pair in me:conditions.
    //But doesn't work with a range of Ps or Ts. So has to have its own section.
    string comment( "Bartis-Widom Phenomenological Rates" );
    PersistPtr ppAnalysis = m_ppIOPtr->XmlWriteMainElement("me:analysis", comment);

    for (calPoint = 0; calPoint < PandTs.size(); ++calPoint){
      m_Env.beta = 1.0 / (boltzmann_RCpK * PandTs[calPoint].get_temperature()) ; //temporary statements
      m_Env.conc = PandTs[calPoint].get_concentration();
      // unit of conc: particles per cubic centimeter
      //clog << "PT Grid " << calPoint << endl;
      cinfo << "PT Grid " << calPoint << endl;
      int precision = PandTs[calPoint].get_precision();
      ctest << "PT Grid " << calPoint << " Condition: conc = " << m_Env.conc << ", temp = " << PandTs[calPoint].get_temperature();
      PersistPtr ppList = ppAnalysis->XmlWriteElement("me:rateList");
      ppList->XmlWriteAttribute("T", toString(PandTs[calPoint].get_temperature()));
      ppList->XmlWriteAttribute("conc", toString(m_Env.conc));
      ppList->XmlWriteAttribute("me:units", "s-1");
      
      switch (precision){
        case 1: ctest << ", diagonalization precision: double-double\n{\n"; break;
        case 2: ctest << ", diagonalization precision: quad-double\n{\n"; break;
        default: ctest << ", diagonalization precision: double\n{\n";
      }
      // Build collison matrix for system.
      {string thisEvent = "Build Collison Operator";
      cinfo << thisEvent << endl;}
      if (!m_pReactionManager->BuildReactionOperator(m_Env, m_Flags)){
        cerr << "Failed building system collison operator.";
        exit(1);
      }

      if (!m_pReactionManager->calculateEquilibriumFractions(m_Env.beta)){
        cerr << "Failed calculating equilibrium fractions.";
        exit(1);
      }

      // Calculate eigenvectors and eigenvalues.
      {string thisEvent = "Diagonlize the Reaction Operator";
      cinfo << thisEvent << " -- Time elapsed: " << timeElapsed << " seconds.\n";
      events.setTimeStamp(thisEvent, timeElapsed);}

      //-------------------------------
      // Reduced raction matrix operation
      //-------------------------------

      if (m_Flags.doBasisSetMethod) {
        m_pReactionManager->constructBasisMatrix();
        //dMatrix reducedMesmerRates(1);
        //m_pReactionManager->BartisWidomRatesFromBasisSetMethod(reducedMesmerRates, m_Flags, ppList);
      }

      //-------------------------------
      // Total raction matrix operation
      //-------------------------------

      // This is where the collision operator being diagonalised.
      m_pReactionManager->diagReactionOperator(m_Flags, precision) ;

      // Time steps loop
      m_pReactionManager->timeEvolution(m_Flags);

      dMatrix mesmerRates(1);
      m_pReactionManager->BartisWidomPhenomenologicalRates(mesmerRates, m_Flags, ppList);


      if (m_Flags.searchMethod){
        vector<conditionSet> expRates;
        PandTs[calPoint].get_experimentalRates(expRates);
        chiSquare += m_pReactionManager->calcChiSquare(mesmerRates, expRates);
      }

      ctest << "}\n";

    }
    //---------------------------------------------

    {string thisEvent = "Finish Calculation";
    events.setTimeStamp(thisEvent, timeElapsed);
    cinfo << endl << thisEvent << " -- Time elapsed: " << timeElapsed << " seconds.\n";
    cwarn << "In total, " << calPoint << " temperature/concentration-pressure points calculated." << endl;}

    if (m_Flags.viewEvents) cinfo << events;
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
    if(m_pTitle)
    {
      PersistPtr ppItem = ppList->XmlWriteElement("metadata");
      ppItem->XmlWriteAttribute("name", "dc:title");
      ppItem->XmlWriteAttribute("content", m_pTitle);
    }

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
      thisEvent = "Write metadata";
      timeString = events.setTimeStamp(thisEvent);
      cinfo << thisEvent << endl;
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

}//namespace
