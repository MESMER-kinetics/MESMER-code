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
      m_pMoleculeManager->addmol(Bgtxt, molType, ppMolList, m_Env);
      m_pMoleculeManager->set_BathGasMolecule(Bgtxt) ;
    }

    //--------------
    //  The concentration/pressure units are of following formats:
    //  units:
    //   0: particles per cubic centimeter
    //   1: Torr
    //
    //  Allowed input formats are shown below (example units in particles per cubic centimeter).
    //
    //  <me:PTs>
    //    <me:PTset me:units="0">
    //      <me:Prange initial="1e8" increment="2e7" final="2e8">
    //      <me:Trange initial="100" increment="20" final="200">
    //    </me:PTset>
    //  </me:PTs>
    //
    //  The above example creates a matrix of concentration/temperature points of the size:
    //      (number of P points) x (number of T points)
    //
    //  Another example of specifying small numbers of PT points (example units in Torr):
    //
    //  <me:PTs>
    //    <me:PTpair me:units="1" me:P="100" me:T="200">
    //    <me:PTpair me:units="0" me:P="1e18" me:T="298">
    //  </me:PTs>
    //
    //  The looping of the PT points are easy, they are first parsed and all the points are stored in pairs in
    //  vector PandTs, and Mesmer simply loop through all its members.
    //--------------

    PersistPtr ppPTs = ppConditions->XmlMoveTo("me:PTs");
    if(ppPTs)
      readPTs(ppPTs);
    if (!PandTs.size())
      cerr << "No concentration/pressure and temperature specified.";

    // read initial population (needs to be normalized later if their sum not equals to 1.0)
    PersistPtr ppInitalPopulation = ppConditions->XmlMoveTo("me:InitalPopulation");
    if (ppInitalPopulation)
      m_pReactionManager->setInitialPopulation(ppInitalPopulation);

    PersistPtr ppControl = ppIOPtr->XmlMoveTo("me:control");
    if(ppControl)
    {
      m_Env.testDOSEnabled              = ppControl->XmlReadBoolean("me:testDOS");
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
      if (!m_Env.useTheSameCellNumber && m_Env.MaximumTemperature != 0.0){
        m_Env.useTheSameCellNumber = true;
      }
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

    //defaults
    const string DefaultUnit = "PPCC";
    const double DefaultConcentration = 1e17;
    const double DefaultTmeperature = 298.;

    //
    // check for grid values of temperatures and concentrations
    //
    PersistPtr ppPTset = pp->XmlMoveTo("me:PTset");
    while (ppPTset){
      txt = ppPTset->XmlReadValue("me:units");
      string this_units = txt;
      if (!txt){
        cerr << "No units provided. Default units " << DefaultUnit << " are used.";
        this_units = DefaultUnit;
      }

      // if user does not input any value for temperature and concentration, give a Default set of concentration and pressure
      // for simulation
      std::vector<double> Pvals, Tvals;
      if(!ReadRange("me:Prange", Pvals, ppPTset)) Pvals.push_back(DefaultConcentration);
      if(!ReadRange("me:Trange", Tvals, ppPTset))   Tvals.push_back(DefaultTmeperature);

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
      double this_P = DefaultConcentration;
      double this_T = DefaultTmeperature;

      txt = ppPTpair->XmlReadValue("me:P");
      if (txt)
        this_P = atof(txt);
      txt = ppPTpair->XmlReadValue("me:T");
      if (txt)
        this_T = atof(txt);

      CandTpair thisPair(getConvertedP(this_units, this_P, this_T), this_T);
      PandTs.push_back(thisPair);

      ppPTpair = ppPTpair->XmlMoveTo("me:PTpair");
    }
  }

  //
  // Begin calculation.
  //
  void System::calculate()
  {
    TimeCount events; unsigned int timeElapsed =0;

    WriteMetadata();

    // Find the highest temperature
    for (unsigned int i = 0; i < PandTs.size(); ++i){
      m_Env.MaximumTemperature = max(m_Env.MaximumTemperature, PandTs[i].temperature);
    }

    //---------------------------------------------
    // looping over temperatures and concentrations
    unsigned int calPoint = 0;
    for (calPoint = 0; calPoint < PandTs.size(); ++calPoint){
      m_Env.beta = 1.0 / (boltzmann_RCpK * PandTs[calPoint].temperature) ; //temporary statements
      m_Env.conc = PandTs[calPoint].concentration;
      // unit of conc: particles per cubic centimeter
      cerr << "\nGrid " << calPoint << "\n";
      ctest << "Condition: conc = " << m_Env.conc << ", temp = " << PandTs[calPoint].temperature << "\n{\n";
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
      m_pReactionManager->diagCollisionOperator(m_Env) ;

      // Time steps loop
      int timestep = 160;
      m_pReactionManager->timeEvolution(timestep, m_Env);

      m_pReactionManager->BartisWidomPhenomenologicalRates();

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

}//namespace
