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
#include "TimeCounter.h"
#include "ParseForPlugin.h"
#include "ConditionsManager.h"
#include "ParallelManager.h"
#include "XMLPersist.h"
#include <cassert>

using namespace std;
using namespace Constants;

namespace mesmer
{
  //Global varaiable and function declared in System.h
  std::string libfile;

  PersistPtr GetFromLibrary(const std::string molName, PersistPtr ppMolList)
  {
    //Search the library of molecules, copy to the main XML file and return a pointer to the copy
    // (or library version if ot copied).
    PersistPtr ppNewMol;
    if (molName.empty())
      return ppNewMol;
    static PersistPtr ppLib; //initiallized by system to NULL
    if (!ppLib)
    {
      ppLib = XMLPersist::XmlLoad(libfile, "");
      if (libfile.empty() || !ppLib)
      {
        cwarn << "Could not find Library file to search it for missing molecule(s)." << endl;
        return NULL;
      }
    }
    PersistPtr ppLibMolList = ppLib->XmlMoveTo("moleculeList");
    if (!ppLibMolList)
      ppLibMolList = ppLib; //Can do without <moleculeList>

    PersistPtr ppMol = ppLibMolList->XmlMoveTo("molecule");
    string tmolName(molName);
    const char* libId = NULL;
    while (ppMol)
    {
      if (tmolName == ppMol->XmlReadValue("id", false))
      {
        //ignore library molecules with attribute active="false"
        const char* active = ppMol->XmlReadValue("active", optional);
        if (!active || strcmp(active, "false"))
        {
          //Check whether this match is an alias, e.g.
          // <molecule id="aliasName" ref="libraryName"/> 
          const char* txt = ppMol->XmlReadValue("ref", optional);
          if (txt)
          {
            libId = txt;
            tmolName = libId; //continue looking for real molecule
          }
          else //not an alias
          {
            if (ppMolList) //no copy if not specified
            {
              //Delete a molecule of same name in datafile, if present
              PersistPtr ppOldMol = ppMolList;
              while (ppOldMol = ppOldMol->XmlMoveTo("molecule"))
              {
                if (molName == ppOldMol->XmlReadValue("id", false))
                  break;
              }

              //Copy a matching molecule from library to the main XML file
              //Replace old version if present
              ppMol = ppMolList->XmlCopy(ppMol, ppOldMol);

              cinfo << molName << " copied from " << libfile;
              //Write its provenance, under metadatList if present
              PersistPtr pp = ppMol->XmlMoveTo("metadataList");
              pp = pp ? pp : ppMol;
              pp->XmlWriteMetadata("copiedFrom", libfile);
              if (libId)//originally an alias 
              {
                ppMol->XmlWriteAttribute("id", molName);
                ppMol->XmlWriteAttribute("libId", libId);
                cinfo << " Original library id = " << libId;
              }
              cinfo << endl;
            }
            return ppMol;
          }
        }
      }
      ppMol = ppMol->XmlMoveTo("molecule");
    }
    cinfo << "Could not find " << molName << " in " << libfile << endl;
    return ppMol; // empty, no suitable molecule found
  }


  System::System(const std::string& libraryfilename, ParallelManager* pParallelManager) :
    m_pMoleculeManager(0),
    m_pReactionManager(0),
    m_pConditionsManager(0),
    m_pParallelManager(pParallelManager),
    m_pTitle(NULL),
    m_pDescription(NULL)
  {
    libfile = libraryfilename;
  }

  System::~System() {
    delete m_pReactionManager;
    delete m_pMoleculeManager;
    delete m_pConditionsManager;
  }

  // Initialize the System object.
  bool System::initialize(void) {
    m_pMoleculeManager = new MoleculeManager();
    m_pReactionManager = new ReactionManager(m_pMoleculeManager);
    m_pConditionsManager = new ConditionsManager(this);
    return m_collisionOperator.initialize(m_pMoleculeManager, m_pReactionManager, m_pConditionsManager);
  }

  //
  // Parse an input data file.
  //
  bool System::parse(PersistPtr ppIOPtr)
  {
    initializeConversionMaps();
    m_ppIOPtr = ppIOPtr;

    m_pTitle = ppIOPtr->XmlReadValue("title", false);
    if (!m_pTitle)
      m_pTitle = ppIOPtr->XmlReadValue("me:title", false);
    m_pDescription = ppIOPtr->XmlReadValue("description", false);
    if (!m_pDescription)
      m_pDescription = ppIOPtr->XmlReadValue("me:description", false);

    //Molecule List
    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    m_pMoleculeManager->set_PersistPtr(ppMolList ? ppMolList : ppIOPtr);

    // What are we going to do?

    PersistPtr ppControl = ppIOPtr;//first time at top, then gets next block
    while (ppControl = ppControl->XmlMoveTo("me:control"))
    {
      //Ignore a <control> that has an attribute  active="false","no" or "0"
      const char* txt = ppControl->XmlReadValue("active", optional);
      if (txt && !ppControl->XmlReadBoolean("active"))
        continue;

      m_CalcMethod = ParseForPlugin<CalcMethod>(this, "me:calcMethod", ppControl);
      if (!m_CalcMethod)
        return false;
      m_CalcMethodsForEachControl.push_back(m_CalcMethod); //for each <control> block

      static PersistPtr ppConditions;
      // Parse molecules, reactions, parameters and conditions only
      // in conjunction with the first control block.
      if (m_CalcMethodsForEachControl.size() == 1)
      {
        if (m_CalcMethod->DoesOwnParsing()) {
          // UnitTests,etc. does its own file parsing.
          m_FlagsForEachControl.push_back(m_Flags);
          continue; //return true ;
        }

        // Model Parameters 
        PersistPtr ppParams = ppIOPtr->XmlMoveTo("me:modelParameters");
        bool readModelParams(true);
        if (!ppParams)
        {
          // Add this section if it does not exist, to contain defaults
          ppParams = ppIOPtr->XmlWriteElement("me:modelParameters", "me:control");

          // No <me:modelParameters> found so give calcMethod the opportunity
          // to provide defaults. If it does it returns true.
          // Most calcMethods will return false.
          readModelParams = !m_CalcMethod->DoesOwnParsing(CalcMethod::MODELPARAMS);
        }

        if (readModelParams)
        {
          m_Env.CellSize = ppParams->XmlReadDouble("me:cellSize", optional);
          if (IsNan(m_Env.CellSize)) {
            m_Env.CellSize = 1.0; // Default cell size in cm-1.
          }
          size_t nCells = ppParams->XmlReadInteger("me:numberOfCells", optional);
          if (nCells) m_Env.MaxCell = nCells;

          // The grain size and grain number are linked via the total maximum energy,
          // so only one of then is independent. Look for Max. number of grains first
          // and if that fails look for the grain size which must be there explictly
          // or available as a default. Note these are provisional estimates that may
          // be revised later after the system is parsed.

          m_Env.MaxGrn = ppParams->XmlReadInteger("me:numberOfGrains", optional);
          if (m_Env.MaxGrn == 0) {
            m_Env.GrainSize = ppParams->XmlReadInteger("me:grainSize");
            if (m_Env.GrainSize == 0) {
              throw (std::runtime_error("Neither grain number or grain size is specified. At least one of these us required"));
            }
            else {
              m_Env.MaxGrn = size_t(double(m_Env.MaxCell) * m_Env.CellSize / double(m_Env.GrainSize));
            }
          }
          else {
            m_Env.GrainSize = size_t(double(m_Env.MaxCell) * m_Env.CellSize / double(m_Env.MaxGrn));
          }

          m_Env.MaximumTemperature = ppParams->XmlReadDouble("me:maxTemperature", optional);
          if (IsNan(m_Env.MaximumTemperature))
            m_Env.MaximumTemperature = 0.0;
          m_Env.EAboveHill = ppParams->XmlReadDouble("me:energyAboveTheTopHill");
          m_Env.useBasisSetMethod = ppParams->XmlReadBoolean("me:runBasisSetMethodroutines");
          if (m_Env.useBasisSetMethod) {
            PersistPtr ppBasisSet = ppParams->XmlMoveTo("me:runBasisSetMethodroutines");
            if (ppBasisSet) {
              m_Env.nBasisSet = ppBasisSet->XmlReadInteger("me:numberBasisFunctions");
            }
            else {
              cerr << "Basis set method requested but number of basis functions unspecified.";
              return false;
            }
          }
          PersistPtr pp = ppParams->XmlMoveTo("me:automaticallySetMaxEne");
          if (pp) {
            m_Flags.autoSetMaxEne = true;
            m_Flags.popThreshold = ppParams->XmlReadDouble("me:automaticallySetMaxEne");
          }
        }

        cinfo.flush();

        //Write the energy convention as an attribute on <moleculeList>
        m_pMoleculeManager->WriteEnergyConvention();

        //Reaction Conditions
        ppConditions = ppIOPtr->XmlMoveTo("me:conditions");
        if (!ppConditions)
        {
          cerr << "No conditions specified" << endl;
          //Except with ThermodynamicTable and UnitTests returns false to abort
          return m_CalcMethod->DoesOwnParsing(CalcMethod::NOCONDITIONSOK);
        }

        m_ConditionsForEachControl.push_back(m_pConditionsManager);
        if (!m_pConditionsManager->ParseBathGas(ppConditions))
          return false;;

        //Reaction List
        PersistPtr ppReacList = ppIOPtr->XmlMoveTo("reactionList");
        if (!ppReacList)
        {
          cerr << "No reactions have been specified";
          return false;
        }
        if (!m_pReactionManager->addreactions(ppReacList, m_Env, m_Flags))
          return false;

        if (!m_pConditionsManager->ParseConditions())
          return false;
      }
      else // second and subsequent <control> blocks
      {
        // If there is another <conditions> block parse it and save a pointer to it.
        PersistPtr pp = ppConditions->XmlMoveTo("me:conditions");
        if (pp)
        {
          m_pConditionsManager = new ConditionsManager(this);
          ppConditions = pp;
          if (!m_pConditionsManager->ParseBathGas(pp))
            return false;
          if (!m_pConditionsManager->ParseConditions())
            return false;
          m_ConditionsForEachControl.push_back(m_pConditionsManager);
        }
        else //no additional <control> block found
          m_ConditionsForEachControl.push_back(NULL);
      }

      if (ppControl)
      {
        m_Flags.testDOSEnabled = ppControl->XmlReadBoolean("me:testDOS");
        m_Flags.testRateConstantEnabled = ppControl->XmlReadBoolean("me:testRateConstants");
        m_Flags.microRateEnabled = ppControl->XmlReadBoolean("me:testMicroRates");
        if (m_Flags.microRateEnabled) {
          PersistPtr pptmr = ppControl->XmlMoveTo("me:testMicroRates");
          double CRCTitl = pptmr->XmlReadDouble("Tstep");
          double CRCTMin = pptmr->XmlReadDouble("Tmin");
          double CRCTMax = pptmr->XmlReadDouble("Tmax");
          Reaction::setTestInterval(CRCTMin, CRCTMax, CRCTitl);

          m_Env.MaximumTemperature = max(m_Env.MaximumTemperature, CRCTMax);
        }
        m_Flags.grainBoltzmannEnabled = ppControl->XmlReadBoolean("me:printGrainBoltzmann");
        m_Flags.grainDOSEnabled = ppControl->XmlReadBoolean("me:printGrainDOS");
        m_Flags.grainTSsosEnabled = ppControl->XmlReadBoolean("me:printTSsos");
        m_Flags.cellDOSEnabled = ppControl->XmlReadBoolean("me:printCellDOS");
        m_Flags.reactionOCSEnabled = ppControl->XmlReadBoolean("me:printReactionOperatorColumnSums");
        m_Flags.kfEGrainsEnabled = ppControl->XmlReadBoolean("me:printGrainkfE");
        m_Flags.kbEGrainsEnabled = ppControl->XmlReadBoolean("me:printGrainkbE");
        // Both Tunnelling and Tunneling will work
        m_Flags.TunnellingCoeffEnabled = ppControl->XmlReadBoolean("me:printTunnellingCoefficients");
        if (!m_Flags.TunnellingCoeffEnabled)
          m_Flags.TunnellingCoeffEnabled = ppControl->XmlReadBoolean("me:printTunnelingCoefficients");
        m_Flags.CrossingCoeffEnabled = ppControl->XmlReadBoolean("me:printCrossingCoefficients");
        m_Flags.cellFluxEnabled = ppControl->XmlReadBoolean("me:printCellTransitionStateFlux");
        m_Flags.grainFluxEnabled = ppControl->XmlReadBoolean("me:printGrainTransitionStateFlux");
        m_Flags.rateCoefficientsOnly = ppControl->XmlReadBoolean("me:calculateRateCoefficientsOnly");
        m_Flags.useTheSameCellNumber = ppControl->XmlReadBoolean("me:useTheSameCellNumberForAllConditions");
        m_Flags.grainedProfileEnabled = ppControl->XmlReadBoolean("me:printGrainedSpeciesProfile");
        m_Flags.speciesProfileEnabled = ppControl->XmlReadBoolean("me:printSpeciesProfile");
        m_Flags.printPhenomenologicalEvolution = ppControl->XmlReadBoolean("me:printPhenomenologicalEvolution");
        m_Flags.printSinkFluxes = ppControl->XmlReadBoolean("me:printSinkFluxes");
        m_Flags.InitialDistEnabled = ppControl->XmlReadBoolean("me:printInitialDistribution");
        m_Flags.viewEvents = ppControl->XmlReadBoolean("me:printEventsTimeStamps");
        m_Flags.allowSmallerDEDown = ppControl->XmlReadBoolean("me:allowSmallerDeltaEDown");
        m_Flags.useDOSweightedDT = ppControl->XmlReadBoolean("me:useDOSweighedDownWardTransition");
        m_Flags.bForceMacroDetailedBalance = ppControl->XmlReadBoolean("me:ForceMacroDetailedBalance");
        m_Flags.useOrigFreqForZPECorrection = ppControl->XmlReadBoolean("me:useOrigFreqForZPECorrection");
        m_Flags.printZMatrixRatios = ppControl->XmlReadBoolean("me:printZMatrixRatios");

        m_Flags.useTraceWeighting = ppControl->XmlReadBoolean("me:useTraceWeighting");
        if (m_Flags.useTraceWeighting) {
          PersistPtr ppTraceWeights = ppControl->XmlMoveTo("me:useTraceWeighting");
          m_Flags.updateTraceWeights = ppTraceWeights->XmlReadBoolean("updateTraceWeights");
        }

        // Check to see if me:ForceMacroDetailedBalance is set in defaults.
        const char* text = ppControl->XmlReadValue("me:ForceMacroDetailedBalance", true);
        if (text) {
          string s(text);
          if (s == "1" || s == "yes" || s == "true") {
            m_Flags.bForceMacroDetailedBalance = true;
          }
        }

        // System configuration information
        if (ppControl->XmlReadBoolean("me:runPlatformDependentPrecisionCheck")) configuration();

        const char* txtPCOP = ppControl->XmlReadValue("me:printCollisionOperatorLevel", false);
        if (txtPCOP) {
          istringstream ss(txtPCOP);
          ss >> m_Flags.showCollisionOperator;
        }

        const char* txt = ppControl->XmlReadValue("me:eigenvalues", optional);
        if (txt && strcmp(txt, "all") == 0)
          m_Flags.printEigenValuesNum = -1;
        else
          //no longer uses defaults.xml
          m_Flags.printEigenValuesNum = ppControl->XmlReadInteger("me:eigenvalues", optional);

        // Should eigenvectors be printed?
        if (m_Flags.printEigenValuesNum != 0) {
          m_Flags.printEigenVectors = ppControl->XmlReadBoolean("me:printEigenVectors");
        }

        const char* txtROS = ppControl->XmlReadValue("me:printReactionOperatorSize", optional);
        if (txtROS) {
          istringstream sROS(txtROS);
          sROS >> m_Flags.printReactionOperatorNum;
        }

        const char* txtSTI = ppControl->XmlReadValue("me:shortestTimeOfInterest", optional);
        if (txtSTI) {
          istringstream ss(txtSTI);
          ss >> m_Flags.shortestTimeOfInterest;
        }

        const char* txtMET = ppControl->XmlReadValue("me:MaximumEvolutionTime", optional);
        if (txtMET) {
          istringstream ss(txtMET);
          ss >> m_Flags.maxEvolutionTime;
        }

        PersistPtr pp = ppControl->XmlMoveTo("me:printGrainProfileAtTime");
        if (pp)
          if (!m_collisionOperator.parseDataForGrainProfileAtTime(pp))
            return false;

        //Copy flags to vector and make a fresh m_Flags
        m_FlagsForEachControl.push_back(m_Flags);
        m_Flags = MesmerFlags();
      }
    }

    if (!Rdouble::SetUpLinkedVars())
      return false;

    return true;
  }

  //
  // Main calculation method.
  //
  void System::executeCalculation()
  {
    if (!m_ConditionsForEachControl.empty())
    {
      size_t nConditionBlocks = m_ConditionsForEachControl.size()
        - count(m_ConditionsForEachControl.begin(), m_ConditionsForEachControl.end(), (ConditionsManager*)NULL);
      for (size_t i(0); i < m_CalcMethodsForEachControl.size(); ++i)
      {
        m_CalcMethod = m_CalcMethodsForEachControl[i];
        assert(m_CalcMethod);
        m_Flags = m_FlagsForEachControl[i];
        cinfo << "\n--Execute calcMethod " << m_CalcMethod->getID();

        if (nConditionBlocks > 1)
        {
          if (m_ConditionsForEachControl[i]) //update if non-NULL
            m_pConditionsManager = m_ConditionsForEachControl[i];

          const char a[4][7] = { "first", "second", "third", "fourth" };
          cinfo << " using " << a[min(i, nConditionBlocks)] << " conditions block ";
        }
        cinfo << endl;
        m_CalcMethod->DoCalculation(this);
        m_pConditionsManager->WriteXMLandClear();
      }
    }
    else {
      m_CalcMethod->DoCalculation(this);
      m_pConditionsManager->WriteXMLandClear();
    }

    //Calls Finish() for each Reaction. Usually does nothing except in AssociationReaction
    m_pReactionManager->finish();
  }

  //
  // Begin calculation.
  // over all PT values, constant parameters
  bool System::calculate(double& chiSquare, vector<double>& residuals, bool writeReport)
  {
    // Controls the print-out of grain/cell DOS in each cycle (This is only for source term)
    if (m_Flags.cellDOSEnabled) m_Flags.cyclePrintCellDOS = true;
    if (m_Flags.grainDOSEnabled) m_Flags.cyclePrintGrainDOS = true;

    TimeCount events; unsigned int timeElapsed = 0;
    bool debugOutput = (meErrorLog.GetOutputLevel() == obDebug);

    //
    // Reset microcanonical rate re-calculation flag as parameters, such
    // as reaction threshold may have been altered between invocations of
    // this method.
    //
    for (size_t i(0); i < m_pReactionManager->size(); ++i) {
      (*m_pReactionManager)[i]->resetCalcFlag();
    }

    // Find the highest temperature
    m_Env.MaximumTemperature = m_pConditionsManager->getMaxTemperature();

    chiSquare = 0.0;
    residuals.clear();
    residuals.resize(m_pConditionsManager->getNumPTPoints(), 0.0);

    // Main loop over temperatures and concentrations.

    int rank = m_pParallelManager->rank();
    int nRanks = m_pParallelManager->size();
    m_pParallelManager->setMutex(__FUNCTION__);
    size_t nPoint = m_pConditionsManager->getNumPTPoints();
    size_t addand = (nPoint % nRanks > 0) ? 1 : 0;
    size_t nLoop = nRanks * (nPoint / nRanks + addand);

    for (size_t calPoint(rank); calPoint < nLoop; calPoint += nRanks) {

      int status(-1);
      stringstream error_msg;
      error_msg << "Rank " << rank << ": no throw." << endl;
      try {

        // This condition catches those situations where the number of ranks is not a factor 
        // of the number of (p,T) points and so some ranks have nothing to do but none-the-less 
        // need to commumicate this. 
        if (calPoint < nPoint) {

          m_pConditionsManager->get_analysisData(calPoint)->clear();
          m_pConditionsManager->set_currentData(int(calPoint));

          m_Env.conc = m_pConditionsManager->PTPointConc(calPoint); // Unit of conc: particles per cubic centimetre.
          m_Env.bathGasName = m_pConditionsManager->PTPointBathGas(calPoint);
          m_Env.collisionTemperature = m_pConditionsManager->PTPointTemp(calPoint);
          m_Env.radiationTemperature = m_pConditionsManager->PTPointRadT(calPoint);
          m_Env.radiationAttenuation = m_pConditionsManager->PTPointRadAttn(calPoint);
          m_Env.photolysisFrequency  = m_pConditionsManager->PTPointPhotoFrq(calPoint);
          m_Env.beta = 1.0 / (boltzmann_RCpK * m_Env.collisionTemperature);

          map<Reaction*, double> excessConcs = m_pConditionsManager->PTPointExcessConcs(calPoint);
          //Set excess Concentrations for all reactions
          for (map<Reaction*, double>::iterator it = excessConcs.begin(); it != excessConcs.end(); ++it)
            it->first->set_concExcessReactant(it->second);

          if (writeReport && debugOutput) { cinfo << "PT Grid " << calPoint << endl; }
          Precision precision = m_pConditionsManager->PTPointPrecision(calPoint);
          m_collisionOperator.setPrecision(precision);

          stest << "PT Grid " << calPoint << " Condition: conc = " << m_Env.conc << ", temp = "
            << m_pConditionsManager->PTPointTemp(calPoint);

          switch (precision) {
          case DOUBLE_DOUBLE:
            stest << ", diagonalization precision: double-double\n{\n";
            break;
          case QUAD_DOUBLE:
            stest << ", diagonalization precision: quad-double\n{\n";
            break;
          default:
            stest << ", diagonalization precision: double\n{\n";
          }

          // Build collison matrix for system.
          if (!m_collisionOperator.BuildReactionOperator(m_Env, m_Flags, writeReport))
            throw (std::runtime_error("Failed building system collison operator."));

          if (writeReport && debugOutput) {
            string thisEvent = "Build Collison Operator";
            string mytime = events.setTimeStamp(thisEvent, timeElapsed);
            cinfo << thisEvent;
            if (timeElapsed > 0)
              cinfo << " -- Time elapsed: " << timeElapsed << " seconds.";
            cinfo << endl;
          }

          if (!m_Flags.rateCoefficientsOnly) {

            if (!m_collisionOperator.calculateEquilibriumFractions(m_Env))
              throw (std::runtime_error("Failed calculating equilibrium fractions."));

            // Diagonalise the collision operator.
            m_collisionOperator.diagReactionOperator(m_Flags, m_Env, m_pConditionsManager->get_analysisData(calPoint));

            if (writeReport && debugOutput) {
              string thisEvent = "Diagonalize the Reaction Operator";
              string mytime = events.setTimeStamp(thisEvent, timeElapsed);
              cinfo << thisEvent;
              if (timeElapsed > 0)
                cinfo << " -- Time elapsed: " << timeElapsed << " seconds.";
              cinfo << endl;
            }

            if (!m_Env.useBasisSetMethod) {

              if (m_Flags.bIsSystemSecondOrder) {
                cinfo << "At present it is not possible to print profiles for systems with a second order term." << once << endl;
              }
              else if (writeReport) {
                // Calculate time-dependent properties.
                m_collisionOperator.timeEvolution(m_Flags, m_pConditionsManager->get_analysisData(calPoint));
                if (m_collisionOperator.hasGrainProfileData()) {
                  m_collisionOperator.printGrainProfileAtTime(m_pConditionsManager->get_analysisData(calPoint));
                }
              }

              // Calculate rate coefficients. 
              qdMatrix mesmerRates(1);
              qdMatrix lossRates(1);
              m_collisionOperator.BartisWidomPhenomenologicalRates(mesmerRates, lossRates, m_Flags, m_pConditionsManager->get_analysisData(calPoint));

              // Write out phenomenological rate coefficients.
              m_collisionOperator.PrintPhenomenologicalRates(mesmerRates, lossRates, m_Flags, m_pConditionsManager->get_analysisData(calPoint));

              // For these conditions calculate the for the experimentally observable data types. 

              calcRateCoefficients(mesmerRates, lossRates, calPoint);

              calcYields(calPoint);

              calcEigenvalues(calPoint);

              calcRawData(calPoint);

              stest << "}\n";

            }
            else {

              qdMatrix mesmerRates(1);
              m_collisionOperator.BartisWidomBasisSetRates(mesmerRates, m_Flags);

            }
          }
        }

      }
      catch (std::runtime_error& e) {
        status = rank;
        error_msg.str("");
        error_msg << "Rank " << rank << ": " << e.what() << endl;
      }

      m_pParallelManager->barrier();
      m_pParallelManager->CheckForThrow(status, error_msg.str());
      m_pParallelManager->writeCtest();

    } // End of temperature and concentration loop. 

    // Reduce Chi^2 values and redistribute

    m_pParallelManager->barrier();
    m_pConditionsManager->reconcileTable();
    m_pConditionsManager->AddCalcValToXml();
    m_pParallelManager->clearMutex();


    double ChiSquaredTest(0.0);
    m_pConditionsManager->calculateChiSquared(ChiSquaredTest, residuals);
    chiSquare = ChiSquaredTest;

    if (writeReport) {
      m_pConditionsManager->WriteDataTable();
      m_pConditionsManager->WriteAnalysisXML(m_ppIOPtr);
    }

    stringstream line;
    string thisEvent = "Finish Calculation of P-T case";
    string mytime = events.setTimeStamp(thisEvent, timeElapsed);
    line << thisEvent;
    if (timeElapsed > 0)
      line << " -- Time elapsed: " << timeElapsed << " seconds.";
    line << endl;
    if (writeReport)
      line << m_pConditionsManager->getNumPTPoints() << " temperature/concentration-pressure points calculated." << endl;

    if (m_Flags.viewEvents)
      line << events;
    if (writeReport && debugOutput)
      cinfo << line.str();

    return true;
  }

  //
  // Calculate rate coefficients etc. for a specific condition.
  //
  bool System::calculate(size_t nConds, vector<double>& Quantities, bool writeReport)
  {
    //
    // Reset microcanonical rate re-calculation flag as parameters, such
    // as reaction threshold may have been altered between invocations of
    // this method.
    //
    for (size_t i(0); i < m_pReactionManager->size(); ++i) {
      (*m_pReactionManager)[i]->resetCalcFlag();
    }

    m_Flags.printEigenValuesNum = 0;

    m_Env.collisionTemperature = m_pConditionsManager->PTPointTemp(nConds);
    m_Env.beta = 1.0 / (boltzmann_RCpK * m_Env.collisionTemperature);

    // unit of conc: particles per cubic centimetre
    m_Env.conc = m_pConditionsManager->PTPointConc(nConds);

    // Set excess Concentrations for all reactions.
    map<Reaction*, double> excessConcs = m_pConditionsManager->PTPointExcessConcs(nConds);
    for (map<Reaction*, double>::iterator it = excessConcs.begin(); it != excessConcs.end(); ++it)
      it->first->set_concExcessReactant(it->second);

    m_Env.bathGasName = m_pConditionsManager->PTPointBathGas(nConds);

    Precision precision = m_pConditionsManager->PTPointPrecision(nConds);
    m_collisionOperator.setPrecision(precision);

    // Build collison matrix for system.
    if (!m_collisionOperator.BuildReactionOperator(m_Env, m_Flags))
      throw (std::runtime_error("Failed building system collison operator."));

    if (!m_collisionOperator.calculateEquilibriumFractions(m_Env))
      throw (std::runtime_error("Failed calculating equilibrium fractions."));

    // Diagonalise the collision operator.
    m_collisionOperator.diagReactionOperator(m_Flags, m_Env);

    // Calculate rate coefficients.
    qdMatrix  mesmerRates(1);
    qdMatrix lossRates(1);
    m_collisionOperator.BartisWidomPhenomenologicalRates(mesmerRates, lossRates, m_Flags);

    // For this conditions calculate the experimentally observable data types. 

    stringstream rateCoeffTable;
    calcRateCoefficients(mesmerRates, lossRates, nConds);

    calcYields(nConds);

    calcEigenvalues(nConds);

    double data(0.0);
    while (rateCoeffTable >> data) {
      Quantities.push_back(data);
    }

    return true;
  }

  //
  // Calculate rate coefficients for conditions other than those 
  // supplied directly by user e.g. for analytical representation.
  //
  bool System::calculate(double& Temperature, double& Concentration, Precision precision,
    map<string, double>& phenRates, double& MaxT, const string& bathGas)
  {
    assert(!bathGas.empty());
    m_Env.bathGasName = bathGas;
    m_Flags.printEigenValuesNum = 0;

    m_Env.collisionTemperature = Temperature;
    m_Env.beta = 1.0 / (boltzmann_RCpK * m_Env.collisionTemperature);
    m_Env.MaximumTemperature = MaxT;
    m_Env.conc = Concentration; // Unit of conc: particles per cubic centimeter.

    m_collisionOperator.setPrecision(precision);

    // Build collison matrix for system.
    if (!m_collisionOperator.BuildReactionOperator(m_Env, m_Flags))
      throw (std::runtime_error("Failed building system collison operator."));

    if (!m_collisionOperator.calculateEquilibriumFractions(m_Env))
      throw (std::runtime_error("Failed calculating equilibrium fractions."));

    // Diagonalise the collision operator.
    m_collisionOperator.diagReactionOperator(m_Flags, m_Env);

    // Calculate rate coefficients.
    qdMatrix BWrates(1);
    qdMatrix lossRates(1);
    m_collisionOperator.BartisWidomPhenomenologicalRates(BWrates, lossRates, m_Flags);
    // Write out phenomenological rate coefficients.
    m_collisionOperator.PrintPhenomenologicalRates(BWrates, lossRates, m_Flags, NULL);
    m_collisionOperator.get_phenRates(phenRates);

    return true;
  }

  // Wrapper for single calculation.
  bool System::calculate(size_t nCond, map<string, double>& phenRates) {

    //
    // Reset microcanonical rate re-calculation flag as parameters, such
    // as reaction threshold may have been altered between invocations of
    // this method.
    //
    for (size_t i(0); i < m_pReactionManager->size(); ++i) {
      (*m_pReactionManager)[i]->resetCalcFlag();
    }

    double temp = m_pConditionsManager->PTPointTemp(nCond);
    double conc = m_pConditionsManager->PTPointConc(nCond);
    Precision precision = m_pConditionsManager->PTPointPrecision(nCond);
    string bathGas(m_pConditionsManager->PTPointBathGas(nCond));
    double MaxT(2.0 * temp);
    return calculate(temp, conc, precision, phenRates, MaxT, bathGas);
  };

  void System::calcRateCoefficients(const qdMatrix& mesmerRates, const qdMatrix& lossRates, const size_t calPoint) {

    vector<conditionSet> expRates;
    m_pConditionsManager->get_experimentalRates(calPoint, expRates);
    vector<double> calcRates(expRates.size(), 0.0);
    for (size_t i(0); i < expRates.size(); ++i) {

      string ref1, ref2, refReaction;
      double expRate(0.0), expErr(0.0);
      expRates[i].get_conditionSet(ref1, ref2, refReaction, expRate, expErr);

      // Get the position of ref1 and ref2 inside m_SpeciesSequence vector
      int seqMatrixLoc1(-1), seqMatrixLoc2(-1);
      seqMatrixLoc1 = m_collisionOperator.getSpeciesSequenceIndex(ref1);
      seqMatrixLoc2 = m_collisionOperator.getSpeciesSequenceIndex(ref2);

      if (seqMatrixLoc1 < 0) {
        throw(std::runtime_error("Failed to locate reactant species " + ref1 + " in rate coefficient matrix."));
      }
      else if (seqMatrixLoc2 < 0) {
        //
        // Search for a dissociation reaction.
        // 
        int sinkLocation = m_collisionOperator.getSinkSequenceIndex(refReaction);

        if (sinkLocation < 0) {
          throw(std::runtime_error("Failed to locate sink species in rate coefficient matrix."));
        }

        double rateCoeff = fabs(to_double(lossRates[sinkLocation][seqMatrixLoc1]));

        calcRates[i] = rateCoeff;
      }
      else {

        // 
        // In the following it is assumed that experimental rate coefficients will always 
        // be quoted as absolute values. Since the diagonal values of the BW matrix are
        // negative, their absolute value is required for comparision with experimental
        // values hence the fabs invocation.
        //
        double rateCoeff = fabs(to_double(mesmerRates[seqMatrixLoc2][seqMatrixLoc1]));

        // Is a bimolecular rate coefficient required?

        Reaction* reaction = m_pReactionManager->find(refReaction);
        if (reaction)
          reaction->normalizeRateCoefficient(rateCoeff, ref1);

        calcRates[i] = rateCoeff;
      }
    }

    m_pConditionsManager->set_calculatedRates(calPoint, calcRates);

    return;
  }

  void System::calcYields(const size_t calPoint) {

    vector<conditionSet> expYields;
    m_pConditionsManager->get_experimentalYields(calPoint, expYields);
    if (expYields.size() == 0)
      return;

    vector<double> calcYields(expYields.size(), 0.0);
    //
    // Calculate yields for these conditions. Assume all yields are measured at the same time.
    //
    YieldMap yieldMap;
    try {
      string product, yieldTime, ref;
      double expYield(0.0), expErr(0.0);
      expYields[0].get_conditionSet(product, yieldTime, ref, expYield, expErr);

      double time;
      stringstream s1(yieldTime); s1 >> time;

      m_collisionOperator.calculateYields(yieldMap, time);
    }
    catch (std::runtime_error& e) {
      cerr << "Error: during calculation of Chi^2 for yields:" << endl;
      cerr << e.what() << endl;
      return;
    }

    for (size_t i(0); i < expYields.size(); ++i) {

      string ref1, ref2, refReaction;
      double expYield(0.0), expErr(0.0);
      expYields[i].get_conditionSet(ref1, ref2, refReaction, expYield, expErr);

      // Is there a wild card in ref1?
      size_t found = ref1.find("*");
      bool wild(found != string::npos);
      string base;
      if (wild)
        base = ref1.substr(0, found - 1);

      double yield(0.0);
      YieldMap::const_iterator yielditr = yieldMap.begin();
      for (; yielditr != yieldMap.end(); ++yielditr) {

        Reaction* sinkReaction = yielditr->first;
        vector<Molecule*> pdts;
        sinkReaction->get_products(pdts);
        for (size_t i(0); i < pdts.size(); i++) {
          string pdtname = pdts[i]->getName();
          if (ref1 == pdtname || (wild && pdtname.find(base) != string::npos))
            yield += yielditr->second;
        }

      }

      calcYields[i] = yield;
    }

    m_pConditionsManager->set_calculatedYields(calPoint, calcYields);

    return;
  }

  void System::calcEigenvalues(const size_t calPoint) {

    vector<conditionSet> expEigenvalues;
    m_pConditionsManager->get_experimentalEigenvalues(calPoint, expEigenvalues);
    vector<double> calcEigenvalues(expEigenvalues.size(), 0.0);
    for (size_t i(0); i < expEigenvalues.size(); ++i) {

      string ref1, ref2, strIdEigenvalue;
      double expEigenvalue(0.0), expErr(0.0);
      expEigenvalues[i].get_conditionSet(ref1, ref2, strIdEigenvalue, expEigenvalue, expErr);

      size_t idEigenvalue;
      stringstream s1(strIdEigenvalue); s1 >> idEigenvalue;

      double eigenvalue(0.0);
      try {
        eigenvalue = m_collisionOperator.getEigenvalue(idEigenvalue);
      }
      catch (std::runtime_error& e) {
        cerr << e.what() << endl;
        continue;
      }

      calcEigenvalues[i] = eigenvalue;
    }

    m_pConditionsManager->set_calculatedEigenvalues(calPoint, calcEigenvalues);

    return;
  }

  void System::calcRawData(const size_t calPoint) {

    //
    // With trace data, an indpendent assessment of the precision of the measurements is not possible. 
    // This stems from the fact the signal is based on a number of counts (which are approximately
    // poisson distributed) and is a relative measure. Because of this an independent set of error
    // estimates cannot be obtained, so the Chi^2 statistic is not applicable and the relevant flag
    // needs to be set.
    //

    vector<RawDataSet>& rawDataSets = m_pConditionsManager->get_experimentalrawDataSets(calPoint);
    if (rawDataSets.size() == 0)
      return;

    //
    // Calculate profile for each experimental set.
    //
    for (size_t i(0); i < rawDataSets.size(); ++i) {

      RawDataSet& dataSet = rawDataSets[i];

      // Extract the times for which the trace should be calculated. 
      vector<double> times, expSignal;
      const size_t ntimes = dataSet.data.size();
      for (size_t j(0); j < ntimes; ++j) {
        double time = dataSet.data[j].first;
        if (time > 0.0) {
          times.push_back(time);
          expSignal.push_back(dataSet.data[j].second);
        }
      }

      vector<double> signal(times.size(), 0.0);
      m_collisionOperator.calculateTrace(dataSet.m_ref1, times, signal);

      // Alter the amplitude of the calculated signal by a linear least squares shift.

      double sumef(0.0), sumf2(0.0);
      for (size_t j(0); j < signal.size(); ++j) {
        sumef += expSignal[j] * signal[j];
        sumf2 += signal[j] * signal[j];
      }

      double alpha = sumef / sumf2;
      for (size_t j(0); j < signal.size(); ++j) {
        signal[j] *= alpha;
      }
      dataSet.m_calcTrace = signal;
    }

    return;
  }

  void System::WriteMetadata(const string& infilename)
  {
    PersistPtr ppList = m_ppIOPtr->XmlWriteMainElement("metadataList", "");
    ppList->XmlWriteAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/");
    if (m_pTitle)
      ppList->XmlWriteValueElement("dc:title", m_pTitle);
    if (m_pDescription)
      ppList->XmlWriteValueElement("dc:description", m_pDescription);
    ppList->XmlWriteValueElement("dc:source", infilename);

    ppList->XmlWriteValueElement("dc:creator", "Mesmer v" + string(MESMER_VERSION));
    TimeCount events;
    string timeString = events.setTimeStamp("");
    cinfo << "Write metadata " << timeString << endl;
    ppList->XmlWriteValueElement("dc:date", timeString);

    const char* author = getenv("MESMER_AUTHOR");
    if (!author)
      author = getenv("USERNAME"); //Windows
    if (!author)
      author = getenv("LOGNAME"); //Unix?
    if (!author)
      author = "unknown";
    ppList->XmlWriteValueElement("dc:contributor", author);
  }

  void System::getAllBathGases(set<string>& bathGases)
  {
    for (size_t i(0); i != m_ConditionsForEachControl.size(); ++i)
      if (m_ConditionsForEachControl[i])
        m_ConditionsForEachControl[i]->getAllBathGases(bathGases);

    bathGases.insert(getMoleculeManager()->get_BathGasName());
  }


  void System::configuration(void) {
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

  /*
  System::parse() loops to read all <control> blocks.
  CalcMethod is pushed to a vector<CalcMethod> calcMethods
  Make a MesmerFlags, populate it, push to vector<MesmerFlags> flags

  System::executeCalculation()
  for loop over size of calcMethods
  m_CalcMethod set to current member
  copy from flags vector to m_Flags
  */
} //namespace

