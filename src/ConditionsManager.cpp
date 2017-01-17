#include <functional>
#include "Persistence.h"
#include "System.h"
#include "ConditionsManager.h"

using namespace std;
namespace mesmer
{

  ConditionsManager::ConditionsManager(System* pSys) : m_pSys(pSys) {
    m_pParallelManager = m_pSys->getParallelManager();
  }

  bool ConditionsManager::ParseBathGas(PersistPtr ppConditions)
  {

    m_ppConditions = ppConditions;
    const char* txt = m_ppConditions->XmlReadValue("me:bathGas", optional);//bathgas may be specified in PTs
    if (txt)
    {
      string Bgtxt(txt);
      if (!m_pSys->getMoleculeManager()->addmol(Bgtxt, "bathGas", m_pSys->getEnv(), m_pSys->m_Flags))
        return false;
      m_pSys->getMoleculeManager()->set_BathGasMolecule(Bgtxt);
    }
    return true;
  }

  bool ConditionsManager::ParseConditions()
  {
    //Save excess concs as specified in <Reaction>s
    vector<Reaction*> pReacts = m_pSys->getReactionManager()->getReactionsWithExcessReactant();
    for (vector<Reaction*>::iterator it = pReacts.begin(); it != pReacts.end(); ++it)
      baseExcessConcs[*it] = (*it)->get_concExcessReactant();

    if (!readPTs()) return false;
    if (!PandTs.size())
      cerr << "No pressure and temperature specified." << endl;

    // read initial isomer populations (need to be normalized later if their sum's not equal to 1.0)
    PersistPtr ppInitialPopulation = m_ppConditions->XmlMoveTo("me:InitialPopulation");
    if (ppInitialPopulation)
      m_pSys->getReactionManager()->setInitialPopulation(ppInitialPopulation);
    return true;
  }

  bool ConditionsManager::getConditions(vector<double> &Temperature, vector<double> &Concentration)
  {
    bool status(true);

    for (size_t calPoint = 0; calPoint < PandTs.size(); ++calPoint)
    {
      double temp = PandTs[calPoint].m_temperature;
      Temperature.push_back(temp);

      m_pSys->getEnv().conc = PandTs[calPoint].m_concentration; // unit of conc: particles per cubic centimeter
      Concentration.push_back(m_pSys->getEnv().conc);
    }
    return status;
  }

  // This is a function for reading concentration/pressure and temperature conditions.
  bool ConditionsManager::readPTs()
  {
    PersistPtr pp = m_ppConditions;
    while (pp = pp->XmlMoveTo("me:PTs")) //can have multiple <me:PTs>
    {
      const char* txt;

      //default unit, pressure and temperature are in defaults.xml

      // check for grid values of temperatures and concentrations
      PersistPtr ppPTset = pp->XmlMoveTo("me:PTset");
      while (ppPTset)
      {
        string this_units;
        txt = ppPTset->XmlReadValue("me:units", optional);
        if (!txt)
          txt = ppPTset->XmlReadValue("units");
        if (txt)
          this_units = txt;

        Precision this_precision(UNDEFINED_PRECISION);
        txt = ppPTset->XmlReadValue("me:precision", optional);
        if (!txt)
          txt = ppPTset->XmlReadValue("precision");
        if (txt)
          this_precision = txtToPrecision(txt);

        std::vector<double> Pvals, Tvals;
        if (!ReadRange("me:Prange", Pvals, ppPTset) || !ReadRange("me:Trange", Tvals, ppPTset))
          return false;

        const char* bathGasName = m_pSys->getMoleculeManager()->get_BathGasName().c_str();
        for (size_t i(0); i < Pvals.size(); ++i) {
          for (size_t j(0); j < Tvals.size(); ++j) {
            CandTpair thisPair(getConvertedP(this_units, Pvals[i], Tvals[j]), Tvals[j],
              this_precision, bathGasName, baseExcessConcs);
            PandTs.push_back(thisPair);
            m_pSys->getEnv().MaximumTemperature = max(m_pSys->getEnv().MaximumTemperature, thisPair.m_temperature);
          }
        }
        ppPTset = ppPTset->XmlMoveTo("me:PTset");
      }

      //These attributes can be on <me:PTs> and apply to its child elements (to shorten them).
      const char* common_precision = pp->XmlReadValue("precision", optional);
      const char* common_units = pp->XmlReadValue("units", optional);
      const char* common_bathgas = pp->XmlReadValue("bathGas", optional);
      const char* common_ref1 = pp->XmlReadValue("ref1", optional);
      const char* common_ref2 = pp->XmlReadValue("ref2", optional);
      const char* common_ref = pp->XmlReadValue("ref", optional);
      const char* common_reaction = pp->XmlReadValue("refReaction", optional);
      const char* common_reaction_excess = pp->XmlReadValue("refReactionExcess", optional);
      double common_excessReactantConc = pp->XmlReadDouble("excessReactantConc", optional);

      // Check for individually specified concentration/temperature points.
      PersistPtr ppPTpair = pp->XmlMoveTo("me:PTpair");
      while (ppPTpair)
      {
        // Use default only if there are no common units specified.
        txt = ppPTpair->XmlReadValue("me:units", optional);//deprecated
        if (!txt)
          txt = ppPTpair->XmlReadValue("units", !common_units);
        if (!txt)
          txt = common_units;
        string this_units(txt);

        double this_P, this_T;
        this_P = ppPTpair->XmlReadDouble("me:P", optional);
        this_T = ppPTpair->XmlReadDouble("me:T", optional);
        if (IsNan(this_P))
          this_P = ppPTpair->XmlReadDouble("P", true); //preferred forms
        if (IsNan(this_T))
          this_T = ppPTpair->XmlReadDouble("T", true);

        m_pSys->getEnv().MaximumTemperature = max(m_pSys->getEnv().MaximumTemperature, this_T);

        Precision this_precision(UNDEFINED_PRECISION);
        txt = ppPTpair->XmlReadValue("me:precision", optional); //an element
        if (!txt)
          //use default only if there is no common precision specified
          txt = ppPTpair->XmlReadValue("precision", !common_precision); //an attribute
        if (!txt)
        {
          txt = common_precision;
          ppPTpair->XmlWriteAttribute("precision", common_precision);
        }
        if (txt) {
          this_precision = txtToPrecision(txt);
        }

        // Bath gas specific to this PT 
        const char* bathGasName = ppPTpair->XmlReadValue("me:bathGas", optional);
        if (!bathGasName)
          bathGasName = ppPTpair->XmlReadValue("bathGas", optional); //attribute
        if (!bathGasName)
          bathGasName = common_bathgas;
        if (!bathGasName)// if not specified use the general bath gas molecule name
          bathGasName = m_pSys->getMoleculeManager()->get_BathGasName().c_str();
        ppPTpair->XmlWriteAttribute("bathGas", bathGasName);

        // Excess Reactant Concentration for this PT.
        // If there is more than one reaction with an excessReactant specified, 
        // either they all have to be the same molecule, whose concentration is set here,
        // or a refReaction attribute is needed to specify the reaction to which
        // this excessConc is applied.
        // If it is necessary to specify more than one excessReactantConc for a PTPair,
        // this attribute-based method cannot be used and an alternative element-based
        // form (not yet coded) is needed.

        map<Reaction*, double> thisExcessConcs(baseExcessConcs);
        double excessConc = ppPTpair->XmlReadDouble("excessReactantConc", optional);

        // If no excessCReactantConc here, use the one on PTs (if present).
        if (IsNan(excessConc) && !IsNan(common_excessReactantConc))
          excessConc = common_excessReactantConc;

        if (!IsNan(excessConc))
        {
          const char* idtxt = ppPTpair->XmlReadValue("refReactionExcess", optional);
          if (!idtxt && common_reaction_excess) // If no refReactionExcess here use the one on PTs.
            idtxt = common_reaction_excess;
          if (idtxt)
          {
            Reaction* pReact = m_pSys->getReactionManager()->find(idtxt);
            if (!pReact)
              cerr << "Unknown refReactionExcess (for excess reactant concentration)" << endl;
            else
              thisExcessConcs[pReact] = excessConc;
          }
          else
          {
            //check that all excessReactants are the same molecule
            vector<Reaction*> pReacts = m_pSys->getReactionManager()->getReactionsWithExcessReactant();
            vector<Reaction*>::iterator it = pReacts.begin();
            if (pReacts.size() > 1)
            {
              Molecule* pMol = (*it)->getExcessReactant();
              assert(pMol);
              for (; it != pReacts.end(); ++it)
              {
                if (pMol != (*it)->getExcessReactant())
                {
                  cerr << "The attribute excessReactantConc on PTs or PTPair can be used only "
                    << "if every excess Reactant is the same molecule or if refReactionExcess is specified."
                    << endl;
                  throw std::runtime_error("Erroneous excessReactantConc attribute in PTPair");
                }
              }
            }
            //set all excess reactant concentions to the specified value
            for (it = pReacts.begin(); it != pReacts.end(); ++it)
              thisExcessConcs[*it] = excessConc;
          }
        }

        CandTpair thisPair(getConvertedP(this_units, this_P, this_T), this_T,
          this_precision, bathGasName, thisExcessConcs);
        cinfo << this_P << this_units << ", " << this_T << "K at " << txt
          << " precision" << " with " << bathGasName;
        if (!IsNan(excessConc))
          cinfo << ". Excess Reactant Conc = " << excessConc << " particles per cc";
        cinfo << endl;

        // Extract experimental rate coefficient values for chiSquare calculation.

        PersistPtr ppExpRate = ppPTpair->XmlMoveTo("me:experimentalRate");
        while (ppExpRate) {
          double rateValue(0.0), errorValue(0.0);
          string refReaction;
          txt = ppExpRate->XmlRead();
          stringstream s1(txt); s1 >> rateValue;
          txt = ppExpRate->XmlReadValue("ref1", optional);
          if (!txt)
          {
            txt = common_ref1;
            if (!txt)
            {
              cerr << "An <experimental Rate> has a missing ref1 attribute,"
                << " which alternatively can be on its parent <PTs> element." << endl;
              return false;
            }
            ppExpRate->XmlWriteAttribute("ref1", txt);
          }
          string ref1(txt);

          txt = ppExpRate->XmlReadValue("ref2", optional);
          if (!txt)
          {
            txt = common_ref2;
            ppExpRate->XmlWriteAttribute("ref2", txt);
          }
          string ref2(txt);

          txt = ppExpRate->XmlReadValue("refReaction", optional);
          if (!txt)
            txt = common_reaction;
          if (txt) {
            stringstream s3(txt); s3 >> refReaction;
          }
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalRates(ppExpRate, ref1, ref2, refReaction, rateValue, errorValue);
          ppExpRate = ppExpRate->XmlMoveTo("me:experimentalRate");
        }

        // Extract experimental yield values for chiSquare calculation.

        ppExpRate = ppPTpair->XmlMoveTo("me:experimentalYield");
        while (ppExpRate) {
          double yield(0.0), errorValue(0.0);
          txt = ppExpRate->XmlRead();
          stringstream s1(txt); s1 >> yield;
          txt = ppExpRate->XmlReadValue("ref", optional);
          string ref(txt ? txt : common_ref);
          txt = ppExpRate->XmlReadValue("yieldTime", false);
          string yieldTime;
          if (txt) {
            stringstream s3(txt); s3 >> yieldTime;
          }
          else {
            yieldTime = "-1.0";
          }
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalYields(ppExpRate, ref, yieldTime, yield, errorValue);
          ppExpRate = ppExpRate->XmlMoveTo("me:experimentalYield");
        }

        // Extract experimental eigenvalues for chiSquare calculation.

        ppExpRate = ppPTpair->XmlMoveTo("me:experimentalEigenvalue");
        while (ppExpRate) {
          double eigenValue(0.0), errorValue(0.0);
          txt = ppExpRate->XmlRead();
          stringstream s1(txt); s1 >> eigenValue;
          string EigenvalueID(ppExpRate->XmlReadValue("EigenvalueID"));
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalEigenvalues(ppExpRate, EigenvalueID, eigenValue, errorValue);
          ppExpRate = ppExpRate->XmlMoveTo("me:experimentalEigenvalue");
        }

        // Read in all experimental time-series data for analysis.

        bool rawDataOK = true;
        PersistPtr ppRawData = ppPTpair;
        while (ppRawData = ppRawData->XmlMoveTo("me:rawData"))
        {
          // These vars read more than once: startTime may have been overwritten;
          // avoids irrelevant default log message when timeUnits is outside loop.
          double startTime = ppPTpair->XmlReadDouble("startTime", optional); // attribute on PTPair
          string timeUnits = ppPTpair->XmlReadValue("timeUnits");

          RawDataSet ds;
          ds.m_Name = ppRawData->XmlReadValue("name", optional);
          try {
            ds.m_ref1 = ppRawData->XmlReadValue("ref");
          }
          catch (...) {
            cerr << "Raw data set is missing identity of monitored species." << endl;
            continue;
          };

          double startTime2 = ppRawData->XmlReadDouble("startTime", optional);// attribute on rawData
          if (!IsNan(startTime2))
            startTime = startTime2;

          ds.m_excessConc = ppRawData->XmlReadDouble("excessReactantConc", optional);
          if (IsNan(ds.m_excessConc))
            ds.m_excessConc = excessConc; // from PTPair
          if (IsNan(ds.m_excessConc))
          {
            cerr << "Missing excessReactantConc on rawData (and not on PTPair)" << endl;
            return false;
          }

          if (!(ppRawData->XmlMoveTo("me:times") && ppRawData->XmlMoveTo("me:signals")))
          {
            cerr << "Missing me:times or me:signals element" << endl;
            return false;
          }
          stringstream times(ppRawData->XmlReadValue("me:times", optional));
          stringstream signals(ppRawData->XmlReadValue("me:signals", optional));
          double t, val;
          while (times.good() && signals.good())
          {
            times >> t;
            signals >> val;
            ds.data.push_back(make_pair(getConvertedTime(timeUnits, t), val));
          }

          if (times.good() || signals.good())
          {
            cerr << "In the rawData set " << ds.m_Name
              << " the number of times is not equal to the number of signals.";
            rawDataOK = false; //but check other rawData sets first
          }

          //If startTime has been specified, remove data before startTime 
          if (!IsNan(startTime))
            ds.data.erase(remove_if(ds.data.begin(), ds.data.end(), Before(startTime)), ds.data.end());

          thisPair.m_rawDataSets.push_back(ds);
        }

        if (!rawDataOK) return false;

        PandTs.push_back(thisPair);
        ppPTpair = ppPTpair->XmlMoveTo("me:PTpair");
      }

    }
    return true;
  }

  bool ConditionsManager::ReadRange(const string& name, vector<double>& vals, PersistPtr ppbase, bool MustBeThere)
  {
    PersistPtr pp = ppbase;
    for (;;)
    {
      const char* txt;
      pp = pp->XmlMoveTo(name);
      if (pp)
        txt = pp->XmlRead(); //element may have a value
      else //no more elements
        break;
      if (!txt)
        txt = pp->XmlReadValue("initial"); //or use value of "initial" attribute
      if (!txt)
        return false;
      vals.push_back(atof(txt));

      if ((txt = pp->XmlReadValue("increment", false)))//optional attribute
      {
        double incr = atof(txt);
        txt = pp->XmlReadValue("final"); //if have "increment" must have "final"
        if (!txt)
          return false;
        for (double val = vals.back() + incr; val <= atof(txt); val += incr)
          vals.push_back(val);
      }
    }
    if (MustBeThere && vals.size() == 0)
    {
      cerr << "Must specify at least one value of " << name;
      return false;
    }
    return true;
  }

  double ConditionsManager::getMaxTemperature()
  {
    // Find the highest temperature
    double Tmax = 0;
    for (size_t i(0); i < PandTs.size(); ++i)
      Tmax = max(m_pSys->getEnv().MaximumTemperature, PandTs[i].m_temperature);
    return Tmax;
  }

  void ConditionsManager::getAllBathGases(std::set<std::string>& bathGases)
  {
    for (size_t i(0); i != PandTs.size(); ++i)
      bathGases.insert(PandTs[i].m_pBathGasName);
  }

  // Write data table.
  void ConditionsManager::WriteDataTable() const
  {

    int rank = m_pParallelManager->rank();
    if (rank > 0)
      return;

    stringstream rateCoeffTable;

    rateCoeffTable << endl;
    rateCoeffTable << "    Temperature  Concentration    Exp. Coeff.    Cal. Coeff." << endl;
    rateCoeffTable << endl;

    for (size_t calPoint(0); calPoint < PandTs.size(); calPoint++) {
      stringstream conditions;
      conditions << formatFloat(PandTs[calPoint].m_temperature, 6, 15);
      conditions << formatFloat(PandTs[calPoint].m_concentration, 6, 15);
      const vector<conditionSet>& rates = PandTs[calPoint].m_rates;
      for (size_t i(0); i < rates.size(); ++i) {
        rateCoeffTable << conditions.str() << formatFloat(rates[i].m_value, 6, 15) << formatFloat(rates[i].m_calcValue, 6, 15) << endl;
      }
      const vector<conditionSet>& yields = PandTs[calPoint].m_yields;
      for (size_t i(0); i < yields.size(); ++i) {
        rateCoeffTable << conditions.str() << formatFloat(yields[i].m_value, 6, 15) << formatFloat(yields[i].m_calcValue, 6, 15) << endl;
      }
      const vector<conditionSet>& eigenvalues = PandTs[calPoint].m_eigenvalues;
      for (size_t i(0); i < eigenvalues.size(); ++i) {
        rateCoeffTable << conditions.str() << formatFloat(eigenvalues[i].m_value, 6, 15) << formatFloat(eigenvalues[i].m_calcValue, 6, 15) << endl;
      }
    }

    rateCoeffTable << endl;

    cinfo << rateCoeffTable.str();
  }

  // Reconcile table across processs.
  void ConditionsManager::reconcileTable()
  {
    int rank = m_pParallelManager->rank();
    int nRanks = m_pParallelManager->size();

    for (size_t calPoint(0); calPoint < PandTs.size(); calPoint++) {

      int broadcastRank = int(calPoint) % nRanks;

      vector<conditionSet>& rates = PandTs[calPoint].m_rates;
      for (size_t i(0); i < rates.size(); ++i) {
        double tmp = rates[i].m_calcValue;
        m_pParallelManager->broadcastDouble(&tmp, 1, broadcastRank);
        rates[i].m_calcValue = tmp;
      }
      vector<conditionSet>& yields = PandTs[calPoint].m_yields;
      for (size_t i(0); i < yields.size(); ++i) {
        double tmp = yields[i].m_calcValue;
        m_pParallelManager->broadcastDouble(&tmp, 1, broadcastRank);
        yields[i].m_calcValue = tmp;
      }
      vector<conditionSet>& eigenvalues = PandTs[calPoint].m_eigenvalues;
      for (size_t i(0); i < eigenvalues.size(); ++i) {
        double tmp = eigenvalues[i].m_calcValue;
        m_pParallelManager->broadcastDouble(&tmp, 1, broadcastRank);
        eigenvalues[i].m_calcValue = tmp;
      }
    }
  }

  // Write calculated date to output.
  void ConditionsManager::AddCalcValToXml() const
  {

    for (size_t i(0); i < PandTs.size(); i++) {

      for (size_t j(0); j < PandTs[i].m_expDataPtrs.size(); j++) {

        // Add extra attribute(s) containing calculated value and time-stamp
        // to <me:experimentalRate> (or similar element).
        PersistPtr pp = PandTs[i].m_expDataPtrs[j];
        TimeCount events;
        string timeString;
        pp->XmlWriteAttribute("calculated", events.setTimeStamp(timeString));

        // Add elements for rate coefficients.
        WriteDataToXml(pp, PandTs[i].m_rates);

        // Add elements for yields.
        WriteDataToXml(pp, PandTs[i].m_yields);

        // Add elements for eigenvalues.
        WriteDataToXml(pp, PandTs[i].m_eigenvalues);
      }
    }
  }

  // Write calculated date to XML.
  void ConditionsManager::WriteDataToXml(PersistPtr pp, const vector<conditionSet> &data) const {

    for (size_t k(0); k < data.size(); k++) {
      stringstream ss;
      ss << data[k].m_calcValue;
      pp->XmlWriteAttribute("calcVal", ss.str());
      pp->XmlReadDouble("calcVal", false);
    }
  }

  // Write the Analysis section of the XML output.
  void ConditionsManager::WriteAnalysisXML(PersistPtr m_ppIOPtr)
  {
    int rank = m_pParallelManager->rank();
    int nRanks = m_pParallelManager->size();

    // Reconcile data.

    for (size_t calPoint(0); calPoint < PandTs.size(); calPoint++) {

      AnalysisData& analysisData = PandTs[calPoint].m_analysisData;

      int broadcastRank = int(calPoint) % nRanks;

      // Eigenvalue data.

      int *itmp = &analysisData.m_number;
      m_pParallelManager->broadcastInteger(itmp, 1, broadcastRank);

      string stmp = analysisData.m_selection;
      m_pParallelManager->broadcastString(stmp, broadcastRank);
      analysisData.m_selection = stmp;

      vector<double>& eigenvalues = analysisData.m_eigenvalues;
      m_pParallelManager->broadcastVecDouble(eigenvalues, broadcastRank);

      // Loss rate coefficients.

      vector<string>& lossRef = analysisData.m_lossRef;
      m_pParallelManager->broadcastVecString(lossRef, broadcastRank);

      vector<double>& lossRateCoeff = analysisData.m_lossRateCoeff;
      m_pParallelManager->broadcastVecDouble(lossRateCoeff, broadcastRank);

      // First order rate coefficients.

      vector<string>& firstOrderReactionType = analysisData.m_firstOrderReactionType;
      m_pParallelManager->broadcastVecString(firstOrderReactionType, broadcastRank);

      vector<string>& firstOrderFromRef = analysisData.m_firstOrderFromRef;
      m_pParallelManager->broadcastVecString(firstOrderFromRef, broadcastRank);

      vector<string>& firstOrderToRef = analysisData.m_firstOrderToRef;
      m_pParallelManager->broadcastVecString(firstOrderToRef, broadcastRank);

      vector<double>& firstOrderRateCoeff = analysisData.m_firstOrderRateCoeff;
      m_pParallelManager->broadcastVecDouble(firstOrderRateCoeff, broadcastRank);

      vector<double>& timePoints = analysisData.m_timePoints;
      m_pParallelManager->broadcastVecDouble(timePoints, broadcastRank);

      vector<string>& aveEnergyRef = analysisData.m_aveEnergyRef;
      m_pParallelManager->broadcastVecString(aveEnergyRef, broadcastRank);

      vector<double>& aveEnergy = analysisData.m_aveEnergy;
      m_pParallelManager->broadcastVecDouble(aveEnergy, broadcastRank);

      vector<string>& PopRef = analysisData.m_PopRef;
      m_pParallelManager->broadcastVecString(PopRef, broadcastRank);

      vector<double>& Pop = analysisData.m_Pop;
      m_pParallelManager->broadcastVecDouble(Pop, broadcastRank);

    }

    string comment = "All calculations shown";
    PersistPtr ppAnalysis = m_ppIOPtr->XmlWriteMainElement("me:analysis", comment, true);

    // Write <analysis> section.

    for (size_t i(0); i < PandTs.size(); i++) {

      AnalysisData& analysisData = PandTs[i].m_analysisData;

      // Eigenvalues.
      PersistPtr ppEigenList = ppAnalysis->XmlWriteElement("me:eigenvalueList");
      ppEigenList->XmlWriteAttribute("number", toString(analysisData.m_number));
      ppEigenList->XmlWriteAttribute("selection", analysisData.m_selection);
      const vector<double>& eigenvalues = analysisData.m_eigenvalues;
      for (size_t j(0); j < eigenvalues.size(); ++j) {
        ppEigenList->XmlWriteValueElement("me:eigenvalue", eigenvalues[j], 6);
      }

      // Populations

      PersistPtr ppPopList;
      if (m_pSys->m_Flags.speciesProfileEnabled)
      {
        ppPopList = ppAnalysis->XmlWriteElement("me:populationList");
        ppPopList->XmlWriteAttribute("T", toString(PandTs[i].m_temperature));
        ppPopList->XmlWriteAttribute("conc", toString(PandTs[i].m_concentration));
        for (size_t timestep(0), k(0); timestep < analysisData.m_timePoints.size(); ++timestep) {
          PersistPtr ppPop = ppPopList->XmlWriteElement("me:population");
          ppPop->XmlWriteAttribute("time", toString(analysisData.m_timePoints[timestep]));
          ppPop->XmlWriteAttribute("logTime", toString(log10(analysisData.m_timePoints[timestep])));
          for (size_t j(0); j < analysisData.m_PopRef.size(); ++j, ++k) {
            PersistPtr ppVal = ppPop->XmlWriteValueElement("me:pop", toString(analysisData.m_Pop[k]));
            ppVal->XmlWriteAttribute("ref", analysisData.m_PopRef[j]);
          }
        }
      }

      // Average energies

      PersistPtr ppAvEList;
      if (m_pSys->m_Flags.grainedProfileEnabled)
      {
        ppAvEList = ppAnalysis->XmlWriteElement("me:avEnergyList");
        ppAvEList->XmlWriteAttribute("T", toString(PandTs[i].m_temperature));
        ppAvEList->XmlWriteAttribute("conc", toString(PandTs[i].m_concentration));
        ppAvEList->XmlWriteAttribute("energyUnits", "kJ/mol");
        for (size_t j(0), k(0); j < analysisData.m_aveEnergyRef.size(); ++j) {
          PersistPtr ppAvEnergy = ppAvEList->XmlWriteElement("me:avEnergy");
          ppAvEnergy->XmlWriteAttribute("ref", analysisData.m_aveEnergyRef[j]);
          for (size_t timestep(0); timestep < analysisData.m_timePoints.size(); ++timestep, ++k) {
            PersistPtr ppAv = ppAvEnergy->XmlWriteValueElement("me:Av", toString(analysisData.m_aveEnergy[k]));
            ppAv->XmlWriteAttribute("time", toString(analysisData.m_timePoints[timestep]));
            ppAv->XmlWriteAttribute("logTime", toString(log10(analysisData.m_timePoints[timestep])));
          }
        }
      }

      // BW Rate coefficients.

      PersistPtr ppList = ppAnalysis->XmlWriteElement("me:rateList");
      ppList->XmlWriteAttribute("T", toString(PandTs[i].m_temperature));
      ppList->XmlWriteAttribute("conc", toString(PandTs[i].m_concentration));
      ppList->XmlWriteAttribute("bathGas", toString(PandTs[i].m_pBathGasName.c_str()));
      ppList->XmlWriteAttribute("units", "s-1");

      for (size_t j(0); j < analysisData.m_lossRateCoeff.size(); j++) {
        PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderLoss", analysisData.m_lossRateCoeff[j]);
        ppItem->XmlWriteAttribute("ref", analysisData.m_lossRef[j]);
      }

      for (size_t j(0); j < analysisData.m_firstOrderRateCoeff.size(); j++) {
        PersistPtr ppItem = ppList->XmlWriteValueElement("me:firstOrderRate", analysisData.m_firstOrderRateCoeff[j]);
        ppItem->XmlWriteAttribute("fromRef", analysisData.m_firstOrderFromRef[j]);
        ppItem->XmlWriteAttribute("toRef", analysisData.m_firstOrderToRef[j]);
        ppItem->XmlWriteAttribute("reactionType", analysisData.m_firstOrderReactionType[j]);
      }

    }

    // Write general data.
    PersistPtr ppCovariance = ppAnalysis->XmlWriteElement("me:covariance", "");
    generalAnalysisData.m_covariance.WriteToXML(ppCovariance);

  }


}//namespace

