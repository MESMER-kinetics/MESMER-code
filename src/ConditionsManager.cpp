#include <numeric>
#include <functional>
#include "Persistence.h"
#include "System.h"
#include "ConditionsManager.h"
#include "TimeCounter.h"

using namespace std;
namespace mesmer
{

  ConditionsManager::ConditionsManager(System* pSys) : m_pSys(pSys),
    m_ppConditions(),
    m_ppAnalysis(),
    baseExcessConcs(),
    PandTs(),
    m_pParallelManager(NULL),
    generalAnalysisData(),
    currentSet(-1) {
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

  bool ConditionsManager::getConditions(vector<double>& Temperature, vector<double>& Concentration)
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

        string group("");
        txt = ppPTset->XmlReadValue("me:group", optional);
        if (!txt)
          txt = ppPTset->XmlReadValue("group", optional);
        if (txt)
          group = txt;

        string ref1;
        txt = ppPTset->XmlReadValue("me:ref1", optional);
        if (!txt)
          txt = ppPTset->XmlReadValue("ref1", optional);
        if (txt)
          ref1 = txt;

        string ref2;
        txt = ppPTset->XmlReadValue("me:ref2", optional);
        if (!txt)
          txt = ppPTset->XmlReadValue("ref2", optional);
        if (txt)
          ref2 = txt;

        string refReaction;
        txt = ppPTset->XmlReadValue("me:refReaction", optional);
        if (!txt)
          txt = ppPTset->XmlReadValue("refReaction", optional);
        if (txt)
          refReaction = txt;

        std::vector<double> Pvals, Tvals;
        if (!ReadRange("me:Prange", Pvals, ppPTset) || !ReadRange("me:Trange", Tvals, ppPTset))
          return false;

        const char* bathGasName = m_pSys->getMoleculeManager()->get_BathGasName().c_str();
        for (size_t j(0); j < Tvals.size(); ++j) {
          for (size_t i(0); i < Pvals.size(); ++i) {
            CandTpair thisPair(getConvertedP(this_units, Pvals[i], Tvals[j]), Tvals[j], this_precision, bathGasName, baseExcessConcs, group.c_str());
            thisPair.set_experimentalRates(ppPTset, ref1, ref2, refReaction, 0.0, 0.0);
            PandTs.push_back(thisPair);
            m_pSys->getEnv().MaximumTemperature = max(m_pSys->getEnv().MaximumTemperature, thisPair.m_temperature);
          }
        }
        ppPTset = ppPTset->XmlMoveTo("me:PTset");
      }

      // These attributes can be on <me:PTs> and apply to its child elements (to shorten them).
      const char* common_precision = pp->XmlReadValue("precision", optional);
      const char* common_units = pp->XmlReadValue("units", optional);
      const char* common_bathgas = pp->XmlReadValue("bathGas", optional);
      const char* common_ref1 = pp->XmlReadValue("ref1", optional);
      const char* common_ref2 = pp->XmlReadValue("ref2", optional);
      const char* common_ref = pp->XmlReadValue("ref", optional);
      const char* common_reaction = pp->XmlReadValue("refReaction", optional);
      const char* common_reaction_excess = pp->XmlReadValue("refReactionExcess", optional);
      const char* common_group = pp->XmlReadValue("group", optional);
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

        // Bath gas concentration

        double bathGasConc = getConvertedP(this_units, this_P, this_T);
        m_pSys->getEnv().conc = bathGasConc;

        // Excess Reactant Concentration for this PT.
        // If there is more than one reaction with an excessReactant specified, it is 
        // necessary to use the element-based form which is tested for first. Otherwise
        // either all reactions have the same excess molecule, whose concentration is 
        // set here or a refReaction attribute is needed to specify the reaction to which
        // this excessConc is applied.

        map<Reaction*, double> thisExcessConcs;
        // Save excess concs as specified in <Reaction>s.
        vector<Reaction*> pReacts = m_pSys->getReactionManager()->getReactionsWithExcessReactant();
        for (vector<Reaction*>::iterator it = pReacts.begin(); it != pReacts.end(); ++it)
          thisExcessConcs[*it] = (*it)->get_concExcessReactant();

        // First, see if there is a specification of all excess species concentrations.
        PersistPtr ppXSConcArray = ppPTpair->XmlMoveTo("me:excessReactantConcArray");
        if (ppXSConcArray) {
          while (ppXSConcArray = ppXSConcArray->XmlMoveTo("me:excessReactantConc")) {

            const char* ptxt = ppXSConcArray->XmlReadValue("reactionRef");
            string ref;
            if (ptxt)
              ref = string(ptxt);
            else {
              stringstream msg;
              msg << "Failed to located a reaction reference" << endl;
              throw(std::runtime_error(msg.str()));
            }

            Reaction* pReaction = m_pSys->getReactionManager()->find(ref);
            if (!pReaction) {
              stringstream msg;
              msg << "Failed to located reaction referred to in excess concentration array" << ref << "." << endl;
              throw(std::runtime_error(msg.str()));
            }

            double conc = ppXSConcArray->XmlReadDouble("concentration");
            bool bConcPercent = ppXSConcArray->XmlReadBoolean("percent");
            thisExcessConcs[pReaction] = (bConcPercent) ? conc * bathGasConc : conc;
          }
        }
        else {

          // A common excess reactant is assume. Try to find a generic excess 
          // species concentration and apply it to all excess species. 
          double excessConc = ppPTpair->XmlReadDouble("excessReactantConc", optional);
          bool bPercentExcessConc = ppPTpair->XmlMoveTo("percentExcessReactantConc");

          // If no excessReactantConc here, use the one on PTs (if present).
          if (IsNan(excessConc) && !IsNan(common_excessReactantConc)) {
            excessConc = common_excessReactantConc;
          }

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
              else {
                // Make sure the excess reactant has the same concentration in all reactions in which it participates.
                Molecule* pMol = pReact->getExcessReactant();
                vector<Reaction*>::iterator it = pReacts.begin();
                for (; it != pReacts.end(); ++it) {
                  if (pMol == (*it)->getExcessReactant())
                    thisExcessConcs[*it] = (bPercentExcessConc) ? excessConc * bathGasConc : excessConc;
                }
              }
            }
            else
            {
              // Check that all excessReactants are the same molecule.
              vector<Reaction*>::iterator it = pReacts.begin();
              if (pReacts.size() > 1)
              {
                Molecule* pMol = (*it)->getExcessReactant();
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
              // Set all excess reactant concentions to the specified value.
              for (it = pReacts.begin(); it != pReacts.end(); ++it)
                thisExcessConcs[*it] = (bPercentExcessConc) ? excessConc * bathGasConc : excessConc;
            }
          }
        }

        // Group to which this PT belongs.
        const char* group = ppPTpair->XmlReadValue("me:group", optional);
        if (!group)
          group = ppPTpair->XmlReadValue("group", optional); //attribute
        if (!group)
          group = common_group;
        if (!group) // If not specified leave as blank.
          group = "default";
        ppPTpair->XmlWriteAttribute("group", group);

        CandTpair thisPair(getConvertedP(this_units, this_P, this_T), this_T,
          this_precision, bathGasName, thisExcessConcs, group);
        cinfo << this_P << this_units << ", " << this_T << "K at " << txt
          << " precision" << " with " << bathGasName << endl;
        for (vector<Reaction*>::iterator it = pReacts.begin(); it != pReacts.end(); ++it)
          cinfo << "Excess Reactant Conc. for reaction " << (*it)->getName() << " = " << thisExcessConcs[*it] << " particles per cc" << endl;

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
            if (txt)
              ppExpRate->XmlWriteAttribute("ref2", txt);
          }
          string ref2 = (txt) ? txt : "";

          txt = ppExpRate->XmlReadValue("refReaction", optional);
          if (!txt)
            txt = common_reaction;
          if (txt) {
            stringstream s3(txt); s3 >> refReaction;
          }
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalRates(ppExpRate, ref1, ref2, refReaction, rateValue, errorValue);
          m_pSys->m_Flags.bIndependentErrors &= (errorValue > 0.0); // Assume independent errors.
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
          m_pSys->m_Flags.bIndependentErrors &= (errorValue > 0.0); // Assume independent errors.
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
          m_pSys->m_Flags.bIndependentErrors &= (errorValue > 0.0); // Assume independent errors.
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
          ds.m_pPersistPtr = ppRawData;
          ds.m_Name = ppRawData->XmlReadValue("name", optional);
          ds.m_pPersistPtr = ppRawData;
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

          // Check to see if diffusive loss is included and print a warning.
          PersistPtr ppDiffLoss = ppRawData->XmlMoveTo("diffusiveLoss");
          if (!ppDiffLoss) {
            ppDiffLoss = ppRawData->XmlMoveTo("me:diffusiveLossArray");
          }
          if (ppDiffLoss) {
            cinfo << "Note diffusive loss terms are only included during fitting procedures." << once << endl;
          }

          ds.m_weight = ppRawData->XmlReadDouble("weight", optional);
          if (IsNan(ds.m_weight)) {
            ds.m_weight = 1.0; // Assign a default weight of unity.
            ppRawData->XmlWriteAttribute("weight", toString(ds.m_weight));
          }

          PersistPtr ppTimes = ppRawData->XmlMoveTo("me:times");
          PersistPtr ppSignals = ppRawData->XmlMoveTo("me:signals");
          if (!(ppTimes && ppSignals))
          {
            cerr << "Missing me:times or me:signals element" << endl;
            return false;
          }
          stringstream times(ppRawData->XmlReadValue("me:times", optional));
          stringstream signals(ppRawData->XmlReadValue("me:signals", optional));
          stringstream timesout, signalsout;
          double t, val;
          size_t traceSize(0);
          while (times.good() && signals.good())
          {
            times >> t;
            signals >> val;
            ds.data.push_back(make_pair(getConvertedTime(timeUnits, t), val));
            timesout << fixed << t << ' ';
            signalsout << fixed << val << ' ';
            traceSize++;
          }
          ds.m_traceSize = traceSize;

          if (times.good() || signals.good())
          {
            cerr << "In the rawData set " << ds.m_Name
              << " the number of times is not equal to the number of signals.";
            rawDataOK = false; //but check other rawData sets first
          }

          //Replace the me:times and me:signal elements if they use scientific format by versions that do not.
          //(so XSLT 1.0 can be used).
          if (times.str().find_first_of("eE") != string::npos && ppRawData->XmlDeleteElement(ppTimes))
            ppRawData->XmlWriteValueElement("me:times", timesout.str());
          if (signals.str().find_first_of("eE") != string::npos && ppRawData->XmlDeleteElement(ppSignals))
            ppRawData->XmlWriteValueElement("me:signals", signalsout.str());

          //If startTime has been specified, remove data before startTime 
          if (!IsNan(startTime))
            ds.data.erase(remove_if(ds.data.begin(), ds.data.end(), Before(startTime)), ds.data.end());

          thisPair.m_rawDataSets.push_back(ds);

          // Raw data usually is based on photon counts or similar, for which there is limited
          // information about errors, so we must infer errors from the data itself, so there 
          // can be no independent model assessment. 
          m_pSys->m_Flags.bIndependentErrors = false;
        }

        if (!rawDataOK) return false;

        PandTs.push_back(thisPair);
        ppPTpair = ppPTpair->XmlMoveTo("me:PTpair");
      }

    }

    // Normalize weigths on raw data (if any).
    NormalizeExptWeights();

    return true;
  }

  // Normalize Experimental weights.
  void ConditionsManager::NormalizeExptWeights()
  {
    double sum(0.0);
    for (size_t i(0); i < PandTs.size(); ++i) {
      vector<RawDataSet>& data = PandTs[i].m_rawDataSets;
      for (size_t j(0); j < data.size(); ++j) {
        data[j].m_weight = 1.0 / data[j].m_weight;
        sum += data[j].m_weight;
      }
    }

    sum /= double(PandTs.size());

    if (sum > 0.0) {
      for (size_t i(0); i < PandTs.size(); ++i) {
        vector<RawDataSet>& data = PandTs[i].m_rawDataSets;
        for (size_t j(0); j < data.size(); ++j) {
          data[j].m_weight /= sum;
        }
      }
    }

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
      else if ((txt = pp->XmlReadValue("factor", false)))//optional attribute
      {
        double fctr = atof(txt);
        txt = pp->XmlReadValue("final"); //if have "factor" must have "final"
        if (!txt)
          return false;
        for (double val = vals.back() * fctr; val <= atof(txt); val *= fctr)
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

  size_t  ConditionsManager::getTotalNumPoints() const {
    size_t npnts(0);
    for (size_t i(0); i < PandTs.size(); ++i) {
      npnts += PandTs[i].m_rates.size();
      npnts += PandTs[i].m_yields.size();
      npnts += PandTs[i].m_eigenvalues.size();
      for (size_t j(0); j < PandTs[i].m_rawDataSets.size(); ++j) {
        npnts += PandTs[i].m_rawDataSets[j].data.size();
      }
    }

    return npnts;
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

  // Write data table. Called only when System::calculate called with writeReport=true.
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
      const vector<RawDataSet>& Traces = PandTs[calPoint].m_rawDataSets;
      for (size_t i(0); i < Traces.size(); ++i) {
        rateCoeffTable << conditions.str() << setw(15) << "Time" << setw(16) << "Expt." << setw(16) << "Calc." << endl;
        const RawDataSet& Trace = Traces[i];
        for (size_t j(0); j < Trace.data.size(); ++j) {
          double time = Trace.data[j].first;
          double expt = Trace.data[j].second;
          double calc = Trace.m_calcTrace[j];
          rateCoeffTable << setw(30) << " " << formatFloat(time, 6, 15) << "," << formatFloat(expt, 6, 15) << "," << formatFloat(calc, 6, 15) << endl;
        }
        rateCoeffTable << endl;
        rateCoeffTable << setw(30) << "Weight: " << Trace.m_weight << endl;
        rateCoeffTable << endl;
      }
    }

    rateCoeffTable << endl;

    cinfo << rateCoeffTable.str();

    stringstream calcsignals;
    for (size_t calPoint(0); calPoint < PandTs.size(); calPoint++) //every P,T point
    {
      const vector<RawDataSet>& Traces = PandTs[calPoint].m_rawDataSets; //Traces points to the rawdatasets at a particular P,T
      for (size_t i(0); i < Traces.size(); ++i) //each rawdataset
      {
        calcsignals.str() = "";
        const RawDataSet& Trace = Traces[i];
        for (size_t j(0); j < Trace.data.size(); ++j)
          calcsignals << fixed << Trace.m_calcTrace[j] << ' ';
        calcsignals << endl;
      }
      //    PersistPtr pp = Traces[calPoint].m_pPersistPtr;[calPoint] is wrong and crashes; but do not need this line
      //      pp->XmlWriteValueElement("me:calcSignals", calcsignals.str());Now written in AddCalcValToXml()
      //      cout << calcsignals.str() << endl; //temporary comment
    }
  }

  // Reconcile table across processes.
  void ConditionsManager::reconcileTable()
  {
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
      vector<RawDataSet>& Trace = PandTs[calPoint].m_rawDataSets;
      for (size_t i(0); i < Trace.size(); ++i) {
        vector<double> tmp = Trace[i].m_calcTrace;
        m_pParallelManager->broadcastVecDouble(tmp, broadcastRank);
        Trace[i].m_calcTrace = tmp;
        double tmp2 = Trace[i].m_weight;
        m_pParallelManager->broadcastDouble(&tmp2, 1, broadcastRank);
        Trace[i].m_weight = tmp2;
      }
    }
  }

  // Calculate ChiSquared.
  void ConditionsManager::calculateChiSquared(double& totalChiSquared, vector<double>& residuals) const {

    // Data can be assigned to groups and if data has not been assigned to a group it
    // is placed in a default group. The intention here is to stop the shear size of
    // a data set from, say, one laboratory, out weighing all other data, especially
    // in the situation that it has a systematic error. See T. Turányi and co-workers
    // Z. Phys. Chem. 2020; 234(7–9): 1329–1357.

    // Map of the location of all data by group.
    map<string, vector<size_t> > groupMap;
    for (size_t i(0); i < PandTs.size(); ++i) {
      string group = PandTs[i].m_group;
      groupMap[group].push_back(i);
    }

    bool bIndependentErrors = m_pSys->m_Flags.bIndependentErrors;
    bool bUseTraceWeighting = m_pSys->m_Flags.useTraceWeighting;
    bool bUpdateTraceWeights = m_pSys->m_Flags.updateTraceWeights;

    totalChiSquared = 0.0;
    residuals.clear();
    residuals.resize(PandTs.size(), 0.0);
    map<string, vector<size_t> >::const_iterator itr = groupMap.begin();
    for (; itr != groupMap.end(); itr++) {

      const vector<size_t> location = itr->second;
      double chiSquared = 0.0;

      for (size_t ii(0); ii < location.size(); ii++) {
        size_t calPoint = location[ii];
        const vector<conditionSet>& rates = PandTs[calPoint].m_rates;
        for (size_t i(0); i < rates.size(); ++i) {
          double error = (bIndependentErrors) ? rates[i].m_error : 1.0;
          double tmp = (rates[i].m_value - rates[i].m_calcValue) / error;
          residuals[calPoint] += tmp;
          chiSquared += tmp * tmp;
        }
        const vector<conditionSet>& yields = PandTs[calPoint].m_yields;
        for (size_t i(0); i < yields.size(); ++i) {
          double error = (bIndependentErrors) ? yields[i].m_error : 1.0;
          double tmp = (yields[i].m_value - yields[i].m_calcValue) / error;
          residuals[calPoint] += tmp;
          chiSquared += tmp * tmp;
        }
        const vector<conditionSet>& eigenvalues = PandTs[calPoint].m_eigenvalues;
        for (size_t i(0); i < eigenvalues.size(); ++i) {
          double error = (bIndependentErrors) ? eigenvalues[i].m_error : 1.0;
          double tmp = (eigenvalues[i].m_value - eigenvalues[i].m_calcValue) / error;
          residuals[calPoint] += tmp;
          chiSquared += tmp * tmp;
        }
        const vector<RawDataSet>& Trace = PandTs[calPoint].m_rawDataSets;
        for (size_t i(0); i < Trace.size(); ++i) {
          const RawDataSet& dataSet = Trace[i];

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

          const vector<double>& signal = dataSet.m_calcTrace;

          // Accumulate Chi^2 for trace.

          double diff(0.0), localChi(0.0);
          for (size_t j(0); j < signal.size(); ++j) {
            diff = (signal[j] - expSignal[j]);
            localChi += (diff * diff);
          }

          // Taking the residual to be the euclidian distance between functions (represented as vectors).
          if (bUseTraceWeighting) {
            chiSquared += localChi * dataSet.m_weight;
            residuals[calPoint] += sqrt(localChi) * dataSet.m_weight;
            if (bUpdateTraceWeights) {
              dataSet.m_weight = localChi / double(signal.size() - 1);;
            }
          }
          else {
            chiSquared += localChi;
            residuals[calPoint] += sqrt(localChi);
            dataSet.m_weight = localChi / double(signal.size() - 1);
          }

        }

      } // End Main conditions loop.

      totalChiSquared += chiSquared / double(location.size());

    }

    totalChiSquared *= double(PandTs.size()) / double(groupMap.size());

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

      // Add elements for trace data.
      for (size_t j(0); j < PandTs[i].m_rawDataSets.size(); j++) {
        const RawDataSet& dataSet = PandTs[i].m_rawDataSets[j];
        const vector<double>& calcTrace = dataSet.m_calcTrace;
        stringstream ss;
        // Add baseline if required.
        for (size_t k(0); k < dataSet.m_traceSize - calcTrace.size(); k++) {
          ss << fixed << 0.0 << " ";
        }
        //  Add calculated trace.
        for (size_t k(0); k < calcTrace.size(); k++) {
          string tmp = formatFloatNWS(calcTrace[k], 6, 15);
          ss << fixed << tmp << " ";
        }
        PersistPtr pp = dataSet.m_pPersistPtr;
        pp->XmlWriteAttribute("weight", toString(dataSet.m_weight));
        PersistPtr ppcs = pp->XmlMoveTo("me:calcSignals");
        if (ppcs) {
          ppcs->XmlWrite(ss.str());
        }
        else
          pp->XmlWriteValueElement("me:calcSignals", ss.str());
      }

    } // End Main conditions loop.
  }

  // Write calculated data to XML.
  void ConditionsManager::WriteDataToXml(PersistPtr pp, const vector<conditionSet>& data) const {

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

      int* itmp = &analysisData.m_number;
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

      vector<double>& grnTimes = analysisData.m_grnTimes;
      m_pParallelManager->broadcastVecDouble(grnTimes, broadcastRank);

      // The following piece of code deals with the transfer of the STL map that holds the energy distributions.

      if (analysisData.m_grnTimes.size() != 0) {
        vector<string> molID;
        map<string, vector<vector<double> > >::const_iterator itr = analysisData.m_grnDists.begin();
        for (; itr != analysisData.m_grnDists.end(); itr++) {
          molID.push_back(itr->first);
        }
        m_pParallelManager->broadcastVecString(molID, broadcastRank);
        vector<vector<double> > grnDists;
        for (size_t i(0); i < molID.size(); i++) {
          grnDists = analysisData.m_grnDists[molID[i]];
          if (grnDists.size() == 0)
            grnDists.resize(analysisData.m_grnTimes.size());
          for (size_t j(0); j < grnDists.size(); j++) {
            vector<double>& dtmp = grnDists[j];
            m_pParallelManager->broadcastVecDouble(dtmp, broadcastRank);
          }
          if (rank != broadcastRank)
            analysisData.m_grnDists[molID[i]] = grnDists;
        }
      }

      stmp = analysisData.m_warning;
      m_pParallelManager->broadcastString(stmp, broadcastRank);
      analysisData.m_warning = stmp;

    }

    // There will usually be an <analysis> section for every calculate().
    // When fitting set m_Flags.overwriteXmlAnalysis true.
    string comment = m_pSys->m_Flags.overwriteXmlAnalysis ?
      "Only selected calculations shown here" : "All calculations shown";
    PersistPtr ppAnalysis = m_ppIOPtr->XmlWriteMainElement("me:analysis", comment, true);
    m_ppAnalysis = ppAnalysis;
    if (Rdouble::withRange().size() != 0)
    {
      PersistPtr ppParams = ppAnalysis->XmlWriteElement("me:parameters");
      for (size_t i = 0; i != Rdouble::withRange().size(); ++i)
      {
        string varname = Rdouble::withRange()[i]->get_varname();
        //The varnames contain ':' or '(' or')' which is incompatible with XML. replace by '-'.
        replace(varname.begin(), varname.end(), ':', '-');
        replace(varname.begin(), varname.end(), '(', '-');
        replace(varname.begin(), varname.end(), ')', '-');
        //Tolerate spaces in varnames
        replace(varname.begin(), varname.end(), ' ', '_');
        stringstream ss;
        ss << *Rdouble::withRange()[i];
        ppParams->XmlWriteAttribute(varname, ss.str());
      }
    }

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

      // Grain distributions, Output density to XML (Chris)

      if (analysisData.m_grnTimes.size() != 0)
      {
        PersistPtr ppGrainList = ppAnalysis->XmlWriteElement("me:grainPopulationList");
        ppGrainList->XmlWriteAttribute("T", toString(PandTs[i].m_temperature));
        ppGrainList->XmlWriteAttribute("conc", toString(PandTs[i].m_concentration));

        map<string, vector<vector<double> > >::const_iterator itr = analysisData.m_grnDists.begin();
        for (; itr != analysisData.m_grnDists.end(); itr++) {

          for (size_t iTime(0); iTime < analysisData.m_grnTimes.size(); iTime++) {

            vector<double> density = itr->second[iTime];

            double totalPop = accumulate(density.begin(), density.end(), 0.0);
            PersistPtr ppGrainPop = ppGrainList->XmlWriteElement("me:grainPopulation");
            {
              ppGrainPop->XmlWriteAttribute("ref", itr->first);
              ppGrainPop->XmlWriteAttribute("time", toString(analysisData.m_grnTimes[iTime]));
              //ppGrainPop->XmlWriteAttribute("logTime", toString(log10(Times[iTime])));
              ppGrainPop->XmlWriteAttribute("me:pop", toString(totalPop));
              ppGrainPop->XmlWriteAttribute("units", "cm-1");

              // Output grain population at each grain energy in two forms:
              // normalised - sum of all = 1; and log of unnormalised value
              stringstream ssgpop;
              for (size_t j(0); j < density.size(); ++j)
              {
                if (density[j] >= 1e-11) //ignore point if density is very small
                {
                  ssgpop.str("");
                  ssgpop << fixed << setprecision(6) << density[j] / totalPop;
                  PersistPtr ppGrain = ppGrainPop->XmlWriteElement("me:grain");
                  ppGrain->XmlWriteAttribute("energy", toString((j + 0.5) * m_pSys->getEnv().GrainSize)); //cm-1
                  ppGrain->XmlWriteAttribute("normpop", ssgpop.str());
                  ppGrain->XmlWriteAttribute("logpop", toString(log10(density[j])));
                }
              }
            }
          }
        }
      }

      // BW Rate coefficients.

      PersistPtr ppList = ppAnalysis->XmlWriteElement("me:rateList");
      ppList->XmlWriteAttribute("T", toString(PandTs[i].m_temperature));
      ppList->XmlWriteAttribute("conc", toString(PandTs[i].m_concentration));
      ppList->XmlWriteAttribute("bathGas", toString(PandTs[i].m_pBathGasName.c_str()));
      ppList->XmlWriteAttribute("units", "s-1");
      if (!analysisData.m_warning.empty()) ppList->XmlWriteValueElement("me:warning", analysisData.m_warning);

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

  }

  // Write the general analysis data to <me:analysis> section of the XML output.

  void ConditionsManager::WriteXMLandClear()
  {
    if (m_ppAnalysis) {
      generalAnalysisData.writeCovariance(m_ppAnalysis);

      generalAnalysisData.clear();
    }
  }

}//namespace

