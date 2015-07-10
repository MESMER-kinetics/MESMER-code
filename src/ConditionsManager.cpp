#include <functional>
#include "Persistence.h"
#include "System.h"
#include "ConditionsManager.h"

using namespace std;
namespace mesmer
{

ConditionsManager::ConditionsManager(System* pSys) : m_pSys(pSys) {}

bool ConditionsManager::ParseBathGas(PersistPtr ppConditions)
{

  m_ppConditions = ppConditions;
  const char* txt = m_ppConditions->XmlReadValue("me:bathGas", optional);//bathgas may be specified in PTs
  if(txt)
  {
    string Bgtxt(txt);
    if(!m_pSys->getMoleculeManager()->addmol(Bgtxt, "bathGas", m_pSys->getEnv(), m_pSys->m_Flags))
      return false;
    m_pSys->getMoleculeManager()->set_BathGasMolecule(Bgtxt) ;
  }
  return true;
}

bool ConditionsManager::ParseConditions()
{
  //Save excess concs as specified in <Reaction>s
  vector<Reaction*> pReacts = m_pSys->getReactionManager()->getReactionsWithExcessReactant();
  for(vector<Reaction*>::iterator it=pReacts.begin();it!=pReacts.end();++it)
    baseExcessConcs[*it] = (*it)->get_concExcessReactant();

  if(!readPTs()) return false;
  if (!PandTs.size())
    cerr << "No pressure and temperature specified." << endl;

  // read initial isomer populations (need to be normalized later if their sum's not equal to 1.0)
  PersistPtr ppInitialPopulation = m_ppConditions->XmlMoveTo("me:InitialPopulation");
  if (ppInitialPopulation)
    m_pSys->getReactionManager()->setInitialPopulation(ppInitialPopulation);
  return true;
}

bool ConditionsManager::getConditions (vector<double> &Temperature, vector<double> &Concentration)
{
  bool status(true) ;

  for (size_t calPoint = 0; calPoint < PandTs.size(); ++calPoint)
  {
    double temp = PandTs[calPoint].m_temperature ; 
    Temperature.push_back(temp) ;

    m_pSys->getEnv().conc = PandTs[calPoint].m_concentration; // unit of conc: particles per cubic centimeter
    Concentration.push_back(m_pSys->getEnv().conc) ;
  }
  return status ;
}

  // This is a function for reading concentration/pressure and temperature conditions.
  bool ConditionsManager::readPTs()
  {
    PersistPtr pp = m_ppConditions;
    while(pp = pp->XmlMoveTo("me:PTs")) //can have multiple <me:PTs>
    {
      const char* txt;

      //default unit, pressure and temperature //now in defaults.xml
      //const string default_unit = "PPCC";
      //const double default_P = 1e17;
      //const double default_T = 298.;

      // check for grid values of temperatures and concentrations
      PersistPtr ppPTset = pp->XmlMoveTo("me:PTset");
      while (ppPTset)
      {
        string this_units;
        txt = ppPTset->XmlReadValue("me:units", optional);
        if(!txt)
          txt = ppPTset->XmlReadValue("units");
        if (txt)
          this_units = txt;


				Precision this_precision(UNDEFINED_PRECISION);
        txt = ppPTset->XmlReadValue("me:precision", optional);
        if(!txt)
          txt = ppPTset->XmlReadValue("precision");
        if (txt)
          this_precision = txtToPrecision(txt) ;

        std::vector<double> Pvals, Tvals;
        if(!ReadRange("me:Prange", Pvals, ppPTset) || !ReadRange("me:Trange", Tvals, ppPTset))
          return false;

        const char* bathGasName = m_pSys->getMoleculeManager()->get_BathGasName().c_str();
        for (size_t i(0) ; i < Pvals.size(); ++i){
          for (size_t j(0) ; j < Tvals.size(); ++j){
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
      //
      PersistPtr ppPTpair = pp->XmlMoveTo("me:PTpair");
      while (ppPTpair){
        string this_units;
         //use default only if there are no common units specified
        txt = ppPTpair->XmlReadValue("units", !common_units);
        if(!txt)
          txt = ppPTpair->XmlReadValue("me:units", optional);
        if (txt)
          this_units = txt;
        else if(common_units)
        {
          this_units = common_units;
          ppPTpair->XmlWriteAttribute("units",common_units);
        }

        double this_P, this_T;
        this_P = ppPTpair->XmlReadDouble("me:P", optional);
        this_T = ppPTpair->XmlReadDouble("me:T", optional);
        if(IsNan(this_P))
          this_P = ppPTpair->XmlReadDouble("P", true); //preferred forms
        if(IsNan(this_T))
          this_T = ppPTpair->XmlReadDouble("T", true);

        m_pSys->getEnv().MaximumTemperature = max(m_pSys->getEnv().MaximumTemperature, this_T);

        Precision this_precision(UNDEFINED_PRECISION);
        txt = ppPTpair->XmlReadValue("me:precision", optional); //an element
        if(!txt) 
          //use default only if there is no common precision specified
          txt = ppPTpair->XmlReadValue("precision", !common_precision); //an attribute
        if(!txt)
        {
          txt = common_precision;
          ppPTpair->XmlWriteAttribute("precision",common_precision);
        }
        if (txt) {
          this_precision = txtToPrecision(txt) ;
        }

        // Bath gas specific to this PT 
        const char* bathGasName = ppPTpair->XmlReadValue("me:bathGas", optional);
        if(!bathGasName)
          bathGasName = ppPTpair->XmlReadValue("bathGas", optional); //attribute
        if(!bathGasName)
          bathGasName = common_bathgas;
        if(!bathGasName)// if not specified use the general bath gas molecule name
          bathGasName = m_pSys->getMoleculeManager()->get_BathGasName().c_str();
        ppPTpair->XmlWriteAttribute("bathGas",bathGasName);

        // Excess Reactant Concentration for this PT.
        // If there is more than one reaction with an excessReactant specified, 
        // either they all have to be the same molecule, whose concentration is set here,
        // or a refReaction attribute is needed to specify the reaction to which
        // this excessConc is applied.
        // If it is necessary to specify more than one excessReactantConc for a PTPair,
        // this attribute-based method cannot be used and an alternative element-based
        // form (not yet coded) is needed.

        map<Reaction*,double> thisExcessConcs(baseExcessConcs);
        double excessConc = ppPTpair->XmlReadDouble("excessReactantConc", optional);

        // If no excessCReactantConc here, use the one on PTs (if present).
        if(IsNan(excessConc) && !IsNan(common_excessReactantConc))
          excessConc = common_excessReactantConc;

        if(!IsNan(excessConc))
        {
          const char* idtxt = ppPTpair->XmlReadValue("refReactionExcess", optional);
          if(!idtxt && common_reaction_excess) // If no refReactionExcess here use the one on PTs.
            idtxt=common_reaction_excess;
          if(idtxt)
          {
            Reaction* pReact = m_pSys->getReactionManager()->find(idtxt);
            if(!pReact)
              cerr << "Unknown refReactionExcess (for excess reactant concentration)" << endl;
            else
              thisExcessConcs[pReact] = excessConc;
          }
          else
          {
            //check that all excessReactants are the same molecule
            vector<Reaction*> pReacts = m_pSys->getReactionManager()->getReactionsWithExcessReactant();
            vector<Reaction*>::iterator it=pReacts.begin();
            if(pReacts.size()>1)
            {
              Molecule* pMol = (*it)->getExcessReactant();
              assert(pMol);
              for(;it!=pReacts.end();++it)
              {
                if(pMol != (*it)->getExcessReactant())
                {
                  cerr << "The attribute excessReactantConc on PTs or PTPair can be used only "
                       << "if every excess Reactant is the same molecule or if refReactionExcess is specified."
                       << endl;
                  throw std::runtime_error("Erroneous excessReactantConc attribute in PTPair");
                }
              }
            }
            //set all excess reactant concentions to the specified value
            for(it=pReacts.begin(); it!=pReacts.end();++it)
              thisExcessConcs[*it] = excessConc;
          }
        }

        CandTpair thisPair(getConvertedP(this_units, this_P, this_T), this_T,
          this_precision, bathGasName, thisExcessConcs);
        cinfo << this_P << this_units << ", " << this_T << "K at " << txt 
          << " precision" << " with " << bathGasName;
        if(!IsNan(excessConc))
         cinfo << ". Excess Reactant Conc = " << excessConc << " particles per cc";
        cinfo << endl; 

        // Extract experimental rate coefficient values for chiSquare calculation.
        PersistPtr ppExpRate = ppPTpair->XmlMoveTo("me:experimentalRate");
        while (ppExpRate){
          double rateValue(0.0), errorValue(0.0); 
          string refReaction;
          txt = ppExpRate->XmlRead();
          stringstream s1(txt); s1 >> rateValue;
          txt = ppExpRate->XmlReadValue("ref1", optional);
          if(!txt)
          {
            txt = common_ref1;
            ppExpRate->XmlWriteAttribute("ref1",txt);
          }
          string ref1(txt);

          txt = ppExpRate->XmlReadValue("ref2", optional);
          if(!txt)
          {
            txt = common_ref2;
            ppExpRate->XmlWriteAttribute("ref2",txt);
          }
          string ref2(txt);

          txt = ppExpRate->XmlReadValue("refReaction", optional);
          if(!txt)
            txt=common_reaction;
          if (txt) {
            stringstream s3(txt); s3 >> refReaction ;
          }
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalRates(ppExpRate, ref1, ref2, refReaction, rateValue, errorValue);
          ppExpRate = ppExpRate->XmlMoveTo("me:experimentalRate"); //***do we need to loop here?
        }

        // Extract experimental yield values for chiSquare calculation.
        ppExpRate = ppPTpair->XmlMoveTo("me:experimentalYield");
        while (ppExpRate){
          double yield(0.0), errorValue(0.0); 
          txt = ppExpRate->XmlRead();
          stringstream s1(txt); s1 >> yield;
          txt = ppExpRate->XmlReadValue("ref", optional);
          string ref(txt ? txt : common_ref);
          txt = ppExpRate->XmlReadValue("yieldTime", false);
          string yieldTime ;
          if (txt) {
            stringstream s3(txt); s3 >> yieldTime ;
          } else {
            yieldTime = "-1.0" ;
          }
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalYields(ppExpRate, ref, yieldTime, yield, errorValue);
          ppExpRate = ppExpRate->XmlMoveTo("me:experimentalYield"); //***do we need to loop here?
        }

        // Extract experimental eigenvalues for chiSquare calculation.
        ppExpRate = ppPTpair->XmlMoveTo("me:experimentalEigenvalue");
        while (ppExpRate){
          double eigenValue(0.0), errorValue(0.0); 
          txt = ppExpRate->XmlRead();
          stringstream s1(txt); s1 >> eigenValue;
          string EigenvalueID(ppExpRate->XmlReadValue("EigenvalueID")) ;
          stringstream s4(ppExpRate->XmlReadValue("error")); s4 >> errorValue;
          thisPair.set_experimentalEigenvalues(ppExpRate, EigenvalueID, eigenValue, errorValue);
          ppExpRate = ppExpRate->XmlMoveTo("me:experimentalEigenvalue"); //***do we need to loop here?
        }

        //Read in all experimental time-series data for analysis

        double startTime = ppPTpair->XmlReadDouble("startTime", optional); //attribute on PTPair
        string timeUnits = ppPTpair->XmlReadValue("timeUnits");
        bool rawDataOK = true;
        PersistPtr ppRawData = ppPTpair;
        while (ppRawData = ppRawData->XmlMoveTo("me:rawData"))
        {
          RawDataSet ds;
          ds.m_Name = ppRawData->XmlReadValue("name", optional);
          double startTime2 = ppRawData->XmlReadDouble("startTime", optional);//attribute on rawData
          if (!IsNan(startTime2))
            startTime = startTime2;

          ds.m_excessConc = ppRawData->XmlReadDouble("excessReactantConc", optional);
          if (IsNan(ds.m_excessConc))
            ds.m_excessConc = excessConc; //from PTPair
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

            //ds.data.erase(remove_if(ds.data.begin(), ds.data.end(),
            //  [startTime](pair<double,double> pr){return pr.first < startTime;}), ds.data.end()); //C++11

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

  double ConditionsManager::getMaxTemperature()
  {
    // Find the highest temperature
    double Tmax=0;
    for (size_t i(0); i < PandTs.size(); ++i)
      Tmax = max(m_pSys->getEnv().MaximumTemperature, PandTs[i].m_temperature);
    return Tmax;
  }

  void ConditionsManager::getAllBathGases(std::set<std::string>& bathGases)
  {
    for(size_t i(0) ; i!=PandTs.size(); ++i)
      bathGases.insert(PandTs[i].m_pBathGasName);
  }

}//namespace


