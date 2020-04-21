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

#include <functional>
#include "../System.h"
#include "../calcmethod.h"
#include "../gDensityOfStates.h"
#include "../gStructure.h"
#include "ThermodynamicUtils.h"

namespace mesmer
{
  class ThermodynamicTable : public CalcMethod, private ThermodynamicUtils
  {
  public:

    ThermodynamicTable(const char* id) : ThermodynamicUtils(),
      m_id(id),
      m_nTemp(20),
      m_TempInterval(50.0),
      m_Tmin(50),
      m_Tmax(3000),
      m_Tmid(1000),
      m_Unit("kJ/mol"),
      m_unitFctr(1.0),
      m_makeNasaPoly(true),
      m_outputCellVersion(false)
    {
      Register();
    }

    virtual ~ThermodynamicTable() {}
    virtual const char* getID() { return m_id; }
    virtual ThermodynamicTable* Clone() { return new ThermodynamicTable(*this); }

    //Does not do own parsing (returns false with default call),
    //returns true when asked if NOCONDITIONSOK
    //and provides defaults for GrainSize, etc. with MODELPARAMS.
    virtual bool DoesOwnParsing(parseQuery q);

    // Function to do the work
    virtual bool DoCalculation(System* pSys);

  private:

    virtual bool ParseData(PersistPtr pp);

    const char* m_id;

    int m_nTemp;
    double m_TempInterval, m_Tmin, m_Tmax, m_Tmid;
    string m_Unit;
    double m_unitFctr;
    bool m_makeNasaPoly;
    bool m_outputCellVersion;
  };

  ////////////////////////////////////////////////
  //Global instance
  ThermodynamicTable theThermodynamicTable("ThermodynamicTable");
  ///////////////////////////////////////////////

  bool ThermodynamicTable::ParseData(PersistPtr pp)
  {
    ErrorContext c("ThermodynamicTable");
    //Read in parameters in child elements of CalcMethod or use defaults.xml
    //called from ParseForPlugin in System::Parse()
    m_Unit = pp->XmlReadValue("units");
    m_Tmin = pp->XmlReadDouble("me:Tmin");
    m_Tmax = pp->XmlReadDouble("me:Tmax");
    m_TempInterval = pp->XmlReadDouble("me:Tstep");
    // If Tmid is non-zero the NASA polynomial has two temperature ranges.
    m_Tmid = pp->XmlReadDouble("me:Tmid");
    m_unitFctr = 1.0 / kJPerMol_in_RC;
    if (m_Unit == "kcal/mol")
      m_unitFctr = 1.0 / kCalPerMol_in_RC;
    m_outputCellVersion = pp->XmlReadBoolean("me:withCellDOSCalc");

    System* pSys = getParent();
    MoleculeManager* pMoleculeManager = pSys->getMoleculeManager();
    if (pSys->getReactionManager()->size() == 0)
    {
      //No reactions specified, so read in all the molecules here.
      PersistPtr ppmol = pMoleculeManager->get_PersistPtr();
      while (ppmol = ppmol->XmlMoveTo("molecule"))
      {
        // Get the name of the molecule.
        const char* reftxt = ppmol->XmlReadValue("id");
        if (reftxt) {
          string molname(reftxt);
          // Check to see if species is a transition state. If not
          // use molType="forThermo" to activate DOS properties.
          const char* tstxt = ppmol->XmlReadValue("role", optional);
          string role = (tstxt && string(tstxt) == "transitionState") ? string(tstxt) : string("forThermo");
          try {
            pMoleculeManager->addmol(molname, role, pSys->getEnv(), pSys->m_Flags);
          }
          catch (std::runtime_error & e) {
            cerr << e.what() << endl;
            cerr << "Will attempt to treat " << molname << " as undeclared transition state." << endl;
            try {
              pMoleculeManager->addmol(molname, string("transitionState"), pSys->getEnv(), pSys->m_Flags);
            }
            catch (...) {
              cerr << "Definition of " << molname << " required review." << endl;
            }
          }
        }
      }

      // Determine if DOS test information is to appear.
      PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");
      pSys->m_Flags.testDOSEnabled = ppControl->XmlReadBoolean("me:testDOS");
      if (pSys->m_Flags.testDOSEnabled)
        pSys->getEnv().beta = 1.0 / (boltzmann_RCpK * double(m_nTemp) * m_TempInterval);
    }
    if (MolecularComponent::getEnergyConvention() == "arbitrary")
      m_makeNasaPoly = false;

    if (pMoleculeManager->size() == 0 || MolecularComponent::getEnergyConvention().empty())
    {
      cerr << "No suitable molecules have been specified." << endl;
      return false;
    }
    return true;
  }

  // Called with ALL (the default) returns false.
  // Called with MODELPARAMS       sets parameters here and returns true.
  // Called with anything else     returns true.
  bool ThermodynamicTable::DoesOwnParsing(parseQuery q)
  {
    if (q == CalcMethod::MODELPARAMS)
    {
      //Use own default model parameters
      ErrorContext c("ThermodynamicTable. Own default");
      System* pSys = getParent();

      PersistPtr ppParams = pSys->getPersistPtr()->XmlMoveTo("me:modelParameters");
      assert(ppParams);

      MesmerEnv& Env = pSys->getEnv();
      Env.CellSize = 1.0;
      cinfo << "Cell size " << Env.CellSize << " cm-1" << endl;
      PersistPtr pp = ppParams->XmlWriteValueElement("me:cellSize", Env.CellSize);
      pp->XmlWriteAttribute("units", "cm-1");
      pp->XmlWriteAttribute("default", "true");

      Env.MaxCell = 100000;
      cinfo << "Number of cells " << Env.MaxCell << endl;
      pp = ppParams->XmlWriteValueElement("me:numberOfCells", Env.MaxCell);
      pp->XmlWriteAttribute("default", "true");

      //Get the grain parameters from defaults.xml (The reads will always fail.)
      Env.GrainSize = ppParams->XmlReadInteger("me:grainSize");
      //Env.EAboveHill = ppModelParams->XmlReadDouble("me:energyAboveTheTopHill");

      return true;
    }
    return q != ALL;
  }


  bool ThermodynamicTable::DoCalculation(System* pSys)
  {
    // Make provision for the special case of T = 298.15.
    bool   tempLessThan298;
    double temp289(298.15);

    // Loop over all molecules producing a table for each molecule that
    // has an energy specified.
    // A NASA polynomial will also be produced if 7 or more temperatures
    // have been requested in both the upper and lower polynomials 
    // or, if Tmid=0, the single polynomial.

    int nTemps = static_cast<int>(std::floor((m_Tmax - m_Tmin) / m_TempInterval)) + 1;
    int nTempsLower = static_cast<int>(std::floor((m_Tmid - m_Tmin) / m_TempInterval) + 1);
    bool enoughPoints = (!m_Tmid && (nTemps > 6)) || (m_Tmid && (nTempsLower > 6) && ((nTemps - nTempsLower) >= 6));
    if (!enoughPoints)
    {
      cinfo << "Too few data points to fit NASA polynomials." << endl;
      m_makeNasaPoly = false;
    }
    double R = boltzmann_C * AvogadroC;
    if (m_Unit == "kcal/mol") R *= Calorie_in_Joule;

    MoleculeManager* pMoleculeManager = pSys->getMoleculeManager();
    MoleculeManager::constMolIter molItr = pMoleculeManager->begin();
    MoleculeManager::constMolIter molItrEnd = pMoleculeManager->end();
    for (; molItr != molItrEnd; molItr++)
    {
      vector<double> temperature, Hf /* enthalpy of formation / R */;
      Molecule* pmol = molItr->second;

      // Restrict output for molecules without a specified energy.
      double Hf298local = NaN;
      if (m_makeNasaPoly)
      {
        Hf298local = pmol->getDOS().get_Hf298Thermo();
        if (IsNan(Hf298local))
          cinfo << "Restricted thermo output for " << pmol->getName()
          << " because it has no non-arbitrary energy data." << endl;
      }
      PersistPtr pp = pmol->get_PersistentPointer();
      pp = pp->XmlWriteMainElement("me:thermoTable", "", true); //will replace an existing element
      //pp = pp->XmlWriteElement("me:thermoTable");
      pp->XmlWriteAttribute("unitsT", "K");
      pp->XmlWriteAttribute("unitsH", m_Unit);
      pp->XmlWriteAttribute("unitsS", m_Unit.substr(1, m_Unit.length()) + "/K");
      pp->XmlWriteAttribute("unitsG", m_Unit);
      pp->XmlWriteAttribute("unitsCp", m_Unit.substr(1, m_Unit.length()) + "/K");
      if (!IsNan(Hf298local))
        pp->XmlWriteAttribute("unitsHf", m_Unit);

      double S298; // Always calculated. NOTE kJ/mol/K.
      double enthalpy298;
      thermoDynFns thermos;
      pmol->getDOS().thermodynamicsFunctions(298.15, m_unitFctr, thermos);
      enthalpy298 = thermos.enthalpy;
      S298 = thermos.entropy;
      tempLessThan298 = true;
      for (double temp = m_Tmin; temp <= m_Tmax; temp += m_TempInterval)
      {
        double T = temp;
        if (tempLessThan298 && temp > temp289)
        {
          // Special case of T = 289.15.
          tempLessThan298 = false;
          T = temp289;
          temp -= m_TempInterval;
        }
        temperature.push_back(T);

        pmol->getDOS().thermodynamicsFunctions(T, m_unitFctr, thermos);

        PersistPtr ppVal = pp->XmlWriteElement("me:thermoValue");
        ppVal->XmlWriteAttribute("T", T, 2, true);
        ppVal->XmlWriteAttribute("H", thermos.enthalpy, 4, true);
        ppVal->XmlWriteAttribute("S", thermos.entropy * 1000, 4, true);
        ppVal->XmlWriteAttribute("G", thermos.gibbsFreeEnergy, 4, true);
        ppVal->XmlWriteAttribute("Cp", thermos.heatCapacity * 1000, 4, true);
        if (m_outputCellVersion)
        {
          ppVal->XmlWriteAttribute("cellS", thermos.cellEntropy * 1000, 4, true);
          ppVal->XmlWriteAttribute("cellH", thermos.cellEnthalpy, 4, true);
          ppVal->XmlWriteAttribute("cellG", thermos.cellGibbsFreeEnergy, 4, true);
          ppVal->XmlWriteAttribute("cellCp", thermos.cellHeatCapacity * 1000, 4, true);
        }
        if (!IsNan(Hf298local))
        {
          Hf.push_back((thermos.enthalpy - enthalpy298 + Hf298local) * 1000 / R); //e.g. J/mol
          ppVal->XmlWriteAttribute("Hf", Hf.back() * R / 1000, 4, true); //back to kJ/mol
        }
      }

      if (m_makeNasaPoly && !IsNan(Hf298local))
      {
        //Fit NASA polynomial to enthalpy data
        //H/R =a6 + T*a1 + T^2*a2/2 + a3*T^3/3 + a4*T^4/4 + a5*T^5/5
        vector<double> fits1, fits2; //a1, a2/2, a3/3, etc
        if (m_Tmid == 0) // single range (duplicated)
        {
          fits1 = FitPoly(6, temperature.begin(), temperature.end(), Hf.begin());
          fits1[2] *= 2; fits1[3] *= 3; fits1[4] *= 4; fits1[5] *= 5;
          fits2 = fits1;
        }
        else //two ranges
        {
          vector<double>::iterator itermid = find(temperature.begin(), temperature.end(), m_Tmid);
          if (itermid == temperature.end())
          {
            cerr << "In NASA polynomial fits the middle temperature"
              "must be one of the specified temperatures." << endl;
            return false;
          }
          int nlowerrange = itermid - temperature.begin();
          fits1 = FitPoly(6, itermid, temperature.end(), Hf.begin() + nlowerrange); //upper range
          fits2 = FitPoly(6, temperature.begin(), itermid + 1, Hf.begin()); //lower range
          fits1[2] *= 2; fits1[3] *= 3; fits1[4] *= 4; fits1[5] *= 5;
          fits2[2] *= 2; fits2[3] *= 3; fits2[4] *= 4; fits2[5] *= 5;
        }
        vector<double> coeffs(15);
        copy(fits1.begin() + 1, fits1.end(), coeffs.begin());
        copy(fits2.begin() + 1, fits2.end(), coeffs.begin() + 7);

        coeffs[5] = fits1[0];
        coeffs[12] = fits2[0];
        coeffs[14] = Hf298local / R;

        //Set a14 to match S at 298.15K
        coeffs[13] = 0.0;
        coeffs[13] = S298 * 1000 / R - SdivR(coeffs.begin() + 7, 298.15);

        //Set a7 to match a) S at 298K for one range; b) S at Tmid for two range;
        if (m_Tmid == 0)
          coeffs[6] = coeffs[13];
        else
        {
          pmol->getDOS().thermodynamicsFunctions(m_Tmid, m_unitFctr, thermos);
          coeffs[6] = 0.0;
          coeffs[6] = thermos.entropy * 1000 / R - SdivR(coeffs.begin(), m_Tmid);
        }

        // Output to XML using a CML property for Nasa Polynomials
        // previously used in OpenBabel.
        PersistPtr pp = pmol->get_PersistentPointer();
        PersistPtr ppProp = pp->XmlMoveTo("propertyList");
        //ppProp = (ppProp ? ppProp : pp)->XmlWriteMainElement("property","",true);
        ppProp = (ppProp ? ppProp : pp)->XmlWriteElement("property");
        ppProp->XmlWriteAttribute("dictRef", "NasaPolynomial");

        stringstream ss;
        ss << temperature[0];
        PersistPtr ppScalar = ppProp->XmlWriteValueElement("scalar", ss.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaLowT");

        ss.str("");
        ss << temperature.back();
        ppScalar = ppProp->XmlWriteValueElement("scalar", ss.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaHighT");

        ss.str("");
        ss << m_Tmid ? m_Tmid : temperature.back();
        ppScalar = ppProp->XmlWriteValueElement("scalar", ss.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaMidT");

        ppScalar = ppProp->XmlWriteValueElement("scalar", "G");
        ppScalar->XmlWriteAttribute("dictRef", "Phase");

        stringstream vals;
        std::copy(coeffs.begin(), coeffs.end(), ostream_iterator<double>(vals, " "));
        ppScalar = ppProp->XmlWriteValueElement("array", vals.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaCoeffs");
        ppScalar->XmlWriteAttribute("size", "15");

        string poly = WriteChemKinNASAPoly(pmol, coeffs, temperature[0],
          m_Tmid ? m_Tmid : temperature.back(), temperature.back());
        ppScalar = ppProp->XmlWriteValueElement(
          "scalar", poly, true); //Output polynomial as CDATA
        ppScalar->XmlWriteAttribute("dictRef", "NasaPolynomial");
        ctest << poly << endl; // for QA test
      }
    }

    ////TEST FitPoly
    //vector<double> xdata = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    //vector<double> ydata = { 1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321 };
    ////Solution is 3 x2 + 2 x + 1 //ok

    //vector<double> xdata = { 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 2450, 1500 };
    //vector<double> ydata = { 22.69882061, 24.44133150, 26.19332847, 27.95387862, 29.72215460,
    //  31.49742271, 33.27903189, 35.06640374, 36.85902371, 38.65643325, 40.45822296 };
    //vector<double> result = FitPoly(6, xdata.begin(), xdata.end(), ydata.begin());
    //double t = xdata[10];
    //double val(0);
    //for (int i = result.size()-1; i >= 0; --i)
    //  val = (val*t + result[i]); //ok

    return true;
  }

}//namespace

