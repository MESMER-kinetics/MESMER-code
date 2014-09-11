// MolecularComponents.cpp
//
// Author: Chi-Hsiu Liang
//
//-------------------------------------------------------------------------------------------
#include <stdexcept>
#include <numeric>
#include <cmath>
#include <iomanip>
#include "Molecule.h"
#include "System.h"
#include "ParseForPlugin.h"

using namespace std;
using namespace Constants;
using namespace OpenBabel;
namespace mesmer
{

  //-------------------------------------------------------------------------------------------------
  // Bath gas related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gBathProperties::~gBathProperties()
  {
    /*  if (m_Sigma_chk == 0){
    cinfo << "m_Sigma is provided but not used in " << m_host->getName() << "." << endl;
    }
    if (m_Epsilon_chk == 0){
    cinfo << "m_Epsilon is provided but not used in " << m_host->getName() << "." << endl;
    }
    */
  };

  gBathProperties::gBathProperties(Molecule* pMol)
    :m_Sigma(sigmaDefault),
    m_Epsilon(epsilonDefault),
    m_Sigma_chk(-1),
    m_Epsilon_chk(-1)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; //Be forgiving; we can get by without a propertyList element

    setSigma(ppPropList->XmlReadPropertyDouble("me:sigma"));
    setEpsilon(ppPropList->XmlReadPropertyDouble("me:epsilon"));
  }

  void   gBathProperties::setSigma(double value)          {
    m_Sigma = value;
    m_Sigma_chk = 0;
  };

  double gBathProperties::getSigma()                      {
    if (m_Sigma_chk >= 0){
      ++m_Sigma_chk;
      return m_Sigma;
    }
    else{
      cerr << "m_Sigma was not defined but requested in " << m_host->getName() << ". Default value " << sigmaDefault << " is used.\n";
      return m_Sigma;
    }
  };

  void   gBathProperties::setEpsilon(double value)        {
    m_Epsilon = value;
    m_Epsilon_chk = 0;
  };

  double gBathProperties::getEpsilon()                    {
    if (m_Epsilon_chk >= 0){
      ++m_Epsilon_chk;
      return m_Epsilon;
    }
    else{
      cerr << "m_Epsilon was not defined but requested in " << m_host->getName() << ". Default value " << epsilonDefault << " is used.\n";
      return m_Epsilon;
    }
  };

  //-------------------------------------------------------------------------------------------------
  // Cell density of states related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gDensityOfStates::~gDensityOfStates()
  {
    //Delete the density of state calculators because they are cloned instances
    for (unsigned i = 0; i < m_DOSCalculators.size(); ++i)
      delete m_DOSCalculators[i];

    if (m_Hessian)
      delete m_Hessian;
    if (m_Modes)
      delete m_Modes;
  }

  gDensityOfStates::gDensityOfStates(Molecule* pMol)
    :m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(1.0),
    m_ZPE(NaN),
    m_scaleFactor(1.0),
    m_SpinMultiplicity(1),
    m_RC_chk(-1),
    m_Sym_chk(-1),
    m_ZPE_chk(-1),
    m_scaleFactor_chk(-1),
    m_SpinMultiplicity_chk(-1),
    m_EnergyConvention("arbitary"),
    m_eleExc(),
    m_VibFreq(),
    m_Hessian(NULL),
    m_Modes(NULL),
    m_nModes(0),
    m_HessianUnits(),
    m_MaximumCell(0),
    m_cellDOS(),
    m_grainEne(),
    m_grainDOS() {
      m_host = pMol;
  }

  bool gDensityOfStates::initialization() {

    // Define basic molecular parameters. 

    Molecule* pMol = m_host;

    ErrorContext c(pMol->getName());

    PersistPtr pp = pMol->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    // Vibrational frequencies. Test for Hessain first and, if absent,
    // try to read freqeuncies.

    bool hasVibFreq(true);
    const char *txt;
    if (m_Hessian = ReadPropertyMatrix<double>("me:hessian", ppPropList)) {
      PersistPtr pMtrx = ppPropList->XmlMoveToProperty("me:hessian");
      txt = pMtrx->XmlReadValue("units", false);
      m_HessianUnits = (txt) ? string(txt) : "kJ/mol/amu/Ang^2";
      FrqsFromHessian();
    }
    else if (txt = ppPropList->XmlReadProperty("me:vibFreqs", optional)) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_VibFreq.push_back(x);
    }
    else {
      hasVibFreq = false;
      if (!pMol->getStruc().IsAtom())
        cinfo << "Cannot find me:vibFreqs or me:hessian. Assuming an atom or a sink molecule." << endl;
    }

    m_scaleFactor = ppPropList->XmlReadPropertyDouble("me:frequenciesScaleFactor");
    m_scaleFactor_chk = 0;

    // Rotational constants.

    std::vector<double> rCnst(3, 0.0);
    bool hasRotConst(false);
    txt = ppPropList->XmlReadProperty("me:rotConsts", optional);
    if (txt) {
      //data from <me:rotConsts>
      istringstream idata(txt);
      idata >> rCnst[0] >> rCnst[1] >> rCnst[2];
      txt = ppPropList->XmlReadPropertyAttribute("me:rotConsts", "units", optional);
      if (string(txt) == "amuA^2") {
        // The supplied data are moments of inertia with units amuAng^2.
        rCnst[0] = conMntInt2RotCnt / rCnst[0];
        rCnst[1] = conMntInt2RotCnt / rCnst[1];
        rCnst[2] = conMntInt2RotCnt / rCnst[2];
      }
      hasRotConst = true;
      m_RC_chk = 0;
    }
    else {
      // Attempt to calculate rotational constants from atomic coordinates.
      gStructure& gs = pMol->getStruc();
      if (gs.ReadStructure()) {
        rCnst = gs.CalcRotConsts();
        cinfo << "Rotational constants were calculated from atom coordinates: "
          << rCnst[2] << ' ' << rCnst[1] << ' ' << rCnst[0] << " cm-1" << endl;
        hasRotConst = true;
      }
    }
    if (hasRotConst) {
      // Check rotational constants are valid.
      if (rCnst[0] < 0.0 || rCnst[1] < 0.0 || rCnst[2] < 0.0) {
        cerr << "A negative rotational constant found " << endl;
        throw (std::runtime_error("Fatal error"));
      }
      // Make sure the rotational constants are in ascending order.
      std::sort(rCnst.begin(), rCnst.end());
      m_RotCstA = rCnst[2];
      m_RotCstB = rCnst[1];
      m_RotCstC = rCnst[0];
      m_RC_chk = 0;
    }
    else if (!pMol->getStruc().IsAtom()){
      cinfo << "No rotational constants from <me:rotConsts> or structure. "
        "Assuming an atom or a sink molecule." << endl;
    }
    else
      m_RC_chk = 0;

    m_Sym = ppPropList->XmlReadPropertyDouble("me:symmetryNumber");
    m_Sym_chk = 0;

    if (hasVibFreq != hasRotConst){
      cerr << "Improper setting on vibrational frequencies or rotational constants." << endl;
    }

    // Electronic excitations.

    txt = ppPropList->XmlReadProperty("me:electronicExcitation", optional);
    if (txt) {
      istringstream idata(txt);
      m_eleExc.clear();
      double _iele;
      while (idata >> _iele) m_eleExc.push_back(_iele);
    }

    // Spin multiplicity.

    m_SpinMultiplicity = ppPropList->XmlReadPropertyInteger("me:spinMultiplicity", optional);
    if (m_SpinMultiplicity == 0)
      m_SpinMultiplicity = pp->XmlReadInteger("spinMultiplicity");
    m_SpinMultiplicity_chk = 0;

    if (!ReadDOSMethods())
      throw(std::runtime_error(""));

    // Check whether the molecule has the correct number of degrees of freedom.
    if (!m_host->checkDegOfFreedom()){
      string errorMsg = "Incorrect number of degrees of freedom compared with atom count for " + m_host->getName();
      throw (std::runtime_error(errorMsg));
    }

    return ReadZeroPointEnergy(ppPropList);
  }

  bool gDensityOfStates::ReadZeroPointEnergy(PersistPtr &ppPropList) {

    /* For molecular energy me:ZPE is used if it is present. If it is not, a value
    calculated from me:Hf0 or HAT0 is used and a converted value is written back
    to the datafile as a me:ZPE property. It is consequently used in the next run
    and available to be varied or optimized. The original source is recorded in an
    attribute.
    */
    const char* txt;
    string unitsInput;
    double tempzpe = ppPropList->XmlReadPropertyDouble("me:ZPE", optional);
    if (!IsNan(tempzpe))
    {
      // me:ZPE is present
      txt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "units");
      if (txt) //No units specified.
        unitsInput = txt;
      txt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "convention", optional);
      m_EnergyConvention = txt ? txt : "arbitary";

      double zpCorrection = 0.0; //cm-1
      txt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "zeroPointVibEnergyAdded", optional);
      if (txt && strcmp(txt, "false") == 0)
      {
        // Mesmer assumes that the ZPE is the true zero point energy, unless
        // an attribute zeroPointVibEnergyAdded="false" is present. This indicates
        // that the value provided is a computed energy at the bottom of the well
        // and it is corrected here by adding 0.5*Sum(vib freqs).
        if (m_VibFreq.size() > 0)
        {
          zpCorrection = accumulate(m_VibFreq.begin(), m_VibFreq.end(), 0.0);
          zpCorrection *= 0.5;
        }
        //Write back a corrected value and change attribute to zeroPointVibEnergyAdded="true"
        PersistPtr ppScalar = ppPropList->XmlMoveToProperty("me:ZPE");
        ppScalar->XmlWrite(toString(tempzpe + ConvertFromWavenumbers(unitsInput, zpCorrection)));
        ppScalar->XmlWriteAttribute("zeroPointVibEnergyAdded", "true");
      }

      bool rangeSet;
      PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:ZPE");
      ReadRdoubleRange(string(m_host->getName() + ":ZPE"), ppProp, m_ZPE, rangeSet,
        getConvertedEnergy(unitsInput, 1.0), zpCorrection);
      set_zpe(getConvertedEnergy(unitsInput, tempzpe) + zpCorrection);
      m_ZPE_chk = 0;

    }
    else {

      // me:ZPE not present; try me:Hf0 and HAT0 (enthalpy of atomization at 0K)
      // and Hf298 (enthalpy of formation at 298K)
      double Hf0 = ppPropList->XmlReadPropertyDouble("me:Hf0", optional);
      double HfAT0 = ppPropList->XmlReadPropertyDouble("me:HfAT0", optional);
      double Hf298 = ppPropList->XmlReadPropertyDouble("me:Hf298", optional);
      //Do not provide a default ZPE
      //if (IsNan(Hf0) && IsNan(HfAT0) && IsNan(Hf298))
      //{
      //  //None of me:ZPE, me:Hf0, me:HAT0 and meHf298 are present; use a default ZPE.
      //  tempzpe = ppPropList->XmlReadPropertyDouble("me:ZPE", true);
      //}
      //else
      //{
        // If there is not a convention="thermodynamic" attribute on Hf298,
        // convert Hf0,etc. and write back a me:ZPE property which will be used in future
        // Currently, Hf0,etc. cannot be used as a range variables
        /*
        Atomize species X into atoms (at 0K)
        delta H  = Sum(Hf0(atom i)) - Hf0(X)
        = Sum(E(atom i)) - E(X) where E is a compchem energy
        =                - HAT0 (enthalpy of atomization at 0K)
        We have E and Hf0 for each element as gas phase atom in librarymols.xml,
        so E(X) = Hf0(X) + Sum over all atoms( E - Hf0 )
        */
        string origElement = !IsNan(Hf0) ? "me:Hf0" : (!IsNan(HfAT0) ? "me:HfAT0" : "me:Hf298");
        const char* utxt = ppPropList->XmlReadPropertyAttribute(origElement, "units", optional);
        utxt = utxt ? utxt : "kJ/mol";
        if (!IsNan(Hf0))
        {
          //Use Hf0 if provided
          Hf0 = getConvertedEnergy(utxt, Hf0); //cm-1
          tempzpe = Hf0 + getConvertedEnergy("kJ/mol", getHost()->getStruc().CalcSumEMinusHf0(false, false));//cm-1
        }
        else if (!IsNan(HfAT0))
        {
          //Use HfAT0 (atom-based thermochemistry, see DOI: 10.1002/chem.200903252)
          HfAT0 = getConvertedEnergy(utxt, HfAT0); //cm-1
          tempzpe = HfAT0 + getConvertedEnergy("kJ/mol", getHost()->getStruc().CalcSumEMinusHf0(true, false));//cm-1
        }
        else
        {
          if(!IsNan(Hf0))
          {
            //Use Hf298
          const char* convention = ppPropList->XmlReadPropertyAttribute(origElement, "convention", optional);
          if (convention && strcmp(convention, "thermodynamic") == 0)
            {
              m_EnergyConvention = convention;
              tempzpe = getConvertedEnergy(utxt, Hf298); //Raw Hf298 to cm-1
              return calcDensityOfStates(); //necessary here for Unit Tests but I don't know why.
            }
            /*Atomize species X at 298K
            deltaH  = Sum over atoms(Hf298)) - Hf298(X)
            = Sum(E + Sum(H(298K)) - (E(X) + H(298K))
            E(X) = (Hf298 - H(298K))(X) + Sum over atoms(E - Hf298 + H(298K))
            */
            // H is the enthalpy in cm-1 and 298K calculated with m_ZPE=0.
            // Hf0 is the real enthalpy of formation at 0K in cm-1.
            double H, S, G;
            set_zpe(0.0);
            thermodynamicsFunctions(298, 1.0, H, S, G);
            //Hf298 = getConvertedEnergy(utxt, Hf298) + H; //cm-1 sign changed
            tempzpe = getConvertedEnergy(utxt, Hf298) - H
              + getConvertedEnergy("kJ/mol", getHost()->getStruc().CalcSumEMinusHf0(false, true));//cm-1        
          }
        }
        if(IsNan(tempzpe))
        {
          cinfo << "No ZPE (or an alternative) was provided." << endl;
          return true;
        }
        set_zpe(tempzpe);
        m_ZPE_chk = 0;
        //Write the converted value back to a me:ZPE element in the XML file
        stringstream ss;
        ss.precision(9);
        ss << ConvertFromWavenumbers(utxt, tempzpe);
        PersistPtr ppScalar = ppPropList->XmlWriteProperty("me:ZPE", ss.str(), utxt);
        ppScalar->XmlWriteAttribute("source", origElement);
        ppScalar->XmlWriteAttribute("convention", "computational");//orig units
        m_EnergyConvention = "computational";
        cinfo << "New me:ZPE element written with data from " << origElement << endl;
      //}
    }
    return true;
  }

  //
  // Get the number of degrees of freedom for this species.
  //
  unsigned int gDensityOfStates::getNoOfDegOfFreeedom() {
    unsigned int nDOF(0);
    for (vector<DensityOfStatesCalculator*>::size_type j = 0; j < m_DOSCalculators.size(); ++j) {
      nDOF += m_DOSCalculators[j]->NoDegOfFreedom(this);
    }
    return nDOF;
  }

  //
  // Get cell density of states.
  //
  bool gDensityOfStates::getCellDensityOfStates(vector<double> &cellDOS, int startingCell, bool bcalc) {
    // If density of states have not already been calcualted then do so.
    if (bcalc && !calcDensityOfStates())
    {
      cerr << "Failed calculating DOS" << endl;
      return false;
    }
    if (startingCell == 0)
      cellDOS = m_cellDOS;
    else{
      int MaximumCell = m_host->getEnv().MaxCell;
      for (int i(startingCell); i < MaximumCell; ++i){
        cellDOS.push_back(m_cellDOS[i]);
      }
    }
    return true;
  }

  bool gDensityOfStates::ReadDOSMethods() {
    // There must be only one <me:DOSCMethod> specifying a method which
    // includes the rotations. If not present the default is used.
    // There can be multiple additional methods:
    // <me:ExtraDOSCMethod name="..."> optional parameters </me:ExtraDOSCMethod>

    ErrorContext c(getHost()->getName());
    PersistPtr pp = getHost()->get_PersistentPointer();

    //Get the method which includes rotations, or use the default
    const char* name = pp->XmlReadValue("me:DOSCMethod");
    if (!*name) // Must be alt form e.g. <me:DOSCMethod name="QMRotors"/>
      name = pp->XmlMoveTo("me:DOSCMethod")->XmlReadValue("name");

    m_DOSCalculators.push_back(DensityOfStatesCalculator::Find(string(name)));
    if (!m_DOSCalculators[0] || !m_DOSCalculators[0]->ReadParameters(this, pp))
      return false; //error message already output
    if (!m_DOSCalculators[0]->includesRotations())
    {
      cerr << "The calculator specified in <me:DOSCMethod>"
        " should be one that includes the rotations."
        " Use  <me:ExtraDOSCMethod name=\""
        << m_DOSCalculators[0]->getID() << "\"> </me:ExtraDOSCMethod> instead."
        << endl;
      return false;
    }

    // Beyer-Swinehart object added by default at m_DOSCalculators[1]
    m_DOSCalculators.push_back(DensityOfStatesCalculator::Find("BeyerSwinehart"));
    if (!m_DOSCalculators[1])
    {
      cerr << "Beyer-Swinehart algorithm failed to initialize correctly" << endl;
      return false;
    }

    //Read any additional methods
    PersistPtr pp2 = pp;
    while (pp2 = pp2->XmlMoveTo("me:ExtraDOSCMethod"))
    {
      string dosMethod;
      const char* name;
      if ((name = pp2->XmlRead())
        || (name = pp2->XmlReadValue("name", optional))
        || (name = pp2->XmlReadValue("xsi:type", optional)))
      {
        if (name[2] == ':')
          name += 3; //Remove prefix "me:"
        dosMethod = name;
      }

      DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find(dosMethod);
      if (!pDOSCalculator || !pDOSCalculator->ReadParameters(this, pp2))
        return false;
      m_DOSCalculators.push_back(pDOSCalculator);
    }

    //Check there is only one <me:DOSCMethod>
    PersistPtr pp1 = pp->XmlMoveTo("me:DOSCMethod"); //to the element just found
    if (pp1->XmlMoveTo("me:DOSCMethod"))
      cerr << "Too many <me:DOSCMethod> elements on this molecule. "
      << "Only the first is used. Additional methods should be under <me:ExtraDOSCMethod>." << endl;
    return true;
  }

  bool gDensityOfStates::RemoveDOSCalculator(const string& id)
  {
    vector<DensityOfStatesCalculator*>::iterator iter;
    for (iter = m_DOSCalculators.begin(); iter != m_DOSCalculators.end(); ++iter)
    {
      if (id == (*iter)->getID())
      {
        delete *iter; //because plugin was a new instance made with Clone()
        m_DOSCalculators.erase(iter);
        return true;
      }
    }
    return false;
  }
  bool gDensityOfStates::AddDOSCalculator(const string& id)
  {
    m_DOSCalculators.push_back(DensityOfStatesCalculator::Find(id));
    return m_DOSCalculators.back();
  }

  DensityOfStatesCalculator* gDensityOfStates::GetDOSCalculator(const string& id)
  {
    vector<DensityOfStatesCalculator*>::iterator iter;
    for (iter = m_DOSCalculators.begin(); iter != m_DOSCalculators.end(); ++iter)
    {
      if (id == (*iter)->getID())
        return *iter;
    }
    return NULL;
  }

  RotationalTop gDensityOfStates::test_rotConsts()
  {
    std::vector<double> mmtsInt;
    return get_rotConsts(mmtsInt);
  }

  RotationalTop gDensityOfStates::get_rotConsts(std::vector<double> &mmtsInt)
  {
    if (m_RC_chk <= -1){
      ErrorContext e(this->getHost()->getName());
      if (m_RC_chk == -1)
//        cinfo << "Rotational constants were not defined but requested." << endl;
      --m_RC_chk;
      return UNDEFINED_TOP; // treat as a non-rotor
    }

    mmtsInt.clear();
    mmtsInt.push_back(m_RotCstA);
    mmtsInt.push_back(m_RotCstB);
    mmtsInt.push_back(m_RotCstC);
    ++m_RC_chk;

    // The classification of rotors is simplified to only three following types.
    // 3-D rotors may have other attributes, but in ILT they are treated as the same type. 

    if ((mmtsInt[0] + mmtsInt[1] + mmtsInt[2]) == 0.)
      return UNDEFINED_TOP; // not a rotor
    else if ((mmtsInt[0] * mmtsInt[1] * mmtsInt[2]) == 0.)
      return LINEAR;        // 2-D linear
    else
      return NONLINEAR;     // 3-D symmetric/asymmetric/spherical top
  }


  //
  // Calculate the rovibrational density of states.
  //
  bool gDensityOfStates::calcDensityOfStates()
  {
    bool recalc(false);
    const size_t MaximumCell = m_host->getEnv().MaxCell;
    if (MaximumCell > m_MaximumCell) {
      recalc = true ;
      m_MaximumCell = MaximumCell ;
    }

    if (recalc) {
      // Calculate density of states.
      bool ret(true);
      for (size_t i(0); ret && i < m_DOSCalculators.size(); ++i)
        ret = ret && m_DOSCalculators[i]->countCellDOS(this, m_host->getEnv());
      if (!ret)
        return false;
    }

    if(IsNan(m_ZPE))
    {
      //cinfo << "calculation of DOS cutailed because no ZPE" << endl;
    }
    else
    {
      const int cellOffset = get_cellOffset();
      std::vector<double> cellEne;
      getCellEnergies(MaximumCell, m_host->getEnv().CellSize, cellEne);
      calcGrainAverages(m_host->getEnv().MaxGrn, m_host->getEnv().cellPerGrain(), cellOffset, m_cellDOS, cellEne, m_grainDOS, m_grainEne);
    }

    if (recalc) {
      testDensityOfStates();
    }

    recalculateDOScompleted();

    return true;
  }

  // Calculate classical energy
  double gDensityOfStates::getClassicalEnergy(){
    //Basically use the frequencies to calculate the contribution of ZPE from harmonic oscillators approximation
    double ZC = 0.0;
    for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
      ZC += m_VibFreq[i] / 2.0;
    return get_zpe() - ZC;
  }


  //
  // Test the rovibronic density of states
  //
  void gDensityOfStates::testDensityOfStates()
  {
    const int MaximumGrain = m_host->getEnv().MaxGrn;
    const int MaximumCell = m_host->getEnv().MaxCell;
    std::vector<double> cellEne;
    getCellEnergies(MaximumCell, m_host->getEnv().CellSize, cellEne);

    // Partition functions that are higher than the current simulation temperature will not be output.
    const double temperature = 1. / (boltzmann_RCpK * m_host->getEnv().beta);
    const int max_nplus1 = int(temperature / 100.);

    if (m_host->isMolType("modelled") || m_host->isMolType("transitionState")){
      string comment("Rovibronic partition function calculation at various temperatures. qtot : product of QM partition functions for vibrations (1-D harmonic oscillator) and classical partition functions for rotations.  sumc : cell based partition function. sumg : grain based partition function ");

      PersistPtr ppList = m_host->get_PersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment);

      if (m_host->getFlags().testDOSEnabled) ctest << endl << "Test rovibronic density of states for: " << m_host->getName() << "\n{\n";
      if (m_host->getFlags().testDOSEnabled) ctest << "      T           qtot           sumc           sumg\n";


      //loop through predefined test temperatures
      for (int n = 0; n < max_nplus1; ++n) {
        double temp = 100.0*static_cast<double>(n + 2);
        double beta = 1.0 / (boltzmann_RCpK*temp);

        // Calculate rovibronic partition functions based on cells.
        double cellCanPrtnFn = canonicalPartitionFunction(m_cellDOS, cellEne, beta);

        // Calculate rovibronic partition functions based on grains.
        double grainCanPrtnFn = canonicalPartitionFunction(m_grainDOS, m_grainEne, beta);

        // Calculate rovibronic partition functions, using analytical formula where possible.
        double qtot(1.0);
        for (vector<DensityOfStatesCalculator*>::size_type j = 0; j < m_DOSCalculators.size(); ++j) {
          qtot *= m_DOSCalculators[j]->canPrtnFnCntrb(this, beta);
        }

        if (m_host->getFlags().testDOSEnabled) {
          formatFloat(ctest, temp, 6, 7);
          formatFloat(ctest, qtot, 6, 15);
          formatFloat(ctest, cellCanPrtnFn, 6, 15);
          formatFloat(ctest, grainCanPrtnFn, 6, 15);
          ctest << endl;
        }

        //Add to XML document
        PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
        ppItem->XmlWriteValueElement("me:T", temp, 6);
        ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
        ppItem->XmlWriteValueElement("me:sumc", cellCanPrtnFn, 6);
        ppItem->XmlWriteValueElement("me:sumg", grainCanPrtnFn, 6);
      }
      if (m_host->getFlags().testDOSEnabled) ctest << "}" << endl;
    }

    if (m_host->getFlags().cellDOSEnabled){
      ctest << endl << "Cell rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i(0); i < MaximumCell; ++i){
        formatFloat(ctest, cellEne[i], 6, 15);
        formatFloat(ctest, m_cellDOS[i], 6, 15);
        ctest << endl;
      }
      ctest << "}" << endl;
    }

    if (m_host->getFlags().grainDOSEnabled && (m_host->isMolType("modelled") || m_host->isMolType("transitionState"))){
      ctest << endl << "Grain rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i(0); i < MaximumGrain; ++i){
        formatFloat(ctest, m_grainEne[i], 6, 15);
        formatFloat(ctest, m_grainDOS[i], 6, 15);
        ctest << endl;
      }
      ctest << "}" << endl;
    }
  }

  double gDensityOfStates::get_zpe() {
    if (m_ZPE_chk == -1) {
//      cinfo << "m_ZPE was not defined but requested in " << m_host->getName() << ". Default value " << m_ZPE << " is given." << endl;
      --m_ZPE_chk;
    }
    else if (m_ZPE_chk < -1) {
      --m_ZPE_chk;
    }
    else {
      ++m_ZPE_chk;
    }
    return double(m_ZPE);
  }

  double gDensityOfStates::get_scaleFactor() {
    if (m_scaleFactor_chk == -1){
      cinfo << "m_scaleFactor was not defined but requested in " << m_host->getName() << ". Default value " << m_scaleFactor << " is given." << endl;
      --m_scaleFactor_chk;
    }
    else if (m_scaleFactor_chk < -1) {
      --m_scaleFactor_chk;
    }
    else {
      ++m_scaleFactor_chk;
    }
    return m_scaleFactor;
  }

  double gDensityOfStates::get_Sym(void){
    if (m_Sym_chk == -1){
      cinfo << "m_Sym was not defined but requested in " << m_host->getName() << ". Default value " << m_Sym << " is given." << endl;
      --m_Sym_chk;
    }
    else if (m_Sym_chk < -1) {
      --m_Sym_chk;
    }
    else {
      ++m_Sym_chk;
    }
    return m_Sym;
  }

  void gDensityOfStates::get_VibFreq(std::vector<double>& vibFreq){
    const double scalefactor = get_scaleFactor();
    for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
      vibFreq.push_back(m_VibFreq[i] * scalefactor);
  }

  bool gDensityOfStates::removeVibFreq(double freq) {
    vector<double>::iterator pos = find(m_VibFreq.begin(), m_VibFreq.end(), freq);
    if (pos == m_VibFreq.end())
      return false;
    m_VibFreq.erase(pos);
    return true;
  }

  int gDensityOfStates::getSpinMultiplicity(){
    if (m_SpinMultiplicity_chk >= 0) {
      ++m_SpinMultiplicity_chk;
    }
    else {
      cinfo << "m_SpinMultiplicity was not defined but requested in " << m_host->getName() << ". Default value " << m_SpinMultiplicity << " is given." << endl;
    }
    return m_SpinMultiplicity;
  }

  int gDensityOfStates::get_cellOffset(void) {
    double modulus = fmod(get_zpe() - m_host->getEnv().EMin, double(m_host->getEnv().GrainSize))/m_host->getEnv().CellSize ;
    return int(max(modulus,0.0));
  };

  //
  // Get grain density of states.
  //
  void gDensityOfStates::getGrainDensityOfStates(vector<double> &grainDOS, const int startGrnIdx, const int ignoreCellNumber) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates()){
      throw (std::runtime_error("Failed calculating DOS."));
    }
    if (ignoreCellNumber == 0){ // If there is no cells ignored in this grain, the grain DOS dose not need to be recalculated.
      grainDOS = m_grainDOS;
    }
    else{ // Some cells are ignored in this grain, as they do not occur in this part of reaction.
      // first deal with the first grain.
      // const int MaximumCell = m_host->getEnv().MaxCell;
      const int gsz = m_host->getEnv().GrainSize;
      const int cellOffset = get_cellOffset();
      const int grnStartCell = startGrnIdx * gsz - cellOffset;
      double partialDOS(0.0);
      for (int i(ignoreCellNumber); i < gsz; ++i){
        partialDOS += m_cellDOS[i + grnStartCell];
      }
      grainDOS.clear();
      grainDOS.push_back(partialDOS);
      for (int i(startGrnIdx + 1); i < int(m_grainDOS.size()); ++i){
        grainDOS.push_back(m_grainDOS[i]);
      }
    }
  }

  //
  // Get grain energies.
  //
  void gDensityOfStates::getGrainEnergies(vector<double> &grainEne) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";
    grainEne = m_grainEne;
  }

  //
  // Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
  //
  double gDensityOfStates::rovibronicGrnCanPrtnFn() {

    // If density of states have not already been calculated then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";

    return canonicalPartitionFunction(m_grainDOS, m_grainEne, m_host->getEnv().beta);

  }

  //
  // Calculate standard thermodynamic quantities as a function of temperature.
  // The calculation is based on cell densities of states.
  //
  bool gDensityOfStates::thermodynamicsFunctions(double temp, double unitFctr, double& enthalpy, double& entropy, double& gibbsFreeEnergy) {

    std::vector<double> cellEne;
    getCellEnergies(m_host->getEnv().MaxCell, m_host->getEnv().CellSize, cellEne);

    calcDensityOfStates();

    double beta;
    if (temp > 0.0) {
      beta = 1.0 / (boltzmann_RCpK*temp);
    }
    else {
      return false;
    }

    // Calculate rovibronic partition functions based on cells.
    double cellCanPrtnFn = canonicalPartitionFunction(m_cellDOS, cellEne, beta);

    // The following calculates the mean internal molecular energy.
    double internalEnergy = canonicalMeanEnergy(m_cellDOS, cellEne, beta);

    // The rovibronic partition function must be corrected for translation 
    // and (assuming an ideal gas) molecular indistinguishability.
	double molarVol = 1.e+06*idealGasC/(boltzmann_RCpK*beta*atm_in_pascal) ;  // cm3
    gibbsFreeEnergy = unitFctr*(-log(cellCanPrtnFn) 
	                - log(tp_C * pow((m_host->getStruc().getMass() / beta), 1.5)*molarVol) 
					+ log(AvogadroC)) / beta;

    // The enthalpy must be corrected for translation by an additional 3kT/2.
    enthalpy = unitFctr*(internalEnergy + 5.0 / (2.0*beta));

    entropy = (enthalpy - gibbsFreeEnergy) / temp;

    return true;
  }

  // Calculate vibrational frequencies from molecular Hessian. This method 
  // projects out the overall translation and rotation vectors as defined 
  // in Wilson, Decius and Cross. Molecular Vibrations, Dover 1980. See
  // also Miller, Handy and Adams JCP, Vol, 72, 99 (1980).
  bool gDensityOfStates::FrqsFromHessian() {

    const size_t msize = m_Hessian->size();
    m_Modes = new dMatrix(msize, 0.0);

    gStructure& gs = m_host->getStruc();
    bool HasCoords = gs.ReadStructure();
    if (!HasCoords || 3 * size_t(gs.NumAtoms()) != msize) {
      cerr << "The dimension of the defined Hessian for " << getHost()->getName() << " does not match the specified number of atoms.";
      return false;
    }

    // Get atomic massess and coordinates.

    vector<double> atomicMasses, xx, yy, zz;
    gs.getAtomicMasses(atomicMasses);
    gs.getXCoords(xx);
    gs.getYCoords(yy);
    gs.getZCoords(zz);

    // Initialize mass weights.

    vector<double> massWeights;
    size_t i(0), j(0);
    for (; i < atomicMasses.size(); i++) {
      double weight = sqrt(atomicMasses[i]);
      for (j = 0; j < 3; j++) {
        massWeights.push_back(weight);
      }
    }

    // Mass weight Hessian.

    double convFactor(1.0);
    if (m_HessianUnits == "kJ/mol/Ang2") {
      // Nothing to do.
    }
    else if (m_HessianUnits == "kcal/mol/Ang2") {
      convFactor *= Calorie_in_Joule;
    }
    else if (m_HessianUnits == "Hartree/Bohr2") {
      convFactor *= Hartree_In_kJperMol / (bohr_in_angstrom * bohr_in_angstrom);
    }
    else {
      throw (std::runtime_error("Unknown Hessian units."));
    }

    for (i = 0; i < msize; i++) {
      for (j = i; j < msize; j++) {
        (*m_Hessian)[i][j] *= convFactor / (massWeights[i] * massWeights[j]);
        (*m_Hessian)[j][i] = (*m_Hessian)[i][j];
      }
    }

    // Rotate Hessian.

    dMatrix axisAlignment(3, 0.0);
    gs.getAlignmentMatrix(axisAlignment);
    dMatrix Transform(msize, 0.0);
    for (i = 0; i < msize; i += 3) {
      size_t ii(0), jj(0);
      for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {
          Transform[i + ii][i + jj] = axisAlignment[ii][jj];
        }
      }
    }

    dMatrix tmpHessian = (*m_Hessian)*Transform;
    Transform.Transpose();
    *m_Hessian = Transform*tmpHessian;

    // X Translation projector.

    vector<double> mode(msize, 0.0);
    for (i = 0; i < msize; i += 3)
      mode[i] = massWeights[i];

    UpdateProjector(mode);

    // Y Translation projector.

    ShiftTransVector(mode);
    UpdateProjector(mode);

    // Z Translation projector.

    ShiftTransVector(mode);
    UpdateProjector(mode);

    // Rotational modes.
    RotationVector(yy, 2, 1.0, zz, 1, -1.0, massWeights, mode);
    UpdateProjector(mode);

    RotationVector(xx, 2, -1.0, zz, 0, 1.0, massWeights, mode);
    UpdateProjector(mode);

    RotationVector(xx, 1, 1.0, yy, 0, -1.0, massWeights, mode);
    UpdateProjector(mode);

    // Project out translational and rotational modes.

    vector<double> freqs(msize, 0.0);
    calculateFreqs(freqs, m_host->isMolType("transitionState"));

    // Save projected frequencies. Note need to check if configuration is a transtion state.

    size_t firstFrqIdx = (m_host->isMolType("transitionState")) ? 7 : 6;
    m_VibFreq.clear();
    for (i = firstFrqIdx; i < msize; i++)
      m_VibFreq.push_back(freqs[i]);

    return true;
  }

  // Function to calculate the vibrational frequencies from a projected Hessian matrix.
  bool gDensityOfStates::calculateFreqs(vector<double> &freqs, bool projectTransStateMode) {

    const size_t msize = m_Hessian->size();

    dMatrix tmp = *m_Modes;
    tmp.Transpose();
    dMatrix Projector = (*m_Modes)*tmp;
    for (size_t i = 0; i < msize; i++) {
      for (size_t j = 0; j < msize; j++) {
        Projector[i][j] *= -1.0;
      }
      Projector[i][i] += 1.0;
    }

    dMatrix tmpHessian = Projector*(*m_Hessian)*Projector;

    tmpHessian.diagonalize(&freqs[0]);

    double convFactor = conHess2Freq / (2.0*M_PI);
    for (size_t m(0); m < msize; m++) {
      if (freqs[m] > 0.0) {

        // Mostly vibration modes.
        freqs[m] = convFactor*sqrt(freqs[m]);

      }
      else if (projectTransStateMode) {

        // Add Transition state mode to the projected set after orthogonalization.
        // The magic number of "1.0" below is 1 cm-1 and used to filter out any 
        // small -ve frequencies associated with existing projected modes.
        double imFreq = convFactor*sqrt(fabs(freqs[m]));
        if (imFreq > 1.0) {

          vector<double> mode(msize, 0.0);
          for (size_t j(0); j < msize; j++) {
            mode[j] = tmpHessian[j][m];
          }

          orthogonalizeMode(mode);
        }
      }
    }

    return true;
  }

  // This method is used to project a mode from the stored Hessian and
  // re-calculate the remaining frequencies.
  bool gDensityOfStates::projectMode(vector<double> &mode) {

    bool status(true);

    const size_t msize = m_Hessian->size();

    gStructure& gs = m_host->getStruc();
    vector<double> atomicMasses;
    gs.getAtomicMasses(atomicMasses);

    size_t i, j, k;
    vector<double> massWeights(msize, 0.0);
    for (j = 0, i = 0; j < atomicMasses.size(); j++) {
      double weight = sqrt(atomicMasses[j]);
      for (k = 0; k < 3; k++, i++) {
        mode[i] *= weight;
      }
    }

    // Orthogonalize and project out mode.

    orthogonalizeMode(mode);

    vector<double> freqs(msize, 0.0);
    calculateFreqs(freqs);

    // Save projected frequencies.

    size_t nfreq = m_VibFreq.size();
    m_VibFreq.clear();
    for (i = msize - nfreq + 1; i < msize; i++)
      m_VibFreq.push_back(freqs[i]);

    return status;
  }

  // This method is used to orthogonalize a mode against existing
  // projected modes and then add it to the projected set.
  bool gDensityOfStates::orthogonalizeMode(vector<double> &mode) {

    const size_t msize = m_Hessian->size();

    // Orthogonalize against existing modes. 

    size_t i(0), j(0);
    for (i = 0; i < m_nModes; i++) {
      double sum(0.0);
      for (j = 0; j < msize; j++) {
        sum += mode[j] * (*m_Modes)[j][i];
      }
      for (j = 0; j < msize; j++) {
        mode[j] -= sum * (*m_Modes)[j][i];
      }
    }

    UpdateProjector(mode);

    return true;
  }

  // Helper function to shift translation projection vector.
  void gDensityOfStates::ShiftTransVector(vector<double> &mode) {
    const size_t msize = mode.size();
    for (size_t i = msize - 1; i > 0; i--)
      mode[i] = mode[i - 1];
    mode[0] = 0.0;
  }

  // Helper function to create projector.
  void gDensityOfStates::UpdateProjector(vector<double> &mode) {

    // Normalize mode.

    double NormFctr(0.0);
    size_t i(0);
    for (; i < mode.size(); i++) {
      NormFctr += mode[i] * mode[i];
    }

    NormFctr = 1.0 / sqrt(NormFctr);
    for (i = 0; i < mode.size(); i++) {
      mode[i] *= NormFctr;
    }

    // Add mode to existing set.

    for (i = 0; i < m_Modes->size(); i++)
      (*m_Modes)[i][m_nModes] = mode[i];
    m_nModes++;

  }

  // Function to calculate the rotational mode vectors.
  void gDensityOfStates::RotationVector(vector<double> &aa, size_t loca, double sgna, vector<double> &bb, size_t locb, double sgnb, vector<double> &massWeights, vector<double> &mode) {

    mode.clear();
    size_t ncoords = aa.size();
    size_t i(0);
    for (; i < ncoords; i++) {
      vector<double> rcross(3, 0.0);
      rcross[loca] = sgna*aa[i];
      rcross[locb] = sgnb*bb[i];
      for (size_t j(0); j < 3; j++) {
        mode.push_back(rcross[j]);
      }
    }

    // Mass weight vector ;
    for (i = 0; i < mode.size(); i++) {
      mode[i] *= massWeights[i];
    }

  }

  //-------------------------------------------------------------------------------------------------
  // Transition state related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gTransitionState::~gTransitionState()
  {
    //if (m_ImFreq_chk == 0) cinfo << "m_ImFreq is provided but not used in " << m_host->getName() << "." << endl;
  }

  gTransitionState::gTransitionState(Molecule* pMol) :m_ImFreq(0.0),
    m_ImFreq_chk(-1)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;
    txt = ppPropList->XmlReadProperty("me:imFreqs", optional);
    if (!txt) {
      cinfo << "No imaginary vibrational frequency. Check if there is tunneling.\n";
      m_ImFreq_chk = -1;
      // Note: If there is tunneling, an imaginary frequency must also be supplied.
    }
    else {
      istringstream idata(txt); double x;
      // make sure this number is positive
      idata >> x;
      set_imFreq(x);
      bool rangeSet;
      PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:imFreqs");
      ReadRdoubleRange(string(m_host->getName() + ":imFreqs"), ppProp, m_ImFreq, rangeSet);
      m_ImFreq_chk = 0;
    }
  }

  double gTransitionState::get_ImFreq(){
    if (m_ImFreq_chk == -1){
      cerr << "m_ImFreq was not defined but requested." << ". Default value " << m_ImFreq << " is given.";
      --m_ImFreq_chk;
      return m_ImFreq;
    }
    else if (m_ImFreq_chk < -1){
      --m_ImFreq_chk;
      return m_ImFreq;
    }
    ++m_ImFreq_chk;
    return m_ImFreq;
  }


  //-------------------------------------------------------------------------------------------------
  // Population and equilibrium fraction
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //

  // Destructor and initialization, not required.
  // gPopulation::~gPopulation();

  gPopulation::gPopulation(Molecule* pMol)
    :m_initPopulation(0.0),
    m_eqFraction(0.0)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
  }

	gStructure::gStructure(mesmer::Molecule *pMol) : m_MolecularWeight(-1),
    m_PrincipalMI(3, 0.0),
    m_AxisAlignment(NULL),
    Atoms(),
    Bonds(),
    m_atomicOrder(),
    m_HasCoords(false),
    m_verbose(false),
    m_RotBondIDs()
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; //Be forgiving; we can get by without a propertyList element
    double MW = ppPropList->XmlReadPropertyDouble("me:MW", optional);
    if (IsNan(MW))
    {
      ReadStructure();
      if (Atoms.empty())
        cinfo << "Neither chemical structure nor "
        "Molecular Weight as an XML property were provided." << endl;
      else
        MW = CalcMW();
    }
    setMass(MW);
  }

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS
    (gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, std::vector<double>& rctsCellDOS){
      std::vector<double> rct1CellDOS;
      std::vector<double> rct2CellDOS;
      if (!pDOS1.getCellDensityOfStates(rct1CellDOS) || !pDOS2.getCellDensityOfStates(rct2CellDOS))
        return false;
      std::vector<double> rotConsts;
      if (pDOS1.get_rotConsts(rotConsts) == UNDEFINED_TOP){
        rctsCellDOS = rct2CellDOS;
      }
      else if (pDOS2.get_rotConsts(rotConsts) == UNDEFINED_TOP){
        rctsCellDOS = rct1CellDOS;
      }
      else{
        FastLaplaceConvolution(rct1CellDOS, rct2CellDOS, rctsCellDOS);
      }
      return true;
  }

  //Returns true if atoms have coordinates
  bool gStructure::ReadStructure()
  {
    if (!Atoms.empty())
      return m_HasCoords;
    PersistPtr ppMol = getHost()->get_PersistentPointer();
    PersistPtr ppAtom = ppMol->XmlMoveTo("atomArray");
    if (!ppAtom) // there may not be an <atomArray> element
      ppAtom = ppMol;
    while (ppAtom = ppAtom->XmlMoveTo("atom"))
    {
      atom at;
      const char* el = ppAtom->XmlReadValue("elementType");
      if (!el)
      {
        cerr << "<atom> elements must have an elementType attribute" << endl;
        return false;
      }
      at.element = el;
      const char* pId = ppAtom->XmlReadValue("id", optional);
      at.id = (pId) ? pId : at.element;
      double x3, y3, z3;
      x3 = ppAtom->XmlReadDouble("x3", optional);
      y3 = ppAtom->XmlReadDouble("y3", optional);
      z3 = ppAtom->XmlReadDouble("z3", optional);
      if (!IsNan(x3) && !IsNan(y3) && !IsNan(z3))
      {
        at.coords.Set(x3, y3, z3);
        if (x3 != 0 || y3 != 0 || z3 != 0)
          m_HasCoords = true; //at least one atom with non-zero coordinates
      }
      Atoms[at.id] = at;
      m_atomicOrder.push_back(at.id);
    }

    // If coordinates are defined, shift coordinates to the centre of mass/principal axis frame.
    if (m_HasCoords)
      AlignCoords();

    //Read all the bonds. For each bond add a connect to each atom
    PersistPtr ppBond = ppMol->XmlMoveTo("bondArray");
    int ibond = 1;
    if (!ppBond) // there may not be an <bondArray> element
      ppBond = ppMol;
    while (ppBond = ppBond->XmlMoveTo("bond"))
    {
      const char* pId = ppBond->XmlReadValue("id", optional);
      string id;
      if (pId)
        id = pId;
      else
      {
        //id is e.g. "bond3", if not provided
        stringstream ss;
        ss << " bond" << ibond;
        id = ss.str();
      }

      const char* pRefs = ppBond->XmlReadValue("atomRefs2");
      if (!pRefs) return false;
      string refs(pRefs);
      string::size_type pos = refs.find_first_of(" ,");
      string::size_type pos2 = refs.find_last_of(" ,");
      if (pos == string::npos) return false;
      string atomref1 = refs.substr(0, pos);
      string atomref2 = refs.substr(pos2 + 1);
      Bonds[id] = make_pair(atomref1, atomref2);
      Atoms[atomref1].connects.push_back(atomref2);
      Atoms[atomref2].connects.push_back(atomref1);
      ++ibond;
    }

    return m_HasCoords;
  }

  // Method to shift coordinates to the centre of mass/principal axis frame. 
  bool gStructure::AlignCoords() {

    // Calculate centre of mass and subtract from coordinates, and calculate mass weights.

    bool status(true);

    //Determine centre of mass
    map<string, atom>::iterator iter;
    vector3 centreOfMass;
    double mt = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      double mass = atomMass(iter->second.element);
      centreOfMass += iter->second.coords * mass;
      mt += mass;
    }
    centreOfMass /= mt;

    dMatrix MI(3);
    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter)
    {
      vector3 c = iter->second.coords - centreOfMass;
      iter->second.coords.Set(c.x(), c.y(), c.z());
      double  m = atomMass(iter->second.element);
      sxx += m * c.x() * c.x();
      syy += m * c.y() * c.y();
      szz += m * c.z() * c.z();
      sxy += m * c.x() * c.y();
      sxz += m * c.x() * c.z();
      syz += m * c.y() * c.z();
    }

    if (NumAtoms() == 2)
      m_PrincipalMI[0] = szz;
    else {
      MI[0][0] = syy + szz;
      MI[1][1] = sxx + szz;
      MI[2][2] = sxx + syy;
      MI[0][1] = MI[1][0] = -sxy;
      MI[0][2] = MI[2][0] = -sxz;
      MI[1][2] = MI[2][1] = -syz;

      MI.diagonalize(&m_PrincipalMI[0]);

      // Save the Alignment matrix 

      m_AxisAlignment = new dMatrix(MI);

      // Rotate coordinates to principal axis frame.

      MI.Transpose();
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        vector<double> r(3, 0.0);
        iter->second.coords.Get(&r[0]);
        r *= MI;
        iter->second.coords.Set(&r[0]);
      }
    }

    return status;
  }

  double gStructure::CalcMW()
  {
    map<string, atom>::iterator iter;
    double MW = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter)
      MW += atomMass(iter->second.element);
    return MW;
  }

  // Returns in atomset the IDs of all the atoms attached to atomID via bonds,
  // but does not include any atoms already in atomset or atoms beyond them.
  // Handles rings. (Recursive function) 
  void gStructure::GetAttachedAtoms(vector<string>& atomset, const string& atomID)
  {
    atomset.push_back(atomID);
    vector<string>::iterator coniter;
    for (coniter = Atoms[atomID].connects.begin(); coniter != Atoms[atomID].connects.end(); ++coniter)
    {
      if (find(atomset.begin(), atomset.end(), *coniter) != atomset.end())
        continue;
      GetAttachedAtoms(atomset, *coniter);
    }
  }

  double gStructure::CalcMomentAboutAxis(vector<string> atomset, vector3 at1, vector3 at2)
  {
    double sumMoment = 0.0;
    vector<string>::iterator iter;
    for (iter = atomset.begin(); iter != atomset.end(); ++iter)
    {
      vector3 a = Atoms[*iter].coords;
      double d = Point2Line(a, at1, at2);
      sumMoment += atomMass(Atoms[*iter].element) * d * d;
    }
    if (sumMoment > 0.0){
      return sumMoment;
    }
    else {
      string errorMsg = "Problem with calculation of reduced moment of inertia. Either there is something wrong with the atom co-ordinates or the atomic symbol is not present in unitsConversion.cpp ";
      throw (std::runtime_error(errorMsg));
    }
  }

  // Calculates internal rotation eigenvector about an axis define by at1 and at2.
  bool gStructure::CalcInternalRotVec(vector<string> atomset, vector3 at1, vector3 at2, vector<double> &mode, bool ApplyMWeight)
  {
    vector3 diff = at1 - at2;
    diff.normalize();
    vector<string>::iterator iter;
    for (iter = atomset.begin(); iter != atomset.end(); ++iter)
    {
      double massWeight = (ApplyMWeight) ? sqrt(atomMass(Atoms[*iter].element)) : 1.0;
      vector3 a = Atoms[*iter].coords - at1;
      vector3 b = cross(a, diff);
      int atomicOrder = getAtomicOrder(*iter);
      if (atomicOrder >= 0) {
        size_t location = 3 * size_t(atomicOrder);
        for (size_t i(location), j(0); j < 3; i++, j++) {
          mode[i] = massWeight*b[j];
        }
      }
      else {
        string errorMsg = "Problem with calculation of internal rotation eigenvector. Atomic order is not correactly defined.";
        throw (std::runtime_error(errorMsg));
      }
    }

    return true;
  }

  // For a bond between atom at1 and atom at2 find all the atoms connected to at1
  // excluding those connect via at2.
  void gStructure::findRotorConnectedAtoms(vector<string> &atomset, const string at1, const string at2) {
    atomset.push_back(at2); // Will not look beyond this atom on the other side of the bond.
    GetAttachedAtoms(atomset, at1);
    atomset.erase(atomset.begin()); // The other side of the bond is not in this set.
  }

  // Calculate the reduce moment of inertia about axis defined by specifed atoms.
  double gStructure::reducedMomentInertia(pair<string, string>& bondats) {

    vector3 coords1 = GetAtomCoords(bondats.first);
    vector3 coords2 = GetAtomCoords(bondats.second);

    // Calculate moment of inertia about bond axis of atoms on one side of bond...
    vector<string> atomset;
    findRotorConnectedAtoms(atomset, bondats.first, bondats.second);
    double mm1 = CalcMomentAboutAxis(atomset, coords1, coords2);

    //...and the other side of the bond
    atomset.clear();
    findRotorConnectedAtoms(atomset, bondats.second, bondats.first);
    double mm2 = CalcMomentAboutAxis(atomset, coords1, coords2);

    /*
    Is the reduced moment of inertia needed about the bond axis or, separately for the set of
    atoms on each side of the bond, about a parallel axis through their centre of mass?
    See:
    http://www.ccl.net/chemistry/resources/messages/2001/03/21.005-dir/index.html
    http://www.ccl.net/chemistry/resources/messages/2001/03/31.002-dir/index.html
    The bond axis is used here.
    */

    return mm1 * mm2 / (mm1 + mm2); //units a.u.*Angstrom*Angstrom

  }

  // Calculate the reduce moment of inertia about axis defined by specifed atoms.
  void gStructure::reducedMomentInertiaAngular(string bondID, double phase, vector<double>& angles,
    vector<double>& redInvMOI, PersistPtr ppConfigData) {

      // Save coordinates.
      exportToXYZ("orig_coords", NULL, ppConfigData) ;

      map<string, atom>::iterator iter;
      vector<vector3> coordinates;
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        coordinates.push_back(iter->second.coords);
      }

      pair<string, string> bondats = GetAtomsOfBond(bondID);
      atom &at1 = Atoms[bondats.first];
      atom &at2 = Atoms[bondats.second];

      // Move fragment so that at1 is at origin.

      vector3 origin = at1.coords;
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        iter->second.coords -= origin;
      }

      exportToXYZ("origin_at_at1") ;

      // Rotate molecule so that the at1-at2 bond is along z -axis.

      vector3 bond = at2.coords;
      dMatrix rotY(3, 0.0);
      double cosTheta = bond.z() / bond.length();
      double sinTheta = sqrt(1.0 - min(1.0, cosTheta*cosTheta));
      rotY[0][0] = rotY[2][2] = cosTheta;
      rotY[1][1] = 1.0;
      rotY[0][2] = sinTheta;
      rotY[2][0] = -rotY[0][2];

      dMatrix rotZ(3, 0.0);
      double radius = sqrt(bond.x()*bond.x() + bond.y()*bond.y());
      double cosPhi(0.0);
      double sinPhi(0.0);
      if (radius > 0.0) {
        cosPhi = bond.x() / radius;
        sinPhi = bond.y() / radius;
      }
      rotZ[0][0] = rotZ[1][1] = cosPhi;
      rotZ[2][2] = 1.0;
      rotZ[0][1] = sinPhi;
      rotZ[1][0] = -rotZ[0][1];

      dMatrix rotA = rotY*rotZ;

      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        vector<double> r(3, 0.0);
        iter->second.coords.Get(&r[0]);
        r *= rotA;
        iter->second.coords.Set(&r[0]);
      }

      exportToXYZ("z_axis=at1-at2") ;

      // Determine the content of one of the fragments so that it can moved relative to the other.

      vector<string> atomset;
      findRotorConnectedAtoms(atomset, bondats.first, bondats.second);

      // Rotate fragments so that they are in phase with the potential.

      const double Angle(fmod(phase, 360.0)*M_PI / 180.);
      dMatrix rot(3, 0.0);
      rot[0][0] = rot[1][1] = cos(Angle);
      rot[2][2] = 1.0;
      rot[0][1] = sin(Angle);
      rot[1][0] = -rot[0][1];

      for (size_t j(0); j < atomset.size(); j++) {
        atom &at = Atoms[atomset[j]];
        vector<double> r(3, 0.0);
        at.coords.Get(&r[0]);
        r *= rot;
        at.coords.Set(&r[0]);
      }

      // Rotate one fragment relative to the other.

      redInvMOI.clear();
      redInvMOI.push_back(conMntInt2RotCnt*getGRIT(bondID));
      const size_t nAngle(angles.size());
      const double dAngle = 2.0*M_PI / double(nAngle);
      angles[0] = 0.0;
      rot.reset(rot.size());
      rot[0][0] = rot[1][1] = cos(dAngle);
      rot[2][2] = 1.0;
      rot[0][1] = sin(dAngle);
      rot[1][0] = -rot[0][1];

      for (size_t i(1); i < nAngle; i++) {

        angles[i] = angles[i - 1] + dAngle;

        for (size_t j(0); j < atomset.size(); j++) {
          atom &at = Atoms[atomset[j]];
          vector<double> r(3, 0.0);
          at.coords.Get(&r[0]);
          r *= rot;
          at.coords.Set(&r[0]);
        }

        exportToXYZ(("rotate_fragment_" + toString(i)).c_str()) ;
        redInvMOI.push_back(conMntInt2RotCnt*getGRIT(bondID)) ;
      }

      // Restore original coordinates.

      size_t i(0);
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter, i++) {
        iter->second.coords = coordinates[i];
      }

      exportToXYZ("orig_coords", true) ;
  }

  // Calculate the internal rotation eigenvector. Based on the internal rotation 
  // mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010).
  // Typically this vector is used to project out an internal rotational mode
  // from a Hessian.
  void gStructure::internalRotationVector(string bondID, vector<double>& mode, bool ApplyMWeight) {

    pair<string, string> bondats = GetAtomsOfBond(bondID);

    vector3 coords1 = GetAtomCoords(bondats.first);
    vector3 coords2 = GetAtomCoords(bondats.second);

    // Calculate moment of inertia about bond axis of atoms on one side of bond...
    vector<string> atomset1;
    findRotorConnectedAtoms(atomset1, bondats.first, bondats.second);
    double mm1 = CalcMomentAboutAxis(atomset1, coords1, coords2);
    CalcInternalRotVec(atomset1, coords1, coords2, mode, ApplyMWeight);

    //...and the other side of the bond
    vector<string> atomset2;
    findRotorConnectedAtoms(atomset2, bondats.second, bondats.first);
    double mm2 = CalcMomentAboutAxis(atomset2, coords1, coords2);
    CalcInternalRotVec(atomset2, coords2, coords1, mode, ApplyMWeight);

    // In the following weights are applied to relative rotation of each 
    // fragment with respect to each other. In the limit that one fragment
    // is infinitely massive the rotation will be confined to the other 
    // fragment. 

    if (!ApplyMWeight) {
      double fctr1 = mm2 / (mm1 + mm2);
      ApplyInertiaWeighting(atomset1, mode, fctr1);

      double fctr2 = mm1 / (mm1 + mm2);
      ApplyInertiaWeighting(atomset2, mode, fctr2);
    }

  }

  // Calculates the GRIT for the current set of coodinates.
  double gStructure::getGRIT(string bondID) {

    size_t msize(m_RotBondIDs.size() + 3);
    dMatrix GRIT(msize, 0.0);

    // Save coordinates.

    map<string, atom>::iterator iter;
    vector<vector3> coordinates;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      coordinates.push_back(iter->second.coords);
    }

    // Determine centre of mass.

    vector3 centreOfMass;
    double sm(0.0);
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      double mass = atomMass(iter->second.element);
      centreOfMass += iter->second.coords * mass;
      sm += mass;
    }
    centreOfMass /= sm;

    // Determine ther inertial tensor for overall rotation.

    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      vector3 c = iter->second.coords - centreOfMass;
      iter->second.coords = c;
      double  m = atomMass(iter->second.element);
      sxx += m * c.x() * c.x();
      syy += m * c.y() * c.y();
      szz += m * c.z() * c.z();
      sxy += m * c.x() * c.y();
      sxz += m * c.x() * c.z();
      syz += m * c.y() * c.z();
    }

    GRIT[0][0] = syy + szz;
    GRIT[1][1] = sxx + szz;
    GRIT[2][2] = sxx + syy;
    GRIT[0][1] = GRIT[1][0] = -sxy;
    GRIT[0][2] = GRIT[2][0] = -sxz;
    GRIT[1][2] = GRIT[2][1] = -syz;

    // Calculate the velocity vectors (based on the Sharma, Raman and Green vector).

    vector<vector<double> > velocities(m_RotBondIDs.size(), vector<double>(3 * NumAtoms(), 0.0));
    for (size_t i(0); i < m_RotBondIDs.size(); i++) {

      vector<double> velocity(3 * NumAtoms(), 0.0);

      internalRotationVector(m_RotBondIDs[i], velocity, false);

      // Remove centre of mass velocity.
      vector3 centreOfMassVelocity;
      for (size_t j(0); j < m_atomicOrder.size(); j++){
        double mass = atomMass((Atoms.find(m_atomicOrder[j]))->second.element);
        vector3 vtmp;
        vtmp.Set(&velocity[3 * j]);
        centreOfMassVelocity += mass*vtmp;
      }
      centreOfMassVelocity /= sm;

      for (size_t j(0), idx(0); j < m_atomicOrder.size(); j++){
        for (size_t n(0); n < 3; idx++, n++) {
          velocity[idx] -= centreOfMassVelocity[n];
        }
      }

      velocities[i] = velocity;
    }

    // Calculate the internal kinetic energy terms.

    for (size_t i(0), ii(3); i < m_RotBondIDs.size(); i++, ii++) {
      vector<double> &vi = velocities[i];
      for (size_t j(i), jj(ii); j < m_RotBondIDs.size(); j++, jj++) {
        vector<double> &vj = velocities[j];
        double smk(0.0);
        for (size_t m(0), idx(0); m < m_atomicOrder.size(); m++){
          double mass = atomMass((Atoms.find(m_atomicOrder[m]))->second.element);
          for (size_t l(0); l < 3; l++, idx++) {
            smk += mass*vi[idx] * vj[idx];
          }
        }
        GRIT[ii][jj] = smk;
        if (ii != jj)
          GRIT[jj][ii] = GRIT[ii][jj];
      }
    }

    // Calculate the Coriolis terms.

    for (size_t i(0), ii(3); i < m_RotBondIDs.size(); i++, ii++) {
      vector<double> &vi = velocities[i];
      vector3 coriolis;
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        size_t ll = 3 * getAtomicOrder(iter->first);
        vector3 r = iter->second.coords;
        vector3 vtmp;
        vtmp.Set(&vi[ll]);
        double mass = atomMass(iter->second.element);
        coriolis += mass*cross(r, vtmp);
      }
      GRIT[0][ii] = GRIT[ii][0] = coriolis.x();
      GRIT[1][ii] = GRIT[ii][1] = coriolis.y();
      GRIT[2][ii] = GRIT[ii][2] = coriolis.z();
    }

    // Restore original coordinates.

    size_t i(0);
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter, i++) {
      iter->second.coords = coordinates[i];
    }

    if (m_verbose) {
      string MatrixTitle("Generalized rotation inertia tensor:");
      GRIT.print(MatrixTitle, ctest);
      ctest << endl;
    }

    dMatrix invGRIT(GRIT);
    invGRIT.invertLUdecomposition();
    if (m_verbose) {
      string MatrixTitle = "Inverse of Generalized rotation inertia tensor:";
      invGRIT.print(MatrixTitle, ctest);
    }
    // Locate the reduced moment of inertia associated with the bond.

    size_t idx(3);
    for (size_t i(0); i < m_RotBondIDs.size() && m_RotBondIDs[i] != bondID; i++, idx++);

    return invGRIT[idx][idx];
  }

  // Apply inertia weighting to the raw internal rotation velocity vector.
  void gStructure::ApplyInertiaWeighting(vector<string> &atomset, vector<double> &velocity, double fctr) const
  {
    for (size_t j(0); j < atomset.size(); j++)
    {
      size_t i = 3 * size_t(getAtomicOrder(atomset[j]));
      for (size_t n(0); n < 3; i++, n++) {
        velocity[i] *= fctr;
      }
    }
  }

  // Returns the rotational constants (in cm-1) in a vector.
  vector<double> gStructure::CalcRotConsts()
  {
    vector<double> RotConsts(m_PrincipalMI.size(), 0.0); //cm-1
    if (NumAtoms() < 2)
      return RotConsts; //empty
    for (size_t i = 0; i < m_PrincipalMI.size(); ++i) {
      RotConsts[i] = (m_PrincipalMI[i] == 0.0) ? 0.0 : conMntInt2RotCnt / m_PrincipalMI[i];
    }

    return RotConsts;
  }

  double gStructure::CalcSumEMinusHf0(bool UsingAtomBasedThermo, bool useHf298)
  {
    //calculate for each atom (ab initio E - Hf0) and return sum
    if (!ReadStructure())
    {
      cerr << "To use me::Hf0 or Hf298 the molecule needs chemical structure (an atomList at least)" << endl;
      return false;
    }
    double sum = 0.0;
    map<string, double> atomdiffs; //el symbol, diff
    map<string, atom>::iterator iter;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter)
    {
      string el = iter->second.element;
      if (atomdiffs.find(el) == atomdiffs.end())
      {
        //get vals from librarymols.xml
        PersistPtr ppMol = GetFromLibrary(el, PersistPtr());
        double diff(0.0);
        if (ppMol)
        {
          diff = ppMol->XmlReadPropertyDouble("me:ZPE", optional);
          if (useHf298)
          {
            diff -= ppMol->XmlReadPropertyDouble("me:Hf298", optional);
            diff += ppMol->XmlReadPropertyDouble("me:H0-H298", optional);
          }
          else if (!UsingAtomBasedThermo)
            diff -= ppMol->XmlReadPropertyDouble("me:Hf0", optional);
        }
        if (!ppMol || IsNan(diff))
        {
          cerr << "The value of Hf0 for " << getHost()->getName()
            << " will be incorrect because one or more of its elements"
            << " was not in the library, or lacked me:ZPE and me:Hf0 or me:Hf298 properties" << endl;
          return 0.0;
        }
        atomdiffs[el] = diff; //save diff for this el in atomdiffs
        sum += diff;
      }
      else // diff for this el already known
        sum += atomdiffs[el];
    }
    return sum;
  }

  // Returns an ordered array of masses.
  void gStructure::getAtomicMasses(vector<double> &AtomicMasses) const {
    AtomicMasses.clear();
    for (size_t i(0); i < m_atomicOrder.size(); i++){
      double mass = atomMass((Atoms.find(m_atomicOrder[i]))->second.element);
      AtomicMasses.push_back(mass);
    }
  }

  // Returns an ordered array of coordinates.
  void gStructure::getAtomicCoords(vector<double> &coords, AxisLabel cartLabel) const {
    coords.clear();
    for (size_t i(0); i < m_atomicOrder.size(); i++){
      double coord = (Atoms.find(m_atomicOrder[i]))->second.coords[cartLabel];
      coords.push_back(coord);
    }
  }

  // Returns an ordered array of X coordinates.
  void gStructure::getXCoords(vector<double> &coords) const {
    getAtomicCoords(coords, X);
  }

  // Returns an ordered array of Y coordinates.
  void gStructure::getYCoords(vector<double> &coords) const {
    getAtomicCoords(coords, Y);
  }

  // Returns an ordered array of Z coordinates.
  void gStructure::getZCoords(vector<double> &coords) const {
    getAtomicCoords(coords, Z);
  }

  // Export to xmol format.
  void gStructure::exportToXYZ(const char* txt, bool last, PersistPtr ppConfigData){

    // Only write something if the verbosity flag has been set.
    if (!m_verbose)
      return;

    cinfo << Atoms.size() << endl ;
    cinfo << getHost()->getName();
    if(txt)
      cinfo << "-" << txt;
    cinfo << endl ;
    map<string, atom>::const_iterator iter;
    for (iter=Atoms.begin(); iter!=Atoms.end(); ++iter) {
      const atom &at = iter->second ;
      cinfo << setw(10) << at.element ;
      cinfo << formatFloat(at.coords.x(), 6, 15)  
        << formatFloat(at.coords.y(), 6, 15) 
        << formatFloat(at.coords.z(), 6, 15) << endl ;
    }
    exportToCML(txt, last, ppConfigData);
  }

  void gStructure::exportToCML(const char* txt, bool last, PersistPtr ppConfigData){
    static PersistPtr pp;
    if(ppConfigData)
      pp = ppConfigData;
    PersistPtr ppMol = pp->XmlWriteElement("molecule");
    string id = getHost()->getName();
    ppMol->XmlWriteAttribute("id", id  + '-' + txt);
    PersistPtr ppAtList = ppMol->XmlWriteElement("atomArray");
    map<string, atom>::const_iterator iter;
    for (iter=Atoms.begin(); iter!=Atoms.end(); ++iter) {
      const atom &at = iter->second ;
      PersistPtr ppAt = ppAtList->XmlWriteElement("atom");
      ppAt->XmlWriteAttribute("id", at.id);
      ppAt->XmlWriteAttribute("elementType", at.element);
      ppAt->XmlWriteAttribute("x3", at.coords.x(), 4, true);
      ppAt->XmlWriteAttribute("y3", at.coords.y(), 4, true);
      ppAt->XmlWriteAttribute("z3", at.coords.z(), 4, true);
    }
    PersistPtr ppBList = ppMol->XmlWriteElement("bondArray");
    for (map<string, pair<string, string> >::const_iterator iter=Bonds.begin();
      iter!=Bonds.end(); ++iter) {
        PersistPtr ppB = ppBList->XmlWriteElement("bond");
        ppB->XmlWriteAttribute("id", iter->first);
        ppB->XmlWriteAttribute("atomRefs2", iter->second.first + " " + iter->second.second);
        ppB->XmlWriteAttribute("order", "1");
        if(iter->first==m_RotBondIDs.back())
          ppB->XmlWriteAttribute("color", "pink");
    }
  }
}//namespace

//XML written as text; not used.
//void gStructure::exportToCML(const char* txt, bool last, PersistPtr pp){
//  static stringstream ss;
//  static PersistPtr ppConfigData;
//  if(pp)
//    ppConfigData = pp;
//  ss << " <molecule id=\"" << getHost()->getName();
//  if(txt)
//    ss <<'-' << txt;
//  ss << "\">\n" << "  <atomArray>\n";
//  map<string, atom>::const_iterator iter;
//  for (iter=Atoms.begin(); iter!=Atoms.end(); ++iter) {
//    const atom &at = iter->second ;
//    ss  << "   <atom id=" << at.id << " elementType=\"" << at.element << "\" ";
//    ss << setprecision(4) << fixed
//          << "x3=\"" << at.coords.x() <<"\" "
//          << "y3=\"" << at.coords.y() <<"\" "
//          << "z3=\"" << at.coords.z() <<"\"/>\n";
//  }
//  ss << "   </atomArray>\n   <bondArray>";
//  int i(1);
//  for (map<string, pair<string, string> >::const_iterator iter=Bonds.begin();
//         iter!=Bonds.end(); ++iter, ++i) {
//    if(i%2) ss << '\n';
//    ss << "   <bond atomRefs2=\"" << iter->second.first
//          <<' ' << iter->second.second << " order=\"1\"/>";
//  }
//  ss << "\n  </bondArray>\n </molecule>" << endl;
//  if(last)
//    // Written as text (rather than XML) because it is available.
//    // Should be seen ok as valid XML.
//    ppConfigData->XmlWrite(ss.str());
//  ss.str()="";
//}
