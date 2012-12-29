// MolecularComponents.cpp
//
// Author: Chi-Hsiu Liang
//
//-------------------------------------------------------------------------------------------
#include <stdexcept>
#include <numeric>
#include <cmath>
#include "Molecule.h"
#include "System.h"

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
    */};

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
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    setSigma(ppPropList->XmlReadPropertyDouble("me:sigma"));
    setEpsilon(ppPropList->XmlReadPropertyDouble("me:epsilon"));
  }

  void   gBathProperties::setSigma(double value)          {
    m_Sigma = value;
    m_Sigma_chk = 0;
  } ;

  double gBathProperties::getSigma()                      {
    if (m_Sigma_chk >= 0){
      ++m_Sigma_chk;
      return m_Sigma ;
    }
    else{
      cerr << "m_Sigma was not defined but requested in " << m_host->getName() << ". Default value " << sigmaDefault << " is used.\n";
      return m_Sigma ;
    }
  } ;

  void   gBathProperties::setEpsilon(double value)        {
    m_Epsilon = value;
    m_Epsilon_chk = 0;
  } ;

  double gBathProperties::getEpsilon()                    {
    if (m_Epsilon_chk >= 0){
      ++m_Epsilon_chk;
      return m_Epsilon ;
    }
    else{
      cerr << "m_Epsilon was not defined but requested in " << m_host->getName() << ". Default value " << epsilonDefault << " is used.\n";
      return m_Epsilon ;
    }
  } ;

  //-------------------------------------------------------------------------------------------------
  // Cell density of states related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gDensityOfStates::~gDensityOfStates()
  {
    //Delete the density of state calculators because they are cloned instances
    for(unsigned i=0; i<m_DOSCalculators.size(); ++i)
      delete m_DOSCalculators[i];
  }

  gDensityOfStates::gDensityOfStates(Molecule* pMol)
    :m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(1.0),
    m_ZPE(0.0),
    m_scaleFactor(1.0),
    m_SpinMultiplicity(1),
    m_RC_chk(-1),
    m_Sym_chk(-1),
    m_ZPE_chk(-1),
    m_scaleFactor_chk(-1),
    m_SpinMultiplicity_chk(-1),
    m_VibFreq_chk(-1),
    m_EnergyConvention("arbitary"),
    m_eleExc(),
    m_VibFreq(),
    m_grainEne(),
    m_grainDOS() { m_host = pMol; }

  bool gDensityOfStates::initialization() {

    Molecule* pMol = m_host ;

    ErrorContext c(pMol->getName());

    PersistPtr pp = pMol->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //a propertyList element is not essential

    const char* txt;

    bool hasVibFreq(true), hasRotConst(true) ;
    txt= ppPropList->XmlReadProperty("me:vibFreqs", optional);
    if (!txt) {
      hasVibFreq = false;
      if(!pMol->getStruc().IsAtom())
        cinfo << "Cannot find argument me:vibFreqs. Assumes that it is an atom or atomic ion." << endl;
      m_VibFreq_chk = -1;
    } else { 
      istringstream idata(txt);
      double x; 
      while (idata >> x)
        m_VibFreq.push_back(x); m_VibFreq_chk = 0;
    }

    std::vector<double> rCnst(3, 0.0);
    hasRotConst = false;
    txt = ppPropList->XmlReadProperty("me:rotConsts", optional);
    if (txt) {
      //data from <me:rotConsts>
      istringstream idata(txt);
      idata >> rCnst[0] >> rCnst[1] >> rCnst[2];
      txt = ppPropList->XmlReadPropertyAttribute("me:rotConsts", "units", optional) ;
      if (string(txt) == "amuA^2") {
        // The supplied data are moments of inertia with units amuAng^2.
        rCnst[0] = conMntInt2RotCnt/rCnst[0] ;
        rCnst[1] = conMntInt2RotCnt/rCnst[1] ;
        rCnst[2] = conMntInt2RotCnt/rCnst[2] ;
      }
      hasRotConst = true;
      m_RC_chk = 0;
    } else {
      ////data from atom coordinates
      gStructure& gs = pMol->getStruc();
      if(gs.ReadStructure()) {
        rCnst = gs.CalcRotConsts();
        cinfo << "Rotational constants were calculated from atom coordinates: "
          << m_RotCstA << ' ' << m_RotCstB << ' ' << m_RotCstC << " cm-1" << endl;
        hasRotConst = true; 
      }
    }
    if (hasRotConst) {
      // Check rotational constants are valid.
      if (rCnst[0] < 0.0 || rCnst[1] < 0.0 || rCnst[2] < 0.0) {
        cerr << "A negative rotational constant found "<< endl ;
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
        "Assumed to be an atom or atomic ion." << endl;
    }
    else
      m_RC_chk = 0;

    if (hasVibFreq != hasRotConst){
      cerr << "Improper setting on vibrational frequencies or rotational constants." << endl;
    }

    txt = ppPropList->XmlReadProperty("me:electronicExcitation", optional);
    if (txt) {
      istringstream idata(txt);
      m_eleExc.clear();
      double _iele;
      while (idata >> _iele) m_eleExc.push_back(_iele);
    }

    m_Sym = ppPropList->XmlReadPropertyDouble("me:symmetryNumber");
    m_Sym_chk = 0;

    m_scaleFactor = ppPropList->XmlReadPropertyDouble("me:frequenciesScaleFactor");
    m_scaleFactor_chk = 0;

    // Read attribute on the property and then, if not present, on the attribute (with a a default).
    m_SpinMultiplicity = ppPropList->XmlReadPropertyInteger("me:spinMultiplicity", optional);
    if(m_SpinMultiplicity==0)
      m_SpinMultiplicity = pp->XmlReadInteger("spinMultiplicity");
    m_SpinMultiplicity_chk = 0;

    ReadDOSMethods();

    return ReadZeroPointEnergy(ppPropList) ;
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
    if(!IsNan(tempzpe))
    {
      // me:ZPE is present
      unitsInput = ppPropList->XmlReadPropertyAttribute("me:ZPE", "units");
      txt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "convention", optional);
      m_EnergyConvention = txt ? txt : "arbitary";

      double zpCorrection=0.0; //cm-1
      txt= ppPropList->XmlReadPropertyAttribute("me:ZPE", "zeroPointVibEnergyAdded", optional);
      if(txt && strcmp(txt, "false")==0)
      {
        // Mesmer assumes that the ZPE is the true zero point energy, unless
        // an attribute zeroPointVibEnergyAdded="false" is present. This indicates
        // that the value provided is a computed energy at the bottom of the well
        // and it is corrected here by adding 0.5*Sum(vib freqs).
        if(m_VibFreq.size() > 0)
        {
          zpCorrection = accumulate(m_VibFreq.begin(),m_VibFreq.end(), 0.0);
          zpCorrection *= 0.5;
        }
        //Write back a corrected value and change attribute to zeroPointVibEnergyAdded="true"
        PersistPtr ppScalar = ppPropList->XmlMoveToProperty("me:ZPE");
        ppScalar->XmlWrite(toString(tempzpe + ConvertFromWavenumbers(unitsInput, zpCorrection)));
        ppScalar->XmlWriteAttribute("zeroPointVibEnergyAdded", "true");
      }

      bool rangeSet ;
      PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:ZPE"); 
      ReadRdoubleRange(string(m_host->getName()+":ZPE"), ppProp, m_ZPE, rangeSet, getConvertedEnergy(unitsInput, 1.0), zpCorrection) ;
      set_zpe(getConvertedEnergy(unitsInput, tempzpe) + zpCorrection);
      m_ZPE_chk = 0;

    } else {

      // me:ZPE not present; try me:Hf0 and HAT0 (enthalpy of atomization at 0K)
      // and Hf298 (enthalpy of formation at 298K)
      double Hf0 = ppPropList->XmlReadPropertyDouble("me:Hf0", optional);
      double HfAT0 = ppPropList->XmlReadPropertyDouble("me:HfAT0", optional);
      double Hf298 = ppPropList->XmlReadPropertyDouble("me:Hf298", optional);
      if(IsNan(Hf0) && IsNan(HfAT0) && IsNan(Hf298))
      {
        //None of me:ZPE, me:Hf0, me:HAT0 and meHf298 are present; use a default ZPE.
        tempzpe = ppPropList->XmlReadPropertyDouble("me:ZPE", true);
      }
      else
      {
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
        const char* utxt= ppPropList->XmlReadPropertyAttribute(origElement, "units", optional);
        utxt = utxt ? utxt : "kJ/mol";
        if(!IsNan(Hf0))
        {
          //Use Hf0 if provided
          Hf0 = getConvertedEnergy(utxt, Hf0); //cm-1
          tempzpe = Hf0 + getConvertedEnergy("kJ/mol", getHost()->getStruc().CalcSumEMinusHf0(false, false));//cm-1
        }
        else if(!IsNan(HfAT0))
        {
          //Use HfAT0 (atom-based thermochemistry, see DOI: 10.1002/chem.200903252)
          HfAT0 = getConvertedEnergy(utxt, HfAT0); //cm-1
          tempzpe = HfAT0 + getConvertedEnergy("kJ/mol", getHost()->getStruc().CalcSumEMinusHf0(true, false));//cm-1
        }
        else
        {
          //Use Hf298
          const char* convention = ppPropList->XmlReadPropertyAttribute(origElement, "convention", optional);
          if(convention && strcmp(convention,"thermodynamic")==0)
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
      }
    }
    return true ;
  }

  //
  // Get the number of degrees of freedom for this species.
  //
  unsigned int gDensityOfStates::getNoOfDegOfFreeedom() {
    unsigned int nDOF(0) ;
    for ( vector<DensityOfStatesCalculator*>::size_type j = 0 ; j < m_DOSCalculators.size() ; ++j ) {
      nDOF += m_DOSCalculators[j]->NoDegOfFreedom(this) ;
    }
    return nDOF ;
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
  /* Rotational-electronic densities of states objects are specified in XML file as:
     and/or   <me:DOSCMethod name=" name "/>                  preferred
	   multiple <me:DOSCMethod> name </me:DOSCMethod>
     and/or   <me:ExtraDOSCMethod> name </me:ExtraDOSCMethod> deprecated
     and/or   <me:ExtraDOSCMethod name=" name "/>             deprecated
  */
    ErrorContext c(getHost()->getName());

    PersistPtr pp = getHost()->get_PersistentPointer();
    if(!ReadMethodsFromXml("me:DOSCMethod"))
      return false;
    if(!ReadMethodsFromXml("me:ExtraDOSCMethod"))
      return false;
    if(m_DOSCalculators.empty()) //invoke default method
    {
      const char* name = pp->XmlReadValue("me:DOSCMethod");
      if(!name)
        return false;
      DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find(string(name));
      if(!pDOSCalculator)
        return false;
      m_DOSCalculators.push_back(pDOSCalculator);
    }
	//
	// Beyer-Swinehart object added by default.
	//
    DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find("BeyerSwinehart");
    m_DOSCalculators.push_back(pDOSCalculator);
    if(!pDOSCalculator)
    {
      cerr << "Beyer-Swinhart algorithm failed to initialize correctly";
      return false;
    }
    return true;
  }

  bool gDensityOfStates::ReadMethodsFromXml(const string& keyword) {
    PersistPtr pp = getHost()->get_PersistentPointer();
    while(pp = pp->XmlMoveTo(keyword))
    {
      string dosMethod;
      const char* name;
      if( (name = pp->XmlRead()) || (name = pp->XmlReadValue("name", optional)))
        dosMethod = name;
      DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find(dosMethod);
      m_DOSCalculators.push_back(pDOSCalculator);
      if(!pDOSCalculator->ReadParameters(this, pp))
      {
        cerr << dosMethod << " failed to initialize correctly";
        return false;
      }
    }
    return true;
  }

  //
  // Get Electronic excitations
  //
  void gDensityOfStates::getEleExcitation(vector<double> &elecExci){
    elecExci.clear();
    for (vector<double>::size_type i = 0; i < m_eleExc.size(); ++i){
      elecExci.push_back(m_eleExc[i]);
    }
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
        cinfo << "Rotational constants were not defined but requested." << endl;
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

    if      ((mmtsInt[0] + mmtsInt[1] + mmtsInt[2]) == 0.) 
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
    const bool recalc(needReCalculateDOS());
    const int MaximumCell = m_host->getEnv().MaxCell;
    const bool vectorSizeConstant(m_cellDOS.size() == static_cast<unsigned int>(MaximumCell));
    const size_t sizeOfVector(m_cellDOS.size());

    if (sizeOfVector && vectorSizeConstant && !recalc)
      return true;

    // Calculate density of states.
    bool ret(true) ; 
    for(unsigned i=0; ret && i<m_DOSCalculators.size(); ++i)
      ret = ret && m_DOSCalculators[i]->countCellDOS(this, MaximumCell);
    if(!ret)
      return false;
    //-------------------------------------------------------------

    std::vector<double> shiftedCellDOS;
    std::vector<double> shiftedCellEne;
    const int cellOffset = get_cellOffset();
    std::vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);
    shiftCells(MaximumCell, cellOffset, m_cellDOS, cellEne, shiftedCellDOS, shiftedCellEne);

    calcGrainAverages(m_host->getEnv().MaxGrn, m_host->getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, m_grainDOS, m_grainEne);

    testDensityOfStates();

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
    const int MaximumCell  = m_host->getEnv().MaxCell;
    std::vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);

    // Partition functions that are higher than the current simulation temperature will not be output.
    const double temperature = 1. / (boltzmann_RCpK * m_host->getEnv().beta);
    const int max_nplus1 = int(temperature / 100.);

    if (m_host->isMolType("modelled") || m_host->isMolType("transitionState")){
      string comment("Rovibronic partition function calculation at various temperatures. qtot : product of QM partition functions for vibrations (1-D harmonic oscillator) and classical partition functions for rotations.  sumc : cell based partition function. sumg : grain based partition function ");

      PersistPtr ppList = m_host->get_PersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment );

      if (m_host->getFlags().testDOSEnabled) ctest << endl << "Test rovibronic density of states for: " << m_host->getName() << "\n{\n";
      if (m_host->getFlags().testDOSEnabled) ctest << "      T           qtot           sumc           sumg\n";


      //loop through predefined test temperatures
      for ( int n = 0 ; n < max_nplus1 ; ++n ) {
        double temp = 100.0*static_cast<double>(n + 2) ;
        double beta = 1.0/(boltzmann_RCpK*temp) ;

        // Calculate rovibronic partition functions based on cells.
        double cellCanPrtnFn = canonicalPartitionFunction(m_cellDOS, cellEne, beta) ;

        // Calculate rovibronic partition functions based on grains.
        double grainCanPrtnFn = canonicalPartitionFunction(m_grainDOS, m_grainEne, beta) ;

        // Calculate rovibronic partition functions, using analytical formula where possible.
        double qtot(1.0) ;
        for ( vector<DensityOfStatesCalculator*>::size_type j = 0 ; j < m_DOSCalculators.size() ; ++j ) {
          qtot *= m_DOSCalculators[j]->canPrtnFnCntrb(this, beta) ;
        }        

        if (m_host->getFlags().testDOSEnabled) { 
          formatFloat(ctest, temp,  6,  7) ;
          formatFloat(ctest, qtot,  6, 15) ;
          formatFloat(ctest, cellCanPrtnFn,  6, 15) ;
          formatFloat(ctest, grainCanPrtnFn,  6, 15) ;
          ctest << endl ;
        }

        //Add to XML document
        PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
        ppItem->XmlWriteValueElement("me:T",    temp, 6);
        ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
        ppItem->XmlWriteValueElement("me:sumc", cellCanPrtnFn, 6);
        ppItem->XmlWriteValueElement("me:sumg", grainCanPrtnFn, 6);
      }
      if (m_host->getFlags().testDOSEnabled) ctest << "}" << endl;
    }

    if (m_host->getFlags().cellDOSEnabled){
      ctest << endl << "Cell rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i(0); i < MaximumCell; ++i){
        formatFloat(ctest, cellEne[i],  6,  15) ;
        formatFloat(ctest, m_cellDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }

    if (m_host->getFlags().grainDOSEnabled && (m_host->isMolType("modelled") || m_host->isMolType("transitionState"))){
      ctest << endl << "Grain rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i(0); i < MaximumGrain; ++i){
        formatFloat(ctest, m_grainEne[i],  6,  15) ;
        formatFloat(ctest, m_grainDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }
  }

  double gDensityOfStates::get_zpe() {
    if (m_ZPE_chk == -1) {
      cinfo << "m_ZPE was not defined but requested in " << m_host->getName() << ". Default value " << m_ZPE << " is given." << endl;
      --m_ZPE_chk;
    } else if (m_ZPE_chk < -1) {
      --m_ZPE_chk;
    } else {
      ++m_ZPE_chk;
    }
    return double(m_ZPE);
  }

  double gDensityOfStates::get_scaleFactor() {
    if (m_scaleFactor_chk == -1){
      cinfo << "m_scaleFactor was not defined but requested in " << m_host->getName() << ". Default value " << m_scaleFactor << " is given." << endl;
      --m_scaleFactor_chk;
    } else if (m_scaleFactor_chk < -1) {
      --m_scaleFactor_chk;
    } else {
      ++m_scaleFactor_chk;
    }
    return m_scaleFactor ;
  }

  double gDensityOfStates::get_Sym(void){
    if (m_Sym_chk == -1){
      cinfo << "m_Sym was not defined but requested in " << m_host->getName() << ". Default value " << m_Sym << " is given." << endl;
      --m_Sym_chk;
    } else if (m_Sym_chk < -1) {
      --m_Sym_chk;
    } else {
      ++m_Sym_chk;
    }
    return m_Sym ;
  }

  void gDensityOfStates::get_VibFreq(std::vector<double>& vibFreq){
    if (m_VibFreq_chk >=0){
      const double scalefactor = get_scaleFactor();
      for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
        vibFreq.push_back(m_VibFreq[i] * scalefactor);
      ++m_VibFreq_chk;
    }
  }

  bool gDensityOfStates::removeVibFreq(double freq) {
    vector<double>::iterator pos = find(m_VibFreq.begin(),m_VibFreq.end(), freq);
    if(pos==m_VibFreq.end())
      return false;
    m_VibFreq.erase(pos);
    return true;
  }

  int gDensityOfStates::getSpinMultiplicity(){
    if (m_SpinMultiplicity_chk >= 0) {
      ++m_SpinMultiplicity_chk;
    } else {
      cinfo << "m_SpinMultiplicity was not defined but requested in " << m_host->getName() << ". Default value " << m_SpinMultiplicity << " is given." << endl;
    }
    return m_SpinMultiplicity ;
  }

  int gDensityOfStates::get_cellOffset(void) {
    double modulus = fmod(get_zpe() - m_host->getEnv().EMin, m_host->getEnv().GrainSize);
    if(modulus < 0.0)  // presently modulus is only less than 0 for the excess reactant in an association rxn
      modulus = 0.0;   // however, this problem should become obsolete once supermolecule DOS is calculated on the fly
    return int(modulus) ;
  } ;

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
      for (int i(startGrnIdx+1); i < int(m_grainDOS.size()); ++i){
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

    return canonicalPartitionFunction(m_grainDOS, m_grainEne, m_host->getEnv().beta) ;

  }

  //
  // Calculate standard thermodynamic quantities as a function of temperature.
  // The calculation is based on cell densities of states.
  //
  bool gDensityOfStates::thermodynamicsFunctions(double temp, double unitFctr, double& enthalpy, double& entropy, double& gibbsFreeEnergy) {

    std::vector<double> cellEne;
    getCellEnergies(m_host->getEnv().MaxCell, cellEne);

    calcDensityOfStates() ;

    double beta ;
    if ( temp > 0.0 ) { 
      beta = 1.0/(boltzmann_RCpK*temp) ;
    } else {
      return false ;
    }

    // Calculate rovibronic partition functions based on cells.
    double cellCanPrtnFn = canonicalPartitionFunction(m_cellDOS, cellEne, beta) ;

    // The following calculates the mean internal molecular energy.
    double internalEnergy = canonicalMeanEnergy(m_cellDOS, cellEne, beta) ;

    // The rovibronic partition function must be corrected for translation 
    // and (assuming an ideal gas) molecular indistinguishability.
    gibbsFreeEnergy = unitFctr*(-log(cellCanPrtnFn) - log(tp_C * pow((m_host->getStruc().getMass() / beta), 1.5)) + log(AvogadroC))/beta ;

    // The enthalpy must be corrected for translation by an additional 3kT/2.
    enthalpy        = unitFctr*(internalEnergy + 5.0/(2.0*beta)) ;

    entropy         = (enthalpy - gibbsFreeEnergy)/temp  ;

    return true ;
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
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;
    txt= ppPropList->XmlReadProperty("me:imFreqs", optional);
    if (!txt) {
      cinfo << "No imaginary vibrational frequency. Check if there is tunneling.\n";
      m_ImFreq_chk = -1;
      // Note: If there is tunneling, an imaginary frequency must also be supplied.
    } else {
      istringstream idata(txt); double x;
      // make sure this number is positive
      idata >> x ;
      set_imFreq(x);
      bool rangeSet ;
      PersistPtr ppProp = ppPropList->XmlMoveToProperty("me:imFreqs"); 
      ReadRdoubleRange(string(m_host->getName()+":imFreqs"), ppProp, m_ImFreq, rangeSet) ;
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

  //-------------------------------------------------------------------------------------------------
  // Collisional redistribution related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gWellProperties::gWellProperties(Molecule* pMol) : MolecularComponent(),
    m_collisionFrequency(0.0),
    m_ncolloptrsize(0),
    m_lowestBarrier(9e23),
    m_numGroupedGrains(0),
    m_pDistributionCalculator(NULL),
    m_pEnergyTransferModel(NULL),
    m_grainFracBeta(0.),
    m_grainDist(),
    m_egme(NULL),
    m_egvec(NULL),
    m_egval()
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
  }

  gWellProperties::~gWellProperties()
  {
    if (m_egme != NULL) delete m_egme ;
    if (m_grainDist.size()) m_grainDist.clear();
    delete m_pEnergyTransferModel;
  }

  bool gWellProperties::initialization(){

    // Determine the method of DOS calculation.

    PersistPtr pp = m_host->get_PersistentPointer();
    const char* pDistCalcMethodtxt = pp->XmlReadValue("me:DistributionCalcMethod") ;
    if(pDistCalcMethodtxt)
    {
      m_pDistributionCalculator = DistributionCalculator::Find(pDistCalcMethodtxt);
      if(!m_pDistributionCalculator) // if the provided method cannot be found,
      {
        cinfo << "Unknown method " << pDistCalcMethodtxt
          << " for the calculation of distribution fraction. Please check spelling error. Default method <Boltzmann> is used." << endl;
        pDistCalcMethodtxt = "Boltzmann";
        m_pDistributionCalculator = DistributionCalculator::Find(pDistCalcMethodtxt);
      }
    }
    else{ // if no method is provided.
      cinfo << "No method for the calculation of distribution fraction is provided in " << m_host->getName() << ". Default method <Boltzmann> is used." << endl;
      pDistCalcMethodtxt = "Boltzmann"; // must exist
      m_pDistributionCalculator = DistributionCalculator::Find(pDistCalcMethodtxt);
    }

    // Specify the energy transfer probability model.
    // The default value is specified in defaults.xml
    const char* pETPModeltxt = pp->XmlReadValue("me:energyTransferModel") ;
    if(!pETPModeltxt)
      return false;

    m_pEnergyTransferModel = EnergyTransferModel::Find(pETPModeltxt);//new instance, deleted in destructor
    if(!m_pEnergyTransferModel)
    {
      cerr << "Unknown energy transfer model" << pETPModeltxt << " for species" << m_host->getName() << endl;
      return false;
    }

    // Initialize energy transfer model.
    return m_pEnergyTransferModel->ReadParameters(getHost());
  }

  double gWellProperties::get_collisionFrequency() const {
    return m_collisionFrequency ;
  } ;

  void gWellProperties::set_colloptrsize(int ncolloptrsize) {
    m_ncolloptrsize = ncolloptrsize ;
  } ;

  int  gWellProperties::get_colloptrsize() const {
    return m_ncolloptrsize ;
  } ;

  const int gWellProperties::get_grnZPE(){
    double grnZpe = (m_host->getDOS().get_zpe() - m_host->getEnv().EMin) / m_host->getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  //
  // Initialize the Collision Operator.
  //
  bool gWellProperties::initCollisionOperator(double beta, Molecule *pBathGasMolecule)
  {
    // If density of states have not already been calcualted then do so.
    if (!m_host->getDOS().calcDensityOfStates()){
      cerr << "Failed calculating DOS";
      return false;
    }

    // Calculate the collision frequency.
    m_collisionFrequency = collisionFrequency(beta, m_host->getEnv().conc, pBathGasMolecule) ;

    // Calculate the collision operator.
    {
      //-----------------------------------------
      // Treat reservoir grains as a source grain (anything going into this grain will behave Boltzmann)
      // First need to find out the lowest barrier associated with the current well
      vector<double> gEne;
      vector<double> gDOS;
      m_host->getDOS().getGrainEnergies(gEne);
      m_host->getDOS().getGrainDensityOfStates(gDOS);

      // Determine the energy of the reservoir grain
      PersistPtr pp = m_host->get_PersistentPointer();
      PersistPtr ppReservoirSize = pp->XmlMoveTo("me:reservoirSize");

      m_numGroupedGrains = 0; // Reset the number of grains grouped into a reservoir grain to zero.

      while (ppReservoirSize){

        // Check the size of the reservoir.
        const char* pChunkSizeTxt = pp->XmlReadValue("me:reservoirSize") ;
        double tmpvalue(0.0); stringstream s2(pChunkSizeTxt); s2 >> tmpvalue ;

        const char* unitsTxt = ppReservoirSize->XmlReadValue("units", false) ;
        string unitsInput;
        if (unitsTxt){
          unitsInput = unitsTxt;
        } else {
          ctest << "No unit for reservoir size has been supplied, use kJ/mol." << endl;
          unitsInput = "kJ/mol";
        }

        const double value(getConvertedEnergy(unitsInput, tmpvalue));
        int grainLoc(int(value/ double(m_host->getEnv().GrainSize)));
        int lowestBarrier = int(getLowestBarrier() / double(m_host->getEnv().GrainSize));

        if (grainLoc > 0){
          if (grainLoc > lowestBarrier){
            ctest << "The reservoir size provided is too high, corrected according to the lowest barrier height." << endl;
            grainLoc = lowestBarrier;
          }
        } else {
          if (abs(grainLoc) > lowestBarrier){
            ctest << "The reservoir size provided is too low, corrected to zero." << endl;
            grainLoc = 0;
            break;
          } else {
            grainLoc = lowestBarrier + grainLoc;
          }
        }
        ctest << "The reservoir is set to " << grainLoc << " grains, which is about " << grainLoc * m_host->getEnv().GrainSize
          << " cm-1 from the well bottom." << endl;

        // Second find out the partition fraction of active states in the current temperature
        double popAbove(0.0), totalPartition(0.0);
        for (int i(0); i < m_ncolloptrsize; ++i){
          const double ptfn(sqrt(exp(log(gDOS[i]) - beta * gEne[i] + 10.0)));
          totalPartition += ptfn;
          if (i >= grainLoc)
            popAbove += ptfn;
        }

        popAbove /= totalPartition;

        m_numGroupedGrains = grainLoc;
        ctest << popAbove << " of the " << m_host->getName() << " population is in the active states. "
          << "The reservoir size = " << m_numGroupedGrains * m_host->getEnv().GrainSize
          << " cm-1, which is " << m_numGroupedGrains * m_host->getEnv().GrainSize / 83.593 << " kJ/mol." << endl;
        break;
      }

      if (m_numGroupedGrains){
        if (!collisionOperatorWithReservoirState(beta, m_ncolloptrsize - m_numGroupedGrains)){
          cerr << "Failed building collision operator with reservoir state.";
          return false;
        }
      }
      else{
        if (!collisionOperator(beta)){
          cerr << "Failed building collision operator.";
          return false;
        }
      }
    }

    //
    // For the basis set method diagonalize the collision operator to obtain
    // the contracted basis set.
    //
    if (m_host->getEnv().useBasisSetMethod) 
      diagonalizeCollisionOperator();

    return true;
  }

  //
  // Diagonalize collision operator
  //
  void gWellProperties::diagonalizeCollisionOperator()
  {
    // Allocate memory.
    m_egval.clear();
    m_egval.resize(m_ncolloptrsize, 0.0);
    if (m_egvec) delete m_egvec ;                      // Delete the existing matrix.
    m_egvec = new dMatrix(m_ncolloptrsize) ;

    // copy the values over
    for (int i(0); i < m_ncolloptrsize; ++i){
      for (int j(0); j < m_ncolloptrsize; ++j){
        (*m_egvec)[i][j] = (*m_egme)[i][j];
      }
    }

    m_egvec->diagonalize(&m_egval[0]) ;

    bool reportBasisSetDetails(false) ;
    if (reportBasisSetDetails) {
      ctest << "\nEigenvectors of: " << m_host->getName() << endl;
      m_egvec->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);

      // The eigenvalues in this series should contain a zero eigenvalue as 
      // energy transfer operators are conservative.
      ctest << "\nEigenvalues of: " << m_host->getName() << "\n{\n";
      for(int i(0); i < m_ncolloptrsize; ++i){
        ctest << m_egval[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //
  // Calculate collision operator
  //
  bool gWellProperties::collisionOperatorWithReservoirState(double beta, const int reducedCollOptrSize)
  {
    if (m_host->getDOS().test_rotConsts() == UNDEFINED_TOP) return true;
    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    // Initialisation and error checking.
    for (int i(0) ; i < m_ncolloptrsize ; ++i ) {
      if (gDOS[i] <= 0.0) {
        cerr << "Data indicates that grain " << i << " of the current colliding molecule has no states.";
        cerr << "Usually for small GrainSize DOSs calculated by QM rotor method some grain may have no states.";
      }
    }

    // Allocate memory.
    if (m_egme) delete m_egme ;                       // Delete any existing matrix.
    m_egme = new dMatrix(reducedCollOptrSize + 1) ;

    //---------------------------------------------------------------------------------------------------
    //-------------------- The part doing the same jobs as making a whole collision operator ------------
    dMatrix* tempEGME = new dMatrix(m_ncolloptrsize);

    // Calculate raw transition matrix.
    if (!rawTransitionMatrix(beta, gEne, gDOS, tempEGME)) return false ;

    if (m_host->getFlags().showCollisionOperator != 0){
      ctest << "\nCollision operator of " << m_host->getName() << " before normalization:\n";
      tempEGME->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    //Normalisation
    tempEGME->normalizeProbabilityMatrix();

    //-------------------- The part doing the same jobs as making a whole collision operator ------------
    //---------------------------------------------------------------------------------------------------

    if (m_host->getFlags().showCollisionOperator >= 1){
      ctest << "\nCollision operator of " << m_host->getName() << " after normalization:\n";
      tempEGME->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    // print out of column sums to check normalization results
    if (m_host->getFlags().reactionOCSEnabled){
      ctest << endl << "Collision operator column Sums" << endl << "{" << endl ;
      for (int i(0) ; i < m_ncolloptrsize; ++i ) {
        double columnSum(0.0) ;
        for (int j(0) ; j < m_ncolloptrsize; ++j ){
          columnSum += to_double((*tempEGME)[j][i]) ;
        }
        ctest << columnSum << endl ;
      }
      ctest << "}" << endl;
    }

    //--------------------------------
    // COPY, SUMMATION AND SUBSTITUTE
    // Need to copy things over, the active states first.
    for (int i(m_numGroupedGrains) ; i < m_ncolloptrsize ; ++i ) {
      for (int j(m_numGroupedGrains); j < m_ncolloptrsize ; ++j ) {
        (*m_egme)[i - m_numGroupedGrains + 1][j - m_numGroupedGrains + 1] = (*tempEGME)[i][j];
      }
    }

    // Sum up the downward transition terms for the reservoir grain
    double sumOfDeactivation(0.0), ptfReservoir(0.0);
    for (int j(0); j < m_ncolloptrsize ; ++j ) {
      if (j < m_numGroupedGrains){
        // summing up the partition function of reservoir state
        ptfReservoir += exp(log(gDOS[j]) - beta * gEne[j] + 10.0);
      }
      else{
        double downwardSum(0.0);
        for (int i(0) ; i < m_numGroupedGrains ; ++i ) {
          downwardSum += (*tempEGME)[i][j]; // sum of the normalized downward prob.
        }
        double ptfj = exp(log(gDOS[j]) - beta * gEne[j] + 10.0);
        sumOfDeactivation += downwardSum * ptfj;

        (*m_egme)[0][j - m_numGroupedGrains + 1] = downwardSum;
      }
    }
    sumOfDeactivation /= ptfReservoir; // k_a * x_r = k_d(E) * f(E) / Q_a * x_a
    // where Q_a is equal to x_a and cancelled out.
    // So, k_a = k_d(E) * f(E) / x_r;

    (*m_egme)[0][0] = -sumOfDeactivation;

    // Symmetrization of the collision matrix.
    vector<double> popDist; // grained population distribution
    const double firstPop = exp(log(gDOS[0]) - beta * gEne[0] + 10.0);
    popDist.push_back(firstPop);
    for (int idx(1); idx < m_ncolloptrsize; ++idx){
      if (idx < m_numGroupedGrains){
        popDist[0] += exp(log(gDOS[idx]) - beta * gEne[idx] + 10.0);
      }
      else{
        popDist.push_back(sqrt(exp(log(gDOS[idx]) - beta * gEne[idx] + 10.0)));
      }
    }
    popDist[0] = sqrt(popDist[0]); // This is the square root of partition function in the reservoir grain

    for (int i(1) ; i < reducedCollOptrSize + 1; ++i ) {
      for (int j(0) ; j < i ; ++j ){
        (*m_egme)[j][i] *= popDist[i]/popDist[j] ;
        (*m_egme)[i][j]  = (*m_egme)[j][i] ;
      }
    }

    // Account for collisional loss by subrtacting unity from the leading diagonal.
    for (int i(1) ; i <= reducedCollOptrSize ; ++i ) {
      (*m_egme)[i][i] -= 1.0 ;
    }

    if (m_host->getFlags().showCollisionOperator >= 2){
      ctest << "Collision operator of " << m_host->getName() << " after :\n";
      m_egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    delete tempEGME;

    return true;
  }

  double gWellProperties::getBoltzmannWeightedEnergy (int numberOfGrains, const vector<double>& gEne, const vector<double>& gDos, double beta, double& totalDOS)
  {
    double totalEP(0.0), totalP(0.0);
    totalDOS = 0.0;
    for (int i(0); i < numberOfGrains; ++i) {
      totalEP += gEne[i] * exp(log(gDos[i]) - beta * gEne[i] + 10.0);
      totalP += exp(log(gDos[i]) - beta * gEne[i] + 10.0);
      totalDOS += gDos[i];
    }
    return (totalEP / totalP);
  }

  //
  // Calculate collision operator
  //
  bool gWellProperties::collisionOperator(double beta)
  {
    if (m_host->getDOS().test_rotConsts() == UNDEFINED_TOP){   // davidglo, requirement for the diamond work.
      cinfo << "For " << m_host->getName() << ", there no rotational states are defined... only vibrations will be used to construct the molecular DOS" << endl;
    }
    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    // Initialisation and error checking.
    int i, j;
    for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
      if (gDOS[i] <= 0.0) {
        cerr << "Data indicates that grain " << i << " of the current colliding molecule has no states.";
        cerr << "Usually for small GrainSize DOSs calculated by QM rotor method some grain may have no states.";
      }
    }

    // Allocate memory.
    if (m_egme) delete m_egme ;                       // Delete any existing matrix.
    m_egme = new dMatrix(m_ncolloptrsize) ;           // Collision operator matrix.

    // Calculate raw transition matrix.
    if (!rawTransitionMatrix(beta, gEne, gDOS, m_egme)) return false ;

    if (m_host->getFlags().showCollisionOperator != 0){
      ctest << "\nCollision operator of " << m_host->getName() << " before normalization:\n";
      m_egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    //Normalisation
    m_egme->normalizeProbabilityMatrix();

    if (m_host->getFlags().showCollisionOperator >= 1){
      ctest << "\nCollision operator of " << m_host->getName() << " after normalization:\n";
      m_egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    // print out of column sums to check normalization results
    if (m_host->getFlags().reactionOCSEnabled){
      ctest << endl << "Collision operator column sums and energy transfer parameters" << endl << "{" << endl ;
      ctest << " Column Sums           E   <Delta E>  <Delta E>d  <Delta E>u" << endl ;
      for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
        double columnSum(0.0) ;
        double meanEnergyTransfer(0.0) ;
        double meanEnergyTransferDown(0.0) ;
        double meanEnergyTransferUp(0.0) ;
        for ( j = 0 ; j < m_ncolloptrsize ; ++j ){
          columnSum += to_double((*m_egme)[j][i]) ;
          meanEnergyTransfer += (gEne[j] - gEne[i])*to_double((*m_egme)[j][i]) ;
          if (gEne[j] < gEne[i]) {
            meanEnergyTransferDown += (gEne[j] - gEne[i])*to_double((*m_egme)[j][i]) ;
          } else {
            meanEnergyTransferUp += (gEne[j] - gEne[i])*to_double((*m_egme)[j][i]) ;
          }
        }
        ctest << formatFloat(columnSum, 3, 12) 
          << formatFloat(gEne[i], 3, 12)
          << formatFloat(meanEnergyTransfer, 3, 12)
          << formatFloat(meanEnergyTransferDown, 3, 12)
          << formatFloat(meanEnergyTransferUp, 3, 12)
          << endl ;
      }
      ctest << "}" << endl;
    }

    // Symmetrization of the collision matrix.
    vector<double> popDist; // grained population distribution
    for (int idx(0); idx < m_ncolloptrsize; ++idx){
      popDist.push_back(sqrt(exp(log(gDOS[idx]) - beta * gEne[idx] + 10.0)));
    }
    for ( i = 1 ; i < m_ncolloptrsize ; ++i ) {
      for ( j = 0 ; j < i ; ++j ) {
        (*m_egme)[j][i] *= popDist[i]/popDist[j] ;
        (*m_egme)[i][j]  = (*m_egme)[j][i] ;
      }
    }

    //account for collisional loss by subrtacting unity from the leading diagonal.
    for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
      (*m_egme)[i][i] -= 1.0 ;
    }

    if (m_host->getFlags().showCollisionOperator >= 2){
      ctest << "Collision operator of " << m_host->getName() << " after :\n";
      m_egme->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    return true;
  }

  //
  // Calculate raw transition matrix.
  //
  bool gWellProperties::rawTransitionMatrix(double beta, vector<double> &gEne, vector<double> &gDOS, dMatrix *egme)
  {
    // Use number of states to weight the downward transition
    if (m_host->getFlags().useDOSweightedDT){
      // The collision operator.
      for ( int i = 0 ; i < m_ncolloptrsize ; ++i ) {
        double ei = gEne[i] ;
        double ni = gDOS[i] ;
        for ( int j = i ; j < m_ncolloptrsize ; ++j ) {
          double ej = gEne[j];
          double nj = gDOS[j];
          // Transfer to lower Energy -
          // double transferDown = exp(-alpha*(ej - ei)) * (ni/nj);
          // (*m_egme)[i][j] = transferDown;
          double transferDown = m_pEnergyTransferModel->calculateTransitionProbability(ej,ei) ;
          (*egme)[i][j] = transferDown * (ni/nj);

          // Transfer to higher Energy (via detailed balance) -
          // double transferUp = exp(-(alpha + beta)*(ej - ei));
          // (*m_egme)[j][i] = transferUp;
          (*egme)[j][i] = transferDown * exp(-beta*(ej - ei));
        }
      }
    } else {
      // The collision operator.
      for ( int i = 0 ; i < m_ncolloptrsize ; ++i ) {
        double ei = gEne[i] ;
        double ni = gDOS[i] ;
        for ( int j = i ; j < m_ncolloptrsize ; ++j ) {
          double ej = gEne[j];
          double nj = gDOS[j];
          // Transfer to lower Energy -
          double transferDown = m_pEnergyTransferModel->calculateTransitionProbability(ej,ei) ;
          (*egme)[i][j] = transferDown;

          // Transfer to higher Energy (via detailed balance) -
          (*egme)[j][i] = transferDown * (nj/ni) * exp(-beta*(ej - ei)) ;
        }
      }
    }

    return true ;

  }

  //
  // Calculate collision frequency.
  //
  double gWellProperties::collisionFrequency(double beta, const double conc, Molecule *pBathGasMolecule)
  {
    //
    // Lennard-Jones Collision frequency. The collision integral is calculated
    // using the formula of Neufeld et al., J.C.P. Vol. 57, Page 1100 (1972).
    // CONCentration is in molec/cm^3.
    //

    double A = 1.16145 ;
    double B = 0.14874 ;
    double C = 0.52487 ;
    double D = 0.77320 ;
    double E = 2.16178 ;
    double F = 2.43787 ;

    double temp = 1.0/(boltzmann_RCpK*beta) ;

    // Calculate collision parameter averages.
    double bthMass = 0.0;
    bthMass = pBathGasMolecule->getStruc().getMass();

    double bthSigma = 0.0;
    bthSigma = pBathGasMolecule->getBath().getSigma();

    if (!bthSigma)
      cerr << "me:sigma is necessary for " << pBathGasMolecule->getName()
      << ". Correct input file to remove this error." << endl;

    double bthEpsilon = 0.0;
    bthEpsilon = pBathGasMolecule->getBath().getEpsilon();

    if (!bthEpsilon)
      cerr << "me:epsilon is necessary for " << pBathGasMolecule->getName()
      << ". Correct input file to remove this error.";
    double mu   = amu * m_host->getStruc().getMass() * bthMass/(m_host->getStruc().getMass() + bthMass) ;
    double eam  = sqrt(m_host->getBath().getEpsilon() * bthEpsilon) ;
    double sam  = (m_host->getBath().getSigma() + bthSigma) * 0.5;
    double tstr = temp / eam;

    // Calculate collision integral.
    double collFrq = A * exp(-log(tstr) * B) + C * exp(-D * tstr) + E * exp(-F * tstr) ;

    // Calculate molecular collision frequency.
    collFrq *= (M_PI * sam * sam * 1.0e-20 * sqrt(8. * boltzmann_C * temp/(M_PI * mu))) ;
    // Calculate overall collision frequency.
    collFrq *= (conc * 1.0e6) ;

    return collFrq;
  }

  //
  // Calculate a reaction matrix element.
  //
  qd_real gWellProperties::matrixElement(int eigveci, int eigvecj, std::vector<double> &k) const
  {
    // Calculate matrix element starting with the higher energy
    // elements first in order to preserve precision as much as possible.
    qd_real sum = 0.0 ;
    for (int i = m_ncolloptrsize - 1 ; i >= 0 ; --i){
      sum +=  qd_real(k[i]) * ((*m_egvec)[i][eigveci]*(*m_egvec)[i][eigvecj]) ;
    }
    return sum ;
  }

  //
  // Accessor a collision operator eigenvector.
  //
  void gWellProperties::eigenVector(int eigveci, std::vector<double> &evec) const
  {
    // evec.clear() ;
    for (int i(0) ; i < m_ncolloptrsize ; ++i){
      evec[i] =  to_double((*m_egvec)[i][eigveci]) ;
    }
  }

  //
  // Copy collision operator to diagonal block of system matrix.
  //
  void gWellProperties::copyCollisionOperator(qdMatrix *CollOptr,
    const int size,
    const int locate,
    const double RducdOmega) const
  {
    // Find size of system matrix.

    int smsize = static_cast<int>(CollOptr->size()) ;
    //int MaximumGrain = m_host->getEnv().MaxGrn;

    // Check there is enough space in system matrix.

    if (locate + size > smsize)
      throw (std::runtime_error("Error in the size of the system matrix."));

    // Copy collision operator to the diagonal block indicated by "locate"
    // and multiply by the reduced collision frequencey.

    for (int i(0) ; i < size ; ++i) {
      int ii(locate + i) ;
      for (int j(0) ; j < size ; ++j) {
        int jj(locate + j) ;
        (*CollOptr)[ii][jj] = RducdOmega * (*m_egme)[i][j] ;
      }
    }
  }

  //
  // Copy eigenvalues to diagonal elements of the reaction operator
  // matrix in the contracted basis representation.
  //
  void gWellProperties::copyCollisionOperatorEigenValues(qdMatrix *CollOptr,
    const int locate,
    const double Omega) const
  {
    // Check that the contracted basis method has been specifed.

    if (!m_host->getEnv().useBasisSetMethod)
      throw (std::runtime_error("Error: Contracted basis representation not requested."));

    // Find size of system matrix.

    int smsize = static_cast<int>(CollOptr->size()) ;
    int nbasis = static_cast<int>(get_nbasis()) ;

    // Check there is enough space in system matrix.

    if (locate + nbasis > smsize)
      throw (std::runtime_error("Error in the size of the reaction operator matrix in contracted basis representation."));

    // Copy collision operator eigenvalues to the diagonal elements indicated
    // by "locate" and multiply by the reduced collision frequencey.

    for (int i(0) ; i < nbasis ; ++i) {
      int ii(locate + i) ;
      (*CollOptr)[ii][ii] = Omega * m_egval[m_ncolloptrsize - i - 1] ;
    }
  }

  //
  // calculates p(E)*exp(-EB)
  //
  void gWellProperties::grainDistribution(vector<double> &grainFrac, const int totalGrnNumber)
  {
    // If density of states have not already been calcualted then do so.
    if (!m_host->getDOS().calcDensityOfStates())
      cerr << "Failed calculating DOS";

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    if (m_grainDist.size() != gDOS.size() || m_host->getEnv().beta != m_grainFracBeta){
      m_pDistributionCalculator->calculateDistribution(m_host, m_grainDist);
      m_grainFracBeta = m_host->getEnv().beta;
    }

    for (int i = 0; i < totalGrnNumber; ++i){
      grainFrac.push_back(m_grainDist[i]);
    }
  }

  //
  // Get normalized grain distribution.
  //
  void gWellProperties::normalizedInitialDistribution(vector<double> &grainFrac, const int totalGrnNumber, const int numberOfGroupedGrains)
  {
    vector<double> tempGrainFrac;
    grainDistribution(tempGrainFrac, totalGrnNumber);

    const double firstPartition = tempGrainFrac[0];
    double prtfn(firstPartition);
    grainFrac.push_back(firstPartition);
    for (int i = 1; i < totalGrnNumber; ++i){
      prtfn += tempGrainFrac[i];
      if (i < numberOfGroupedGrains){
        grainFrac[0] += tempGrainFrac[i];
      } else {
        grainFrac.push_back(tempGrainFrac[i]);
      }
    }

    for (size_t i(0) ; i < grainFrac.size() ; ++i) {
      grainFrac[i] /= prtfn;
    }

    if (m_host->getFlags().grainBoltzmannEnabled){
      ctest << "\nGrain fraction:\n{\n";
      for (int i = 0; i < totalGrnNumber; ++i){
        ctest << grainFrac[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //
  // Get normalized grain distribution.
  //
  void gWellProperties::normalizedGrnBoltzmannDistribution(vector<double> &grainFrac, const int totalGrnNumber, const int numberOfGroupedGrains)
  {
    // If density of states have not already been calcualted then do so.
    if (!m_host->getDOS().calcDensityOfStates())
      cerr << "Failed calculating DOS";

    vector<double> tempGrnFrac;
    grainFrac.clear();

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    double prtfn(0.0);
    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    const double firstPartition = exp(log(gDOS[0]) - m_host->getEnv().beta * gEne[0] + 10.0);
    tempGrnFrac.push_back(firstPartition);
    prtfn = firstPartition;
    for (int i = 1; i < totalGrnNumber; ++i) {
      const double thisPartition = exp(log(gDOS[i]) - m_host->getEnv().beta * gEne[i] + 10.0);
      prtfn += thisPartition;
      if (i < numberOfGroupedGrains){
        tempGrnFrac[0] += thisPartition;
      }
      else{
        tempGrnFrac.push_back(thisPartition);
      }
    }

    const int tempGrnFracSize = int(tempGrnFrac.size());
    for (int i = 0; i < tempGrnFracSize; ++i){
      tempGrnFrac[i] /= prtfn;
    }

    grainFrac = tempGrnFrac;

  }

  //
  // Accessor for number of basis functions to be used in contracted basis set method.
  //
  size_t gWellProperties::get_nbasis() const { return m_host->getEnv().nBasisSet ; }


  gStructure::gStructure(mesmer::Molecule *pMol) : m_MolecularWeight(-1), m_HasCoords(false)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element
    double MW = ppPropList->XmlReadPropertyDouble("me:MW", optional);
    if(IsNan(MW))
    {
      ReadStructure();
      if(Atoms.empty())
        cerr << "If no chemical structure is provided,"
        "Molecular Weight must be input as an XML property." << endl;
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
      if(! pDOS1.getCellDensityOfStates(rct1CellDOS)|| !pDOS2.getCellDensityOfStates(rct2CellDOS))
        return false;
      std::vector<double> rotConsts;
      if (pDOS1.get_rotConsts(rotConsts) == UNDEFINED_TOP){
        rctsCellDOS = rct2CellDOS;
      }
      else if(pDOS2.get_rotConsts(rotConsts) == UNDEFINED_TOP){
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
    if(!Atoms.empty())
      return m_HasCoords;
    PersistPtr ppMol = getHost()->get_PersistentPointer();
    PersistPtr ppAtom = ppMol->XmlMoveTo("atomArray");
    if(!ppAtom) // there may not be an <atomArray> element
      ppAtom = ppMol; 
    while(ppAtom=ppAtom->XmlMoveTo("atom"))
    {
      atom at;
      const char* el = ppAtom->XmlReadValue("elementType");
      if(!el)
      {
        cerr << "<atom> elements must have an elementType attribute" << endl;
        return false;
      }
      at.element = el;
      const char* pId = ppAtom->XmlReadValue("id", optional);
      at.id = (pId)? pId : at.element ;
      double x3, y3, z3;
      x3 = ppAtom->XmlReadDouble("x3", optional);
      y3 = ppAtom->XmlReadDouble("y3", optional);
      z3 = ppAtom->XmlReadDouble("z3", optional);
      if(!IsNan(x3) && !IsNan(y3) && !IsNan(z3))
      {
        at.coords.Set(x3, y3, z3);
        if(x3!=0 || y3!=0 || z3!=0)
          m_HasCoords = true; //at least one atom with non-zero coordinates
      }
      Atoms[at.id] = at;
    }

    //Read all the bonds. For each bond add a connect to each atom
    PersistPtr ppBond = ppMol->XmlMoveTo("bondArray");
    int ibond=1;
    if(!ppBond) // there may not be an <bondArray> element
      ppBond = ppMol;
    while(ppBond=ppBond->XmlMoveTo("bond"))
    {
      const char* pId = ppBond->XmlReadValue("id", optional);
      string id;
      if(pId)
        id = pId;
      else
      {
        //id is e.g. "bond3", if not provided
        stringstream ss;
        ss << " bond" << ibond;
        id=ss.str();
      }

      const char* pRefs = ppBond->XmlReadValue("atomRefs2");
      if(!pRefs) return false;
      string refs(pRefs);
      string::size_type pos = refs.find_first_of(" ,");
      string::size_type pos2 = refs.find_last_of(" ,");
      if(pos==string::npos) return false;
      string atomref1 = refs.substr(0, pos);
      string atomref2 = refs.substr(pos2+1);
      Bonds[id] = make_pair(atomref1,atomref2);
      Atoms[atomref1].connects.push_back(atomref2);
      Atoms[atomref2].connects.push_back(atomref1);
      ++ibond;
    }

    return m_HasCoords;
  }

  double gStructure::CalcMW()
  {
    map<string, atom>::iterator iter;
    double MW = 0.0;
    for(iter=Atoms.begin(); iter!=Atoms.end(); ++iter)
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
    for(coniter=Atoms[atomID].connects.begin(); coniter!=Atoms[atomID].connects.end();++coniter)
    {
      if(find(atomset.begin(), atomset.end(), *coniter)!=atomset.end())
        continue;
      GetAttachedAtoms(atomset, *coniter);    
    }
  }

  double gStructure::CalcMomentAboutAxis(vector<string> atomset, vector3 at1, vector3 at2)
  {
    double sumMoment = 0.0;
    vector<string>::iterator iter;
    for(iter=atomset.begin(); iter!=atomset.end(); ++iter)
    {
      vector3 a = Atoms[*iter].coords;
      double d = Point2Line(a, at1, at2);
      sumMoment += atomMass(Atoms[*iter].element) * d * d;
    }
    return sumMoment;
  }

  //Returns the rotational constants (in cm-1) in a vector
  //OK for atoms and diatomics but currently no recognition of symmetry
  vector<double> gStructure::CalcRotConsts()
  {
    vector<double> RotConsts(3, 0.0); //cm-1
    if(NumAtoms()<2)
      return RotConsts; //empty
    //Determine centre of mass
    map<string, atom>::iterator iter;
    vector3 centreOfMass; 
    double mt = 0.0;
    for(iter=Atoms.begin(); iter!=Atoms.end(); ++iter)
    {
      double mass = atomMass(iter->second.element);
      centreOfMass += iter->second.coords * mass;
      mt  += mass;
    }
    centreOfMass /= mt ;

    dMatrix MI(3);
    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
    for(iter=Atoms.begin(); iter!=Atoms.end(); ++iter)
    {
      vector3 c = iter->second.coords - centreOfMass;
      double  m = atomMass(iter->second.element);
      sxx += m * c.x() * c.x();
      syy += m * c.y() * c.y();
      szz += m * c.z() * c.z();
      sxy += m * c.x() * c.y();
      sxz += m * c.x() * c.z();
      syz += m * c.y() * c.z();
    }

    vector<double> PrincipalMI(3, 0.0);//initially amuAng2, eventually gcm2
    if(NumAtoms()==2)
      PrincipalMI[0] = szz;
    else
    {
      MI[0][0] = syy+szz;
      MI[0][1] = -sxy;
      MI[0][2] = -sxz; 
      MI[1][0] = -sxy;
      MI[1][1] = sxx+szz;
      MI[0][2] = -syz; 
      MI[2][0] = -sxz;
      MI[2][1] = -syz; 
      MI[2][2] = sxx + syy;

      MI.diagonalize(&PrincipalMI[0]);
    }

    const double amuA2TOgcm2 = 1.0E-16/AvogadroC;
    for(unsigned i=0; i<3;++i)
    {
      RotConsts[i] = PrincipalMI[i]==0.0 ? 0.0 : conMntInt2RotCnt/PrincipalMI[i];
      PrincipalMI[i] *= amuA2TOgcm2;
    }

    return RotConsts;

  }

  double gStructure::CalcSumEMinusHf0(bool UsingAtomBasedThermo, bool useHf298)
  {
    //calculate for each atom (ab initio E - Hf0) and return sum
    if(!ReadStructure())
    {
      cerr << "To use me::Hf0 or Hf298 the molecule needs chemical structure (an atomList at least)" << endl;
      return false;
    }
    double sum = 0.0;
    map<string, double> atomdiffs; //el symbol, diff
    map<string, atom>::iterator iter;
    for(iter=Atoms.begin(); iter!=Atoms.end(); ++iter)
    {
      string el =iter->second.element;
      if(atomdiffs.find(el)==atomdiffs.end())
      {
        //get vals from librarymols.xml
        PersistPtr ppMol = GetFromLibrary(el, PersistPtr());
        double diff;
        if(ppMol)
        {
          diff = ppMol->XmlReadPropertyDouble("me:ZPE",optional);
          if(useHf298)
          {
            diff -= ppMol->XmlReadPropertyDouble("me:Hf298",optional);
            diff += ppMol->XmlReadPropertyDouble("me:H0-H298",optional);
          }
          else if(!UsingAtomBasedThermo)
            diff -= ppMol->XmlReadPropertyDouble("me:Hf0",optional);
        }
        if(!ppMol || IsNan(diff))
        {
          cerr << "The value of Hf0 for " << getHost()->getName()
            << " will be incorrect because one or more of its elements"
            << " was not in the library, or lacked me:ZPE and me:Hf0 or me:Hf298 properties" << endl;
          return 0.0;
        }
        atomdiffs[el] = diff; //save diff for this el in atomdiffs
        sum += diff;
      }
      else //diff for this el already known
        sum += atomdiffs[el];
    }
    return sum;
  }

}//namespace
