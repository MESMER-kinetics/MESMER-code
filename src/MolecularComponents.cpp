// MolecularComponents.cpp
//
// Author: Chi-Hsiu Liang
//
//-------------------------------------------------------------------------------------------
#include <stdexcept>
#include "Molecule.h"

using namespace std ;
using namespace Constants ;

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
    if (m_Sigma_chk == 0){
      cinfo << "m_Sigma is provided but not used in " << m_host->getName() << "." << endl;
    }
    if (m_Epsilon_chk == 0){
      cinfo << "m_Epsilon is provided but not used in " << m_host->getName() << "." << endl;
    }
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
      //exit(1);
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
      //exit(1);
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
    if (m_RC_chk == 0) cinfo << "Rotational constants are provided but not used in " << m_host->getName() << "." << endl;
    if (m_Sym_chk == 0) cinfo << "m_Sym is provided but not used in " << m_host->getName() << "." << endl;
    if (m_ZPE_chk == 0) cinfo << "m_ZPE is provided but not used in " << m_host->getName() << "." << endl;
    if (m_scaleFactor_chk == 0) cinfo << "m_scaleFactor is provided but not used in " << m_host->getName() << "." << endl;
    if (m_SpinMultiplicity_chk == 0) cinfo << "m_SpinMultiplicity is provided but not used in " << m_host->getName() << "." << endl;
    if (m_VibFreq_chk == 0) cinfo << "m_VibFreq is provided but not used in " << m_host->getName() << "." << endl;

    // Free any memory assigned for calculating densities of states. (must be in reverse order)
    if (m_grainDOS.size()) m_grainDOS.clear();
    if (m_grainEne.size()) m_grainEne.clear();
    if (m_cellDOS.size()) m_cellDOS.clear();
    if (m_VibFreq.size()) m_VibFreq.clear();
    if (m_eleExc.size()) m_eleExc.clear();
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
    m_grainDOS()
  {
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    bool hasVibFreq = true; bool hasRotConst = true;
    txt= ppPropList->XmlReadProperty("me:vibFreqs", optional);
    if(!txt){
      hasVibFreq = false;
      cinfo << "Cannot find argument me:vibFreqs. Assumes that it is an atom or atomic ion." << endl;
      m_VibFreq_chk = -1;
    }
    else { istringstream idata(txt); double x; while (idata >> x) m_VibFreq.push_back(x); m_VibFreq_chk = 0;}

    txt= ppPropList->XmlReadProperty("me:rotConsts", optional);
    if(!txt){
      hasRotConst = false;
      cinfo << "Cannot find argument me:rotConsts. Assumes that it is an atom or atomic ion." << endl;
      m_RC_chk = -1;
    }
    else {
      istringstream idata(txt);
      std::vector<double> rCnst(3);
      idata >> rCnst[0]
      >> rCnst[1]
      >> rCnst[2];
      rCnst[0] = abs(rCnst[0]);
      rCnst[1] = abs(rCnst[1]);
      rCnst[2] = abs(rCnst[2]);
      std::sort(rCnst.begin(), rCnst.end()); // make sure the rotational constants are in ascending order.
      m_RotCstA = rCnst[2]; // order the rotational constants in descending order (largest in the beginning)
      m_RotCstB = rCnst[1];
      m_RotCstC = rCnst[0];
      m_RC_chk = 0;
    }

    if (hasVibFreq != hasRotConst){
      cerr << "Improper setting on vibrational frequencies or rotational constants. Check input file to remove this error.";
    }

    txt= ppPropList->XmlReadProperty("me:electronicExcitation", optional);
    if(txt) {
      istringstream idata(txt);
      m_eleExc.clear();
      double _iele;
      while (idata >> _iele) m_eleExc.push_back(_iele);
    }

    m_Sym = ppPropList->XmlReadPropertyDouble("me:symmetryNumber");
    m_Sym_chk = 0;

    double tempzpe = ppPropList->XmlReadPropertyDouble("me:ZPE");
    string unitsInput = ppPropList->XmlReadPropertyAttribute("me:ZPE", "units");
    txt= ppPropList->XmlReadPropertyAttribute("me:ZPE", "convention", optional);
    m_EnergyConvention = txt ? txt : "arbitary";

    const char* pLowertxt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "lower", optional);
    const char* pUppertxt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "upper", optional);
    const char* pStepStxt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "stepsize", optional);
    double value(getConvertedEnergy(unitsInput, tempzpe));
    if (pLowertxt && pUppertxt){
      double tempLV(0.0), tempUV(0.0), tempSS(0.0);
      stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt);
      s3 >> tempLV; s4 >> tempUV; s5 >> tempSS;
      double valueL(getConvertedEnergy(unitsInput, tempLV)),
        valueU(getConvertedEnergy(unitsInput, tempUV)),
        stepsize(getConvertedEnergy(unitsInput, tempSS));
      set_zpe(valueL, valueU, stepsize);
    }
    else{
      set_zpe(value);
    }
    m_ZPE_chk = 0;

    // The reason why me:frequenciesScaleFactor stands out to be a separate property in the propertyList is that
    // this value is not usually necessary. The default value is 1.0 and it is usually the case.
    m_scaleFactor = ppPropList->XmlReadPropertyDouble("me:frequenciesScaleFactor");
    m_scaleFactor_chk = 0;



    m_SpinMultiplicity = ppPropList->XmlReadPropertyInteger("me:spinMultiplicity");
    m_SpinMultiplicity_chk = 0;

    /* For molecular energy me:ZPE is used if it is present. If it is not, a value
    calculated from meHf298 is used and written back to the datafile as a me:ZPE
    property. It is consequently used in the next run and available to be varied
    or optimized. The original calculated value remains recorded in an attribute.
    */
    txt = ppPropList->XmlReadProperty("me:Hf298", false);
    if(txt && m_ZPE_chk < 0){ //ignore if there is a me:ZPE
      double Hf298;
      istringstream idata(txt);
      idata >> Hf298; //orig units
      const char* utxt= ppPropList->XmlReadPropertyAttribute("me:ZPE", "units");
      utxt = utxt ? utxt : "kJ/mol";
      Hf298 = getConvertedEnergy(utxt, Hf298); //cm-1

      //Calculate ZPE from Thermodynamic Heat of Formation

      //*** This is INCOMPLETE and will calculate an incorrect result. NEEDS REVISITING ***

      //The general way is to calculate dln(rot/vib partition function)/dBeta + 1.5kT
      //calcDensityOfStates(); //but this calls get_zpe and ZPE hasn't been set yet
      //double Z = max(canonicalPartitionFunction(m_grainDOS, m_grainEne, boltzmann_RCpK * 298), 1.0);

      //Vibrations are treated classically, rotational constants are considered small:
      // 0.5kT for linear mols, 1.5kT for non-linear polyatomics
      std::vector<double> rotConsts;
      int ret = get_rotConsts(rotConsts);
      double Hf0 = Hf298 - 0.5 * boltzmann_RCpK * 298 *(3 + (ret==0) + 2*(ret==2)); //***Vib TODO
      set_zpe(Hf0); //cm-1
      m_ZPE_chk=0;

      //Write the converted value back to a me:ZPE element in the XML file
      stringstream ss;
      ss << ConvertFromWavenumbers(utxt, Hf0);
      PersistPtr ppScalar = ppPropList->XmlWriteProperty("me:ZPE", ss.str(), utxt);
      ppScalar->XmlWriteAttribute("convention", "thermodynamic");//orig units
      ppScalar->XmlWriteAttribute("origValue", ss.str());
      m_EnergyConvention = "thermodynamic";
      cinfo << "New me:ZPE element written with data from me:Hf298" << endl;
    }
    else if(m_ZPE_chk < 0)
      cwarn << "No energy specified (as me:ZPE or me:Hf298 properties)" << endl;
  }

  //
  // Get cell density of states.
  //
  void gDensityOfStates::getCellDensityOfStates(vector<double> &cellDOS, int startingCell) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";
    if (startingCell == 0)
      cellDOS = m_cellDOS;
    else{
      int MaximumCell = m_host->getEnv().MaxCell;
      for (int i(startingCell); i < MaximumCell; ++i){
        cellDOS.push_back(m_cellDOS[i]);
      }
    }
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

  int gDensityOfStates::test_rotConsts()
  {
    std::vector<double> mmtsInt;
    return get_rotConsts(mmtsInt);
  }

  int gDensityOfStates::get_rotConsts(std::vector<double> &mmtsInt)
  {
    if (m_RC_chk == -1){
      cinfo << "Rotational constants were not defined but requested." << endl;
      --m_RC_chk;
      return -4; // treat as a non-rotor
    }
    else if (m_RC_chk < -1){
      --m_RC_chk;
      return -4;
    }
    mmtsInt.clear();
    mmtsInt.push_back(m_RotCstA);
    mmtsInt.push_back(m_RotCstB);
    mmtsInt.push_back(m_RotCstC);
    /* now the classification of rotors is simplified to only three following types. 3-D rotors may have other
    attributes different from one another but in ILT they are treated as the same type. The function return values
    are temporary shorthand representations. */
    ++m_RC_chk;
    if      ((mmtsInt[0] + mmtsInt[1] + mmtsInt[2]) == 0.) return -4; // not a rotor
    else if ((mmtsInt[0] * mmtsInt[1] * mmtsInt[2]) == 0.) return  0; // 2-D linear
    else                                                   return  2; // 3-D symmetric/asymmetric/spherical top
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

    // Determine the method of rotational DOS calculation.
    PersistPtr pp = m_host->get_PersistentPointer();
    PersistPtr ppRotDOSCMethod = pp->XmlMoveTo("me:DOSCMethod");
    DensityOfStatesCalculator* pRotationalDOSCalculator(NULL);

    while (ppRotDOSCMethod){

      const char* pRotDOSCMethodTypetxt = ppRotDOSCMethod->XmlReadValue("me:DOSCType") ;

      if (strcmp(pRotDOSCMethodTypetxt, "external") == 0){
        const char* pRotDOSCMethodtxt = ppRotDOSCMethod->XmlRead();
        pRotationalDOSCalculator = DensityOfStatesCalculator::Find(pRotDOSCMethodtxt);
        if(!pRotationalDOSCalculator) // if the provided method cannot be found,
        {
          cinfo << "Unknown method " << pRotDOSCMethodtxt  << " for the calculation of rotational DOS.\n"
            << "Please check spelling error. Default method <Classical rotors> is used." << endl;
        }
        break;
      }

      ppRotDOSCMethod = ppRotDOSCMethod->XmlMoveTo("me:DOSCMethod");
    }

    if (!pRotationalDOSCalculator){
      pRotationalDOSCalculator = DensityOfStatesCalculator::Find("Classical rotors");
    }

    if (!pRotationalDOSCalculator->countCellDOS(this, MaximumCell)){
      return false;
    }

    //-------------------------------------------------------------
    // Find and calculate any other DOS calculators in the molecule
    // (internal rotations and other vibrational treatments etc.)
    PersistPtr ppDOSCMethod = pp->XmlMoveTo("me:DOSCMethod");
    DensityOfStatesCalculator* pInternalDOSCalculator(NULL);

    while (ppDOSCMethod){
      const char* pDOSCMethodTypetxt = ppDOSCMethod->XmlReadValue("me:DOSCType", false);
      const char* pDOSCMethodtxt = ppDOSCMethod->XmlReadValue("me:DOSCMethod", false);

      if (strcmp(pDOSCMethodTypetxt, "external") == 0){
        ppDOSCMethod = ppDOSCMethod->XmlMoveTo("me:DOSCMethod");
        continue;
      }

      pInternalDOSCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);

      if (!pInternalDOSCalculator->countCellDOS(this, MaximumCell, ppDOSCMethod)){
        return false;
      }

      ppDOSCMethod = ppDOSCMethod->XmlMoveTo("me:DOSCMethod");
    }
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

    string comment("Rovibronic partition function calculation at various temperatures. qtot : partition function as a product of quantum mechanical partition functions for vibrations (1-D harmonic oscillator) and classifical partition functions for rotations.  sumc : (user calculated) cell based partition function. sumg : (user calculated) grain based partition function ");

    PersistPtr ppList = m_host->get_PersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment );

    if (m_host->getFlags().testDOSEnabled) ctest << endl << "Test rovibronic density of states for: " << m_host->getName() << "\n{\n";
    if (m_host->getFlags().testDOSEnabled) ctest << "      T           qtot           sumc           sumg\n";

    //loop through predefined test temperatures
    for ( int n = 0 ; n < 29 ; ++n ) {
      double temp = 100.0*static_cast<double>(n + 2) ;
      double beta = 1.0/(boltzmann_RCpK*temp) ;

      // Calculate rovibronic partition functions based on cells.
      // The following catches the case where the molecule is a single atom
      double cellCanPrtnFn = max(canonicalPartitionFunction(m_cellDOS, cellEne, beta), 1.0) ;
      if (cellCanPrtnFn == 1.0){
        // Electronic partition function for atom is accounted here.
        cellCanPrtnFn = double(getSpinMultiplicity()) ;
      }

      // Calculate rovibronic partition functions based on grains.
      // The following catches the case where the molecule is a single atom
      double grainCanPrtnFn = max(canonicalPartitionFunction(m_grainDOS, m_grainEne, beta), 1.0) ;
      if (grainCanPrtnFn == 1.0){
        // Electronic partition function for atom is accounted here.
        grainCanPrtnFn = double(getSpinMultiplicity()) ;
      }

      // Calculate rovibronic partition functions using analytical formula (treat vibrations classically).
      double qtot = 1.0 ;
      vector<double> rotConst; int rotorType;
      rotorType = get_rotConsts(rotConst);

      vector<double> vibFreq; get_VibFreq(vibFreq);

      switch(rotorType){
        case 2://3-D symmetric/asymmetric/spherical top
          for ( vector<double>::size_type j = 0 ; j < vibFreq.size() ; ++j ) {
            qtot /= (1.0 - exp(-beta*vibFreq[j])) ;
          }
          qtot *= (sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/get_Sym()) ;
          break;
        case 0://2-D linear
          for ( vector<double>::size_type j = 0 ; j < vibFreq.size() ; ++j ) {
            qtot /= (1.0 - exp(-beta*vibFreq[j])) ;
          }
          qtot /= rotConst[0]* get_Sym()*beta ;
          break;
        default:
          qtot = 0.;
      }
      qtot *= double(getSpinMultiplicity());
      qtot = max(qtot, 1.0);
      if (qtot == 1.0){
        // Electronic partition function for atom is accounted here.
        qtot = double(getSpinMultiplicity()) ;
      }

      if (m_host->getFlags().testDOSEnabled) formatFloat(ctest, temp,  6,  7) ;
      if (m_host->getFlags().testDOSEnabled) formatFloat(ctest, qtot,  6, 15) ;
      if (m_host->getFlags().testDOSEnabled) formatFloat(ctest, cellCanPrtnFn,  6, 15) ;
      if (m_host->getFlags().testDOSEnabled) formatFloat(ctest, grainCanPrtnFn,  6, 15) ;
      if (m_host->getFlags().testDOSEnabled) ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
      ppItem->XmlWriteValueElement("me:T",    temp, 6);
      ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
      ppItem->XmlWriteValueElement("me:sumc", cellCanPrtnFn, 6);
      ppItem->XmlWriteValueElement("me:sumg", grainCanPrtnFn, 6);
    }
    if (m_host->getFlags().testDOSEnabled) ctest << "}" << endl;

    if (m_host->getFlags().cellDOSEnabled){
      ctest << endl << "Cell rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i = 0; i < MaximumCell; ++i){
        formatFloat(ctest, cellEne[i],  6,  15) ;
        formatFloat(ctest, m_cellDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }

    if (m_host->getFlags().grainDOSEnabled){
      ctest << endl << "Grain rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i = 0; i < MaximumGrain; ++i){
        formatFloat(ctest, m_grainEne[i],  6,  15) ;
        formatFloat(ctest, m_grainDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }
  }


  double gDensityOfStates::get_zpe() {
    if (m_ZPE_chk == -1) {
      cinfo << "m_ZPE was not defined but requested in " << m_host->getName() << ". Default value " << m_ZPE.get_value() << " is given." << endl;
      --m_ZPE_chk;
      double zpe = m_ZPE.get_value();
      return zpe;
    }
    else if (m_ZPE_chk < -1){
      --m_ZPE_chk;
      double zpe = m_ZPE.get_value();
      return zpe;
    }
    ++m_ZPE_chk;
    double zpe = m_ZPE.get_value();
    return zpe;
  }



  double gDensityOfStates::get_scaleFactor() {
    if (m_scaleFactor_chk == -1){
      cinfo << "m_scaleFactor was not defined but requested in " << m_host->getName() << ". Default value " << m_scaleFactor << " is given." << endl;
      --m_scaleFactor_chk;
      return m_scaleFactor;
    }
    else if (m_scaleFactor_chk < -1){
      --m_scaleFactor_chk;
      return m_scaleFactor;
    }
    ++m_scaleFactor_chk;
    return m_scaleFactor ;
  }

  double gDensityOfStates::get_Sym(void){
    if (m_Sym_chk == -1){
      cinfo << "m_Sym was not defined but requested in " << m_host->getName() << ". Default value " << m_Sym << " is given." << endl;
      --m_Sym_chk;
      return m_Sym;
    }
    else if (m_Sym_chk < -1){
      --m_Sym_chk;
      return m_Sym;
    }
    ++m_Sym_chk;
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

  int gDensityOfStates::getSpinMultiplicity(){
    if (m_SpinMultiplicity_chk >= 0){
      ++m_SpinMultiplicity_chk;
      return m_SpinMultiplicity ;
    }
    else{
      cinfo << "m_SpinMultiplicity was not defined but requested in " << m_host->getName() << ". Default value " << m_SpinMultiplicity << " is given." << endl;
      return m_SpinMultiplicity;
    }
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
      cerr << "Failed calculating DOS.";
      exit(1);
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

    // Calculate the rovibronic partition function based on the grain DOS
    // The following catches the case where the molecule is a single atom
    double CanPrtnFn = max(canonicalPartitionFunction(m_grainDOS, m_grainEne, m_host->getEnv().beta), 1.0) ;
    if (CanPrtnFn == 1.0){
      // Electronic partition function for atom is accounted here.
      CanPrtnFn = double(getSpinMultiplicity()) ;
    }

    return CanPrtnFn ;
  }

  //-------------------------------------------------------------------------------------------------
  // Transition state related properties
  //-------------------------------------------------------------------------------------------------

  //
  // Constructor, destructor and initialization
  //
  gTransitionState::~gTransitionState()
  {
    if (m_ImFreq_chk == 0) cinfo << "m_ImFreq is provided but not used in " << m_host->getName() << "." << endl;
  }

  gTransitionState::gTransitionState(Molecule* pMol)
    :m_ImFreq(0.0),
    m_ImFreq_chk(-1)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    bool hasImFreq = false;
    const char* txt;
    txt= ppPropList->XmlReadProperty("me:imFreqs", optional);
    if(!txt){
      cinfo << "No imaginary vibrational frequency.\n";
      m_ImFreq_chk = -1;
    }
    else {
      istringstream idata(txt); double x;
      while (idata >> x) m_ImFreq = x;
      hasImFreq = true;
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
  gWellProperties::~gWellProperties()
  {
    if (m_DeltaEdown_chk == 0){
      cinfo << "m_DeltaEdown is provided but not used in " << m_host->getName() << "." << endl;
    }
    if (m_egme != NULL) delete m_egme ;
    if (m_grainDist.size()) m_grainDist.clear();
  }

  gWellProperties::gWellProperties(Molecule* pMol)
    :m_DeltaEdownExponent(0.0),
    m_DeltaEdownRefTemp(298.0),
    m_DeltaEdown(0.0),
    m_collisionFrequency(0.0),
    m_ncolloptrsize(0),
    m_lowestBarrier(9e23),
    m_numGroupedGrains(0),
    m_pDistributionCalculator(NULL),
    m_DeltaEdown_chk(-1),
    m_grainFracBeta(0.),
    m_grainDist(0),
    m_egme(NULL),
    m_egvec(NULL),
    m_egval(0)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    //TODO make essential for appropriate molecules
    txt= ppPropList->XmlReadProperty("me:deltaEDown");
    if(!txt){
      cinfo << "No me:deltaEDown provided." << endl;
      // deltaEDown is not always necessary. Hoever, it is not wise to provide a default value.
    }
    else {
      istringstream idata(txt);
      double value(0.0);
      idata >> value;
      const char* pLowertxt    = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "lower", optional);
      const char* pUppertxt    = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "upper", optional);
      const char* pStepStxt    = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "stepsize" ,optional);
      const char* pRefTemptxt  = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "referenceTemperature", optional );
      const char* pExponenttxt = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "exponent", optional);
      if (pLowertxt && pUppertxt){
        double valueL(0.0), valueU(0.0), stepsize(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> valueL; s4 >> valueU; s5 >> stepsize;
        setDeltaEdown(valueL, valueU, stepsize);
      }
      else{
        setDeltaEdown(value);
      }
      if(pRefTemptxt){
        double ref_t(298.);
        stringstream s_temp(pRefTemptxt); s_temp >> ref_t;
        setDeltaEdownRefTemp(ref_t);
      }
      if(pExponenttxt){
        double ref_exp(0.0);
        stringstream s_exp(pExponenttxt); s_exp >> ref_exp;
        setDeltaEdownExponent(ref_exp);
      }
      m_DeltaEdown_chk = 0;
    }

    // Determine the method of DOS calculation.
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

  double gWellProperties::getDeltaEdown()                    {
    if (m_DeltaEdown_chk >= 0){
      ++m_DeltaEdown_chk;
      const double refTemp = getDeltaEdownRefTemp();
      const double dEdExp = getDeltaEdownExponent();
      const double dEdRef = m_DeltaEdown.get_value();
      const double temperature = 1. / (boltzmann_RCpK * m_host->getEnv().beta);
      return dEdRef * pow((temperature/refTemp),dEdExp);
    }
    else{
      cerr << "m_DeltaEdown was not defined but requested.";
      exit(1);
    }
  } ;

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

      int lowestBarrier = int(getLowestBarrier() / double(m_host->getEnv().GrainSize));

      // Second find out in the current temperature under which energy grain where 99% of the population is located when
      // the system is in equilibrium for this well.
      double popUnder(0.0), totalPartition(0.0);
      int idx(0);
      for (int i(0); i < m_ncolloptrsize; ++i){
        totalPartition += sqrt(exp(log(gDOS[i]) - beta * gEne[i] + 10.0));
      }
      for (; idx < m_ncolloptrsize; ++idx){
        popUnder += sqrt(exp(log(gDOS[idx]) - beta * gEne[idx] + 10.0));
        if (popUnder/totalPartition > .999) break;
      }

      m_numGroupedGrains = 0; // Reset the number of grains grouped into a reservoir grain to zero.
      if (lowestBarrier > idx && m_host->getFlags().doReservoirStateMethod){
        m_numGroupedGrains = idx;
      }
      int reducedCollOptrSize = m_ncolloptrsize - m_numGroupedGrains + 1;
      //-----------------------------------------

      if (m_host->getFlags().doReservoirStateMethod && m_numGroupedGrains){
        if (!collisionOperatorWithReservoirState(beta, reducedCollOptrSize)){
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
    if (m_host->getFlags().doBasisSetMethod) diagonalizeCollisionOperator();

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

    ctest << "\nEigenvectors of: " << m_host->getName() << endl;
    m_egvec->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);

    // The eigenvalues in this series should contain a zero eigenvalue (lambda_0) which addresses the conserveness
    // of the system, the second smallest eigenvalue in magnitude is lamda_1.
    ctest << "\nEigenvalues of: " << m_host->getName() << "\n{\n";
    for(int i(0); i < m_ncolloptrsize; ++i){
      ctest << m_egval[i] << endl;
    }
    ctest << "}\n";
  }

  //
  // Calculate collision operator
  //
  bool gWellProperties::collisionOperatorWithReservoirState(double beta, const int reducedCollOptrSize)
  {
    if (m_host->getDOS().test_rotConsts() < 0) return true;
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

    double DEDown = getDeltaEdown();
    double alpha = 1.0/DEDown ;

    // issue a warning message and exit if delta_E_down is smaller than grain size.
    if (DEDown < double(m_host->getEnv().GrainSize) && !m_host->getFlags().allowSmallerDEDown){
      cerr << "Delta E down is smaller than grain size: the solution may not converge.";
      return false;
    }

    // Allocate memory.
    if (m_egme) delete m_egme ;                       // Delete any existing matrix.
    m_egme = new dMatrix(reducedCollOptrSize) ;

    //---------------------------------------------------------------------------------------------------
    //-------------------- The part doing the same jobs as making a whole collision operator ------------
    dMatrix* tempEGME = new dMatrix(m_ncolloptrsize);

    // Use number of states to weigh the downward transition
    if (m_host->getFlags().useDOSweighedDT){
      // The collision operator.
      for (int i(0) ; i < m_ncolloptrsize ; ++i ) {
        double ei = gEne[i] ;
        double ni = gDOS[i] ;
        for (int j(i) ; j < m_ncolloptrsize ; ++j ) {
          double ej = gEne[j];
          double nj = gDOS[j];
          // Transfer to lower Energy -
          double transferDown = exp(-alpha*(ej - ei)) * (ni/nj);
          (*tempEGME)[i][j] = transferDown;

          // Transfer to higher Energy (via detailed balance) -
          double transferUp = exp(-(alpha + beta)*(ej - ei));
          (*tempEGME)[j][i] = transferUp;
        }
      }
    }
    else{
      // The collision operator.
      for (int i(0) ; i < m_ncolloptrsize ; ++i ) {
        const double ei = gEne[i] ;
        const double ni = gDOS[i] ;
        for (int j(i); j < m_ncolloptrsize ; ++j ) {
          const double ej = gEne[j];
          const double nj = gDOS[j];
          // Transfer to lower Energy -
          const double transferDown = exp(-alpha*(ej - ei)) ;
          (*tempEGME)[i][j] = transferDown;

          // Transfer to higher Energy (via detailed balance) -
          const double transferUp = exp(-alpha*(ej - ei)) * (nj/ni) * exp(-beta*(ej - ei)) ;
          (*tempEGME)[j][i] = transferUp;
        }
      }
    }

    if (m_host->getFlags().showCollisionOperator != 0){
      ctest << "\nCollision operator of " << m_host->getName() << " before normalization:\n";
      tempEGME->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    //Normalisation
    tempEGME->normalizeProbabilityMatrix();

    if (m_host->getFlags().showCollisionOperator >= 1){
      ctest << "\nCollision operator of " << m_host->getName() << " after normalization:\n";
      tempEGME->showFinalBits(0, m_host->getFlags().print_TabbedMatrices);
    }

    //-------------------- The part doing the same jobs as making a whole collision operator ------------
    //---------------------------------------------------------------------------------------------------



    //--------------------------------
    // COPY, SUMMATION AND SUBSTITUTE
    // Need to copy things over, the active states first.
    for (int i(m_numGroupedGrains) ; i < m_ncolloptrsize ; ++i ) {
      for (int j(m_numGroupedGrains); j < m_ncolloptrsize ; ++j ) {
        (*m_egme)[i - m_numGroupedGrains + 1][j - m_numGroupedGrains + 1] = (*tempEGME)[i][j];
      }
    }

    // Sum up the downward transition terms for the reservoir grain
    for (int j(m_numGroupedGrains); j < m_ncolloptrsize ; ++j ) {
      double downwardSum(0.0);
      for (int i(0) ; i < m_numGroupedGrains ; ++i ) {
        downwardSum += (*tempEGME)[i][j];
      }
      (*m_egme)[0][j - m_numGroupedGrains + 1] = downwardSum;
    }

    // Summing the contribution of each grain according to its Boltzmann fraction to active state and reservoir state

    // Get Boltzmann distribution
    vector<double> BoltzmannPopFrac ; // Population fraction of the adduct
    normalizedGrnBoltzmannDistribution(BoltzmannPopFrac, m_numGroupedGrains) ;

    for (int j(0); j < m_numGroupedGrains; ++j){
      for (int i(0); i < m_ncolloptrsize; ++i){
        if (i < m_numGroupedGrains){
          (*m_egme)[0][0] += BoltzmannPopFrac[j] * (*tempEGME)[i][j];
        }
        else{
          (*m_egme)[i - m_numGroupedGrains + 1][0] += BoltzmannPopFrac[j] * (*tempEGME)[i][j];
        }
      }
    }

    vector<double> tempVector;
    for (int i(0); i < reducedCollOptrSize; ++i){
      tempVector.push_back((*m_egme)[i][0]);
    }
    // COPY, SUMMATION AND SUBSTITUTE
    //--------------------------------



    // print out of column sums to check normalization results
    if (m_host->getFlags().reactionOCSEnabled){
      ctest << endl << "Collision operator column Sums" << endl << "{" << endl ;
      for (int i(0) ; i < reducedCollOptrSize ; ++i ) {
        double columnSum(0.0) ;
        for (int j(0) ; j < reducedCollOptrSize ; ++j ){
          columnSum += to_double((*m_egme)[j][i]) ;
        }
        ctest << columnSum << endl ;
      }
      ctest << "}" << endl;
    }

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
    popDist[0] = sqrt(popDist[0]);

    for (int i(1) ; i < reducedCollOptrSize ; ++i ) {
      for (int j(0) ; j < i ; ++j ){
        (*m_egme)[j][i] *= popDist[i]/popDist[j] ;
        (*m_egme)[i][j]  = (*m_egme)[j][i] ;
      }
    }

    //account for collisional loss by subrtacting unity from the leading diagonal.
    for (int i(0) ; i < reducedCollOptrSize ; ++i ) {
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
    if (m_host->getDOS().test_rotConsts() < 0) return true;
    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //

    int i, j;

    double DEDown = getDeltaEdown();
    double alpha = 1.0/DEDown ;

    // issue a warning message and exit if delta_E_down is smaller than grain size.
    if (DEDown < double(m_host->getEnv().GrainSize) && !m_host->getFlags().allowSmallerDEDown){
      cerr << "Delta E down is smaller than grain size: the solution may not converge.";
      return false;
    }

    // Allocate memory.
    if (m_egme) delete m_egme ;                       // Delete any existing matrix.
    m_egme = new dMatrix(m_ncolloptrsize) ;           // Collision operator matrix.

    vector<double> gEne;
    vector<double> gDOS;
    m_host->getDOS().getGrainEnergies(gEne);
    m_host->getDOS().getGrainDensityOfStates(gDOS);

    // Initialisation and error checking.
    for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
      if (gDOS[i] <= 0.0) {
        cerr << "Data indicates that grain " << i << " of the current colliding molecule has no states.";
        cerr << "Usually for small GrainSize DOSs calculated by QM rotor method some grain may have no states.";
      }
    }

    // Use number of states to weigh the downward transition
    if (m_host->getFlags().useDOSweighedDT){
      // The collision operator.
      for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
        double ei = gEne[i] ;
        double ni = gDOS[i] ;
        for ( j = i ; j < m_ncolloptrsize ; ++j ) {
          double ej = gEne[j];
          double nj = gDOS[j];
          // Transfer to lower Energy -
          double transferDown = exp(-alpha*(ej - ei)) * (ni/nj);
          (*m_egme)[i][j] = transferDown;

          // Transfer to higher Energy (via detailed balance) -
          double transferUp = exp(-(alpha + beta)*(ej - ei));
          (*m_egme)[j][i] = transferUp;
        }
      }
    }
    else{
      // The collision operator.
      for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
        double ei = gEne[i] ;
        double ni = gDOS[i] ;
        for ( j = i ; j < m_ncolloptrsize ; ++j ) {
          double ej = gEne[j];
          double nj = gDOS[j];
          // Transfer to lower Energy -
          double transferDown = exp(-alpha*(ej - ei)) ;
          (*m_egme)[i][j] = transferDown;

          // Transfer to higher Energy (via detailed balance) -
          double transferUp = exp(-alpha*(ej - ei)) * (nj/ni) * exp(-beta*(ej - ei)) ;
          (*m_egme)[j][i] = transferUp;
        }
      }
    }

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
      ctest << endl << "Collision operator column Sums" << endl << "{" << endl ;
      for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
        double columnSum(0.0) ;
        for ( j = 0 ; j < m_ncolloptrsize ; ++j ){
          columnSum += to_double((*m_egme)[j][i]) ;
        }
        ctest << columnSum << endl ;
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
  double gWellProperties::matrixElement(int eigveci, int eigvecj, std::vector<double> &k, int ndim)
  {
    // Calculate matrix element starting with the higher energy
    // elements first in order to preserve precision as much as possible.
    double sum = 0.0 ;
    for (int i = m_ncolloptrsize - 1 ; i >= 0 ; --i){
      sum +=  k[i]* to_double((*m_egvec)[i][eigveci]*(*m_egvec)[i][eigvecj]) ;
    }
    return sum ;
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

    if (locate + size > smsize) {
      cerr << "Error in the size of the system matrix.";
      exit(1) ;
    }

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

    if (!m_host->getFlags().doBasisSetMethod) {
      cerr << "Error: Contracted basis representation not requested.";
      exit(1) ;
    }

    // Find size of system matrix.

    int smsize = static_cast<int>(CollOptr->size()) ;

    // Check there is enough space in system matrix.

    if (locate + get_nbasis() > smsize) {
      cerr << "Error in the size of the reaction operator matrix in contracted basis representation.";
      exit(1) ;
    }

    // Copy collision operator eigenvalues to the diagonal elements indicated
    // by "locate" and multiply by the reduced collision frequencey.

    for (int i(0) ; i < get_nbasis() ; ++i) {
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
      m_pDistributionCalculator->calculateDistribution(gDOS, gEne, m_host->getEnv().beta, m_grainDist);
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
      }
      else{
        grainFrac.push_back(tempGrainFrac[i]);
      }
    }

    int grainFracSize = int(grainFrac.size());
    for (int i = 0; i < grainFracSize; ++i){
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

    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    const double firstPartition = exp(log(gDOS[0]) - m_host->getEnv().beta * gEne[0] + 10.0);
    tempGrnFrac.push_back(firstPartition);
    double prtfn(firstPartition);
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

  gStructure::gStructure(mesmer::Molecule *pMol) : m_MolecularWeight(-1)
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element
    setMass(ppPropList->XmlReadPropertyDouble("me:MW"));
  }

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS
    (gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, std::vector<double>& rctsCellDOS){
      std::vector<double> rct1CellDOS;
      std::vector<double> rct2CellDOS;
      pDOS1.getCellDensityOfStates(rct1CellDOS);
      pDOS2.getCellDensityOfStates(rct2CellDOS);
      std::vector<double> rotConsts;
      if (pDOS1.get_rotConsts(rotConsts) == -4){
        rctsCellDOS = rct2CellDOS;
      }
      else if(pDOS2.get_rotConsts(rotConsts) == -4){
        rctsCellDOS = rct1CellDOS;
      }
      else{
        FastLaplaceConvolution(rct1CellDOS, rct2CellDOS, rctsCellDOS);
      }
      return true;
  }


}//namespace
