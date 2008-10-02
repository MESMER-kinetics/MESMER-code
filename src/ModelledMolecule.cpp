//
// ModelledMolecule.cpp
//
// Author: Struan Robertson
//-------------------------------------------------------------------------------------------
#include "Molecule.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{

  //
  //Constructor
  //
  ModelledMolecule::ModelledMolecule(const MesmerEnv& Env, MesmerFlags& Flags): Molecule(Env, Flags),
    m_Mass(0.0),
    m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(1.0),
    m_ZPE(0.0),
    m_scaleFactor(1.0),
    m_SpinMultiplicity(1),
    m_initPopulation(0.0),
    m_eqFraction(0.0),
    m_pDensityOfStatesCalculator(NULL),
    m_Mass_chk(-1),
    m_RC_chk(-1),
    m_Sym_chk(-1),
    m_ZPE_chk(-1),
    m_scaleFactor_chk(-1),
    m_SpinMultiplicity_chk(-1),
    m_VibFreq_chk(-1),
    m_eleExc(),
    m_VibFreq(),
    m_cellDOS(),
    m_grainEne(),
    m_grainDOS()
  {}

  ModelledMolecule::~ModelledMolecule()
  {
    if (m_Mass_chk == 0) cinfo << "m_Mass is provided but not used in " << getName() << endl;
    if (m_RC_chk == 0) cinfo << "Rotational constants are provided but not used in " << getName() << endl;
    if (m_Sym_chk == 0) cinfo << "m_Sym is provided but not used in " << getName() << endl;
    if (m_ZPE_chk == 0) cinfo << "m_ZPE is provided but not used in " << getName() << endl;
    if (m_scaleFactor_chk == 0) cinfo << "m_scaleFactor is provided but not used in " << getName() << endl;
    if (m_SpinMultiplicity_chk == 0) cinfo << "m_SpinMultiplicity is provided but not used in " << getName() << endl;
    if (m_VibFreq_chk == 0) cinfo << "m_VibFreq is provided but not used in " << getName() << endl;

    // Free any memory assigned for calculating densities of states. (must be in reverse order)
    if (m_grainDOS.size()) m_grainDOS.clear();
    if (m_grainEne.size()) m_grainEne.clear();
    if (m_cellDOS.size()) m_cellDOS.clear();
    if (m_VibFreq.size()) m_VibFreq.clear();
    if (m_eleExc.size()) m_eleExc.clear();
  }

  //
  //Initialization
  //
  bool ModelledMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    PersistPtr oldpp = pp;

    if(!Molecule::InitializeMolecule(pp)){
      cerr << "InitializeMolecule for " << getName() << " before constructing ModelledMolecule with errors.";
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      cerr << "Cannot find argument me:MW in " << getName();
      setFlag(true); // later put a function to calculate the molecular weight if the user forgot to provide it.
    }
    else { istringstream idata(txt); double mass(0.); idata >> mass; setMass(mass);}

    bool hasVibFreq = true; bool hasRotConst = true;
    txt= ppPropList->XmlReadProperty("me:vibFreqs");
    if(!txt){
      hasVibFreq = false;
      cinfo << "Cannot find argument me:vibFreqs. Maybe an atom or atomic ion." << endl;
      m_VibFreq_chk = -1;
      //setFlag(true); // it maybe an atom so not necessary to set this flag. Just produce warning.
    }
    else { istringstream idata(txt); double x; while (idata >> x) m_VibFreq.push_back(x); m_VibFreq_chk = 0;}

    txt= ppPropList->XmlReadProperty("me:rotConsts");
    if(!txt){
      hasRotConst = false;
      cinfo << "Cannot find argument me:rotConsts. Maybe an atom or atomic ion." << endl;
      m_RC_chk = -1;
      //setFlag(true); // it maybe an atom so not necessary to set this flag. Just produce warning.
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
      std::sort(rCnst.begin(), rCnst.end());
      m_RotCstA = rCnst[2];
      m_RotCstB = rCnst[1];
      m_RotCstC = rCnst[0];
      m_RC_chk = 0;
    }

    if (hasVibFreq != hasRotConst){
      cerr << getName()
        << " has improper setting on vibrational frequencies or rotational constants. Check input file to remove this error.";
      setFlag(true);
    }

    txt= ppPropList->XmlReadProperty("me:eletronicExcitation");
    if(!txt){
      cinfo << "Cannot find argument me:eletronicExcitation for " << getName() << endl;
    }
    else {
      istringstream idata(txt); double _iele = 0.; m_eleExc.clear();
      while (idata >> _iele) m_eleExc.push_back(_iele);
    }

    txt= ppPropList->XmlReadProperty("me:symmetryNumber");
    if(!txt){
      cinfo << "Cannot find argument me:symmetryNumber. Default value " << m_Sym << " is used." << endl;
      m_Sym_chk = -1;
      //setFlag(true);
    }
    else { istringstream idata(txt); idata >> m_Sym; m_Sym_chk = 0;}

    txt = ppPropList->XmlReadProperty("me:ZPE");
    if(!txt){
      cinfo << "Cannot find argument me:ZPE" << endl;
      m_ZPE_chk = -1;
    }
    else {
      istringstream idata(txt);
      double tempzpe = 0.0;
      idata >> tempzpe;
      txt= ppPropList->XmlReadPropertyAttribute("me:ZPE", "units");
      string unitsInput;
      if (txt){
        unitsInput = txt;
      }
      else{
        unitsInput = "kJ/mol";
      }
      const char* pLowertxt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "lower");
      const char* pUppertxt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "upper");
      const char* pStepStxt = ppPropList->XmlReadPropertyAttribute("me:ZPE", "stepsize");
      double value(getConvertedEnergy(unitsInput, tempzpe));
      if (pLowertxt && pUppertxt){
        double tempLV(0.0), tempUV(0.0), tempSS(0.0);
        stringstream s3(pLowertxt), s4(pUppertxt), s5(pStepStxt); s3 >> tempLV; s4 >> tempUV; s5 >> tempSS;
        double valueL(getConvertedEnergy(unitsInput, tempLV)),
          valueU(getConvertedEnergy(unitsInput, tempUV)),
          stepsize(getConvertedEnergy(unitsInput, tempSS));
        set_zpe(valueL, valueU, stepsize);
      }
      else{
        set_zpe(value);
      }
      m_ZPE_chk = 0;
    }

    // The reason why me:frequenciesScaleFactor stands out to be a separate property in the propertyList is that
    // this value is not usually necessary. The default value is 1.0 and it is usually the case.
    txt= ppPropList->XmlReadProperty("me:frequenciesScaleFactor");
    if(!txt){
      cinfo << "Cannot find argument me:frequenciesScaleFactor" << endl;
      m_scaleFactor_chk = -1;
    }
    else { istringstream idata(txt); idata >> m_scaleFactor ; m_scaleFactor_chk = 0;}

    // Determine the method of DOS calculation.
    const char* pDOSCMethodtxt = pp->XmlReadValue("me:DOSCMethod", false) ;
    if(pDOSCMethodtxt)
    {
      m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
      if(!m_pDensityOfStatesCalculator) // if the provided method cannot be found,
      {
        cinfo << "Unknown method " << pDOSCMethodtxt
          << " for the calculation of DOS in " << getName()
          << ". Please check spelling error. Default method <Classical rotors> is used." << endl;
        pDOSCMethodtxt = "Classical rotors";
        m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
      }
    }
    else{ // if no method is provided.
      cinfo << "No method for the calculation of DOS in " << getName()
        << " is provided. Default method <Classical rotors> is used." << endl;
      pDOSCMethodtxt = "Classical rotors"; // must exist
      m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
    }

    txt= ppPropList->XmlReadProperty("me:spinMultiplicity");
    if(!txt){
      cinfo << "Cannot find argument me:spinMultiplicity in " << getName()
        << ". Default value "<< m_SpinMultiplicity << " is used." << endl;
    }
    else
    {
      istringstream idata(txt);
      idata >> m_SpinMultiplicity;
      m_SpinMultiplicity_chk = 0;
    }

    if (getErrorFlag()){
      cerr << "Error(s) while initializing: " << getName();
      return false;
    }

    return true;
  }

  //
  // Get cell density of states.
  //
  void ModelledMolecule::getCellDensityOfStates(vector<double> &cellDOS, int startingCell) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";
    if (startingCell == 0)
      cellDOS = m_cellDOS;
    else{
      int MaximumCell = getEnv().MaxCell;
      for (int i(startingCell); i < MaximumCell; ++i){
        cellDOS.push_back(m_cellDOS[i]);
      }
    }
  }

  //
  // Get grain density of states.
  //
  void ModelledMolecule::getGrainDensityOfStates(vector<double> &grainDOS, const int startGrnIdx, const int ignoreCellNumber) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates()){
      cerr << "Failed calculating DOS for " << getName();
      exit(1);
    }
    if (ignoreCellNumber == 0){ // If there is no cells ignored in this grain, the grain DOS dose not need to be recalculated.
      grainDOS = m_grainDOS;
    }
    else{ // Some cells are ignored in this grain, as they do not occur in this part of reaction.
      // first deal with the first grain.
      const int MaximumCell = getEnv().MaxCell;
      const int gsz = getEnv().GrainSize;
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
  void ModelledMolecule::getGrainEnergies(vector<double> &grainEne) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";
    grainEne = m_grainEne;
  }

  //
  // Get Electronic excitations
  //
  void ModelledMolecule::getEleExcitation(vector<double> &elecExci){
    elecExci.clear();
    for (vector<double>::size_type i = 0; i < m_eleExc.size(); ++i){
      elecExci.push_back(m_eleExc[i]);
    }
  }

  //
  // Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
  //
  double ModelledMolecule::rovibronicGrnCanPrtnFn() {
    // If density of states have not already been calculated then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";

    // Calculate the rovibronic partition function based on the grain DOS
    // The following catches the case where the molecule is a single atom
    double CanPrtnFn = max(canonicalPartitionFunction(m_grainDOS, m_grainEne, getEnv().beta), 1.0) ;
    if (CanPrtnFn == 1.0){
      // Electronic partition function for atom is accounted here.
      CanPrtnFn = double(getSpinMultiplicity()) ;
    }

    return CanPrtnFn ;
  }

  int ModelledMolecule::test_rotConsts()
  {
    std::vector<double> mmtsInt;
    return get_rotConsts(mmtsInt);
  }

  int ModelledMolecule::get_rotConsts(std::vector<double> &mmtsInt)
  {
    //if (m_RC_chk = -1){ // replace the line below by this line _2007_12_07__16_01_51_ you will encounter a problem somewhere else
    if (m_RC_chk == -1){
      cinfo << "Rotational constants were not defined but requested in " << getName() << endl;
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
  bool ModelledMolecule::calcDensityOfStates()
  {
    const bool recalc(needReCalculateDOS());
    const bool vectorSizeConstant(m_cellDOS.size() == static_cast<unsigned int>(getEnv().MaxCell));
    const size_t sizeOfVector(m_cellDOS.size());

    if (sizeOfVector && vectorSizeConstant && !recalc)
      return true;
    if (!get_DensityOfStatesCalculator()->countCellDOS(this)){
      return false;
    }

    std::vector<double> shiftedCellDOS;
    std::vector<double> shiftedCellEne;
    const int MaximumCell = getEnv().MaxCell;
    const int cellOffset = get_cellOffset();
    std::vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);
    shiftCells(MaximumCell, cellOffset, m_cellDOS, cellEne, shiftedCellDOS, shiftedCellEne);

    calcGrainAverages(getEnv().MaxGrn, getEnv().GrainSize, shiftedCellDOS, shiftedCellEne, m_grainDOS, m_grainEne, getName());

    testDensityOfStates();

    recalculateDOScompleted();
    return true;
  }

  // Calculate classical energy
  double ModelledMolecule::getClassicalEnergy(){
    //Basically use the frequencies to calculate the contribution of ZPE from harmonic oscillators approximation
    double ZC = 0.0;
    for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
      ZC += m_VibFreq[i] / 2.0;
    return get_zpe() - ZC;
  }


  //
  // Test the rovibronic density of states for ModelledMolecule.
  //
  void ModelledMolecule::testDensityOfStates()
  {
    const int MaximumGrain = getEnv().MaxGrn;
    const int MaximumCell  = getEnv().MaxCell;
    std::vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);

    string comment("Rovibronic partition function calculation at various temperatures. qtot : partition function as a product of quantum mechanical partition functions for vibrations (1-D harmonic oscillator) and classifical partition functions for rotations.  sumc : (user calculated) cell based partition function. sumg : (user calculated) grain based partition function ");

    PersistPtr ppList = getPersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment );

    if (getFlags().testDOSEnabled) ctest << endl << "Test rovibronic density of states for: " << getName() << "\n{\n";
    if (getFlags().testDOSEnabled) ctest << "      T           qtot           sumc           sumg\n";

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
      // times the scale factor
      for (vector<double>::size_type i = 0; i < vibFreq.size(); ++i){
        vibFreq[i] *= get_scaleFactor();
      }

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

      if (getFlags().testDOSEnabled) formatFloat(ctest, temp,  6,  7) ;
      if (getFlags().testDOSEnabled) formatFloat(ctest, qtot,  6, 15) ;
      if (getFlags().testDOSEnabled) formatFloat(ctest, cellCanPrtnFn,  6, 15) ;
      if (getFlags().testDOSEnabled) formatFloat(ctest, grainCanPrtnFn,  6, 15) ;
      if (getFlags().testDOSEnabled) ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
      ppItem->XmlWriteValueElement("me:T",    temp, 6);
      ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
      ppItem->XmlWriteValueElement("me:sumc", cellCanPrtnFn, 6);
      ppItem->XmlWriteValueElement("me:sumg", grainCanPrtnFn, 6);
    }
    if (getFlags().testDOSEnabled) ctest << "}" << endl;

    if (getFlags().cellDOSEnabled){
      ctest << endl << "Cell rovibronic density of states of " << getName() << endl << "{" << endl;
      for (int i = 0; i < MaximumCell; ++i){
        formatFloat(ctest, cellEne[i],  6,  15) ;
        formatFloat(ctest, m_cellDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }

    if (getFlags().grainDOSEnabled){
      ctest << endl << "Grain rovibronic density of states of " << getName() << endl << "{" << endl;
      for (int i = 0; i < MaximumGrain; ++i){
        formatFloat(ctest, m_grainEne[i],  6,  15) ;
        formatFloat(ctest, m_grainDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }
  }

  void   ModelledMolecule::setMass(double value)           {
    m_Mass = value;
    m_Mass_chk = 0;
  } ;

  double ModelledMolecule::getMass()                       {
    if (m_Mass_chk >= 0){
      ++m_Mass_chk;
      return m_Mass ;
    }
    else{
      cerr << "m_Mass was not defined but requested in " << getName();
      exit(1);
    }
  } ;

  double ModelledMolecule::get_zpe() {
    if (m_ZPE_chk == -1){
      cinfo << "m_ZPE was not defined but requested in " << getName() << ". Default value " << m_ZPE.get_value() << " is given." << endl;
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

  double ModelledMolecule::get_scaleFactor() {
    if (m_scaleFactor_chk == -1){
      cinfo << "m_scaleFactor was not defined but requested in " << getName() << ". Default value " << m_scaleFactor << " is given." << endl;
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

  void ModelledMolecule::set_scaleFactor(double value) {
    m_scaleFactor = value;
    m_scaleFactor_chk = 0;
  }

  double ModelledMolecule::get_Sym(void){
    if (m_Sym_chk == -1){
      cinfo << "m_Sym was not defined but requested in " << getName() << ". Default value " << m_Sym << " is given." << endl;
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

  void ModelledMolecule::get_VibFreq(std::vector<double>& vibFreq){
    if (m_VibFreq_chk >=0){
      for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
        vibFreq.push_back(m_VibFreq[i]);
      ++m_VibFreq_chk;
    }
  }

  int ModelledMolecule::getSpinMultiplicity(){
    if (m_SpinMultiplicity_chk >= 0){
      ++m_SpinMultiplicity_chk;
      return m_SpinMultiplicity ;
    }
    else{
      cinfo << "m_SpinMultiplicity was not defined but requested in " << getName() << ". Default value " << m_SpinMultiplicity << " is given." << endl;
      return m_SpinMultiplicity;
    }
  }

  void   ModelledMolecule::setSpinMultiplicity(int value){
    m_SpinMultiplicity = value;
  }


}//namespace

