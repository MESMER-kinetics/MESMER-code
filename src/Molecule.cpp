//
// Molecule.cpp
//
// Author: Struan Robertson
// Date:   5/Jan/2003
//
// This file contains the implementation of the Molecule class and its derived classes.
//
//-------------------------------------------------------------------------------------------
#include "Molecule.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  Molecule::Molecule(const MesmerEnv& Env):
    m_Env(Env),
    m_flag(0),
    m_ppPersist(NULL),
    m_Name(),
    m_Mass(0.0),
    m_Sigma(sigmaDefault),
    m_Epsilon(epsilonDefault),
    m_Mass_chk(-1),
    m_Sigma_chk(-1),
    m_Epsilon_chk(-1)
  {}

  Molecule::~Molecule()
  {
    if (m_Mass_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_Mass is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_Sigma_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_Sigma is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_Epsilon_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_Epsilon is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
  }

  ModelledMolecule::ModelledMolecule(const MesmerEnv& Env): Molecule(Env),
    m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(1.0),
    m_ZPE(0.0),
    m_SpinMultiplicity(1),
    m_grnZpe(0),
    m_pDensityOfStatesCalculator(NULL),
    m_RC_chk(-1),
    m_Sym_chk(-1),
    m_ZPE_chk(-1),
    m_SpinMultiplicity_chk(-1),
    m_VibFreq_chk(-1),
    m_eleExc(),
    m_VibFreq(),
    m_cellEne(),
    m_cellDOS(),
    m_grainEne(),
    m_grainDOS()
  {}

  ModelledMolecule::~ModelledMolecule()
  {
    if (m_RC_chk == 0){
      stringstream errorMsg;
      errorMsg << "Rotational constants are provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_Sym_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_Sym is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_ZPE_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_ZPE is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_SpinMultiplicity_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_SpinMultiplicity is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_VibFreq_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_VibFreq is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    // Free any memory assigned for calculating densities of states.
    if (m_grainDOS.size()) m_grainDOS.clear();
    if (m_grainEne.size()) m_grainEne.clear();
    if (m_cellDOS.size()) m_cellDOS.clear();
    if (m_cellEne.size()) m_cellEne.clear();
    if (m_VibFreq.size()) m_VibFreq.clear();
  }

  CollidingMolecule::CollidingMolecule(const MesmerEnv& Env) : ModelledMolecule(Env),
    m_DeltaEdown(0.0),
    m_collisionFrequency(0.0),
    m_ncolloptrsize(0),
    m_DeltaEdown_chk(-1),
    m_egme(NULL)
  {}

  CollidingMolecule::~CollidingMolecule()
  {
    if (m_DeltaEdown_chk == 0){
      stringstream errorMsg;
      errorMsg << "m_DeltaEdown is provided but not used in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    if (m_egme != NULL) delete m_egme ;
  }

  TransitionState::TransitionState(const MesmerEnv& Env) : ModelledMolecule(Env)
  {}

  SuperMolecule::SuperMolecule(const MesmerEnv& Env) : ModelledMolecule(Env),
    m_mol1(NULL),
    m_mol2(NULL)
  {}

    SuperMolecule::~SuperMolecule()
  {}

  /* Will need Clone() functions
  Molecule::Molecule(const Molecule& molecule) {
  // Copy constructor - define later SHR 23/Feb/2003
  }

  Molecule& Molecule::operator=(const Molecule&) {
  // Assignment operator - define later SHR 23/Feb/2003

  return *this ;
  }
  */

  //
  // Read the Molecular data from I/O stream.
  bool Molecule::InitializeMolecule(PersistPtr pp)
  {
    m_ppPersist = pp;
    const char* id= m_ppPersist->XmlReadValue("id");
    if (id) m_Name = id;
    if (m_Name.empty()) {
      stringstream errorMsg;
      errorMsg << "Molecular name is absent.\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      m_Name = "unknown";
      setFlag(true);
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:MW\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true); // later put a function to calculate the molecular weight if the user forgot to provide it.
    }
    else { istringstream idata(txt); double mass(0.); idata >> mass; setMass(mass);}


    if (getFlag()){
      stringstream errorMsg;
      errorMsg << "Error(s) while initializing molecule: " << m_Name;
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    return true;
  }

  void   Molecule::setMass(double value)           {
    m_Mass = value;
    m_Mass_chk = 0;
  } ;

  double Molecule::getMass()                       {
    if (m_Mass_chk >= 0){
      ++m_Mass_chk;
      return m_Mass ;
    }
    else{
      stringstream errorMsg;
      errorMsg << "m_Mass was not defined but requested in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1);
    }
  } ;

  void   Molecule::setSigma(double value)          {
    m_Sigma = value;
    m_Sigma_chk = 0;
  } ;

  double Molecule::getSigma()                      {
    if (m_Sigma_chk >= 0){
      ++m_Sigma_chk;
      return m_Sigma ;
    }
    else{
      stringstream errorMsg;
      errorMsg << "m_Sigma was not defined but requested in " << getName()
               << ". Default value " << sigmaDefault << " is used.\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      //exit(1);
      return m_Sigma ;
    }
  } ;

  void   Molecule::setEpsilon(double value)        {
    m_Epsilon = value;
    m_Epsilon_chk = 0;
  } ;

  double Molecule::getEpsilon()                    {
    if (m_Epsilon_chk >= 0){
      ++m_Epsilon_chk;
      return m_Epsilon ;
    }
    else{
      stringstream errorMsg;
      errorMsg << "m_Epsilon was not defined but requested in " << getName()
               << ". Default value " << epsilonDefault << " is used.\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      //exit(1);
      return m_Epsilon ;
    }
  } ;

  void   CollidingMolecule::setDeltaEdown(double value)        {
    m_DeltaEdown = value;
    m_DeltaEdown_chk = 0;
  } ;

  double CollidingMolecule::getDeltaEdown()                    {
    if (m_DeltaEdown_chk >= 0){
      ++m_DeltaEdown_chk;
      return m_DeltaEdown ;
    }
    else{
      stringstream errorMsg;
      errorMsg << "m_DeltaEdown was not defined but requested in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1);
    }
  } ;

  bool BathGasMolecule::InitializeMolecule(PersistPtr pp)
  {
    stringstream errorMsg;

    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for " << getName() << " before constructing BathGasMolecule.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt = ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:sigma.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true);
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt = ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:epsilon.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true);
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    if (getFlag()){
      stringstream errorMsg;
      errorMsg << "Error(s) while initializing: " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    return true;
  }


  bool ModelledMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    PersistPtr oldpp = pp;

    if(!Molecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule for " << getName() << " before constructing ModelledMolecule with errors.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    bool hasVibFreq = true; bool hasRotConst = true;
    txt= ppPropList->XmlReadProperty("me:vibFreqs");
    if(!txt){
      hasVibFreq = false;
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:vibFreqs. Maybe an atom or atomic ion.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      m_VibFreq_chk = -1;
      //setFlag(true); // it maybe an atom so not necessary to set this flag. Just produce warning.
    }
    else { istringstream idata(txt); double x; while (idata >> x) m_VibFreq.push_back(x); m_VibFreq_chk = 0;}

    txt= ppPropList->XmlReadProperty("me:rotConsts");
    if(!txt){
      hasRotConst = false;
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:rotConsts. Maybe an atom or atomic ion.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
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
      stringstream errorMsg;
      errorMsg << getName()
      << " has improper setting on vibrational frequencies or rotational constants. Check input file to remove this error.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true);
    }

    txt= ppPropList->XmlReadProperty("me:eletronicExcitation");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:eletronicExcitation for " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    else {
      istringstream idata(txt); double _iele = 0.; m_eleExc.clear();
      while (idata >> _iele) m_eleExc.push_back(_iele);
    }

    txt= ppPropList->XmlReadProperty("me:symmetryNumber");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:symmetryNumber. Default value " << m_Sym << " is used.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      m_Sym_chk = -1;
      //setFlag(true);
    }
    else { istringstream idata(txt); idata >> m_Sym; m_Sym_chk = 0;}

    txt= ppPropList->XmlReadProperty("me:ZPE");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:ZPE";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      m_ZPE_chk = -1;
    }
    else { istringstream idata(txt); idata >> m_ZPE ; m_ZPE_chk = 0;}

    // Determine the method of DOS calculation.
    const char* pDOSCMethodtxt = pp->XmlReadValue("me:DOSCMethod", false) ;
    if(pDOSCMethodtxt)
    {
      m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
      if(!m_pDensityOfStatesCalculator) // if the provided method cannot be found,
      {
        stringstream errorMsg;
        errorMsg << "Unknown method " << pDOSCMethodtxt
          << " for the calculation of DOS in " << getName()
          << ". Please check spelling error. Default method <Classical rotors> is used.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        pDOSCMethodtxt = "Classical rotors";
        m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
      }
    }
    else{ // if no method is provided.
        stringstream errorMsg;
        errorMsg << "No method for the calculation of DOS in " << getName()
          << " is provided. Default method <Classical rotors> is used.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      pDOSCMethodtxt = "Classical rotors"; // must exist
      m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
    }

    txt= ppPropList->XmlReadProperty("me:spinMultiplicity");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:spinMultiplicity in " << getName()
               << ". Default value "<< m_SpinMultiplicity << " is used.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    else
    {
      istringstream idata(txt);
      idata >> m_SpinMultiplicity;
      m_SpinMultiplicity_chk = 0;
    }

    if (getFlag()){
      stringstream errorMsg;
      errorMsg << "Error(s) while initializing: " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    return true;
  }

  bool CollidingMolecule::InitializeMolecule(PersistPtr pp)
  {
    stringstream errorMsg;
    PersistPtr oldpp = pp;

    //Read base class parameters first
    if(!ModelledMolecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for " << getName()
      << " before constructing CollidingMolecule.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      setFlag(true);
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:sigma in " << getName()
               << ". Default value " << sigmaDefault << " is used.\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      // sigma and epsilon are not always necessary.
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt= ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:epsilon in " << getName()
               << ". Default value " << epsilonDefault << " is used.\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      // sigma and epsilon are not always necessary.
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    txt= ppPropList->XmlReadProperty("me:deltaEDown");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "Cannot find argument me:deltaEDown in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      //setFlag(true);
      // deltaEDown is not always necessary. Hoever, it is not wise to provide a default value.
    }
    else { istringstream idata(txt); idata >> m_DeltaEdown; m_DeltaEdown_chk = 0;}

    if (getFlag()){
      stringstream errorMsg;
      errorMsg << "Error(s) while initializing: " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    return true;
  }

  bool SuperMolecule::InitializeMolecule(PersistPtr pp)
  {
    //the construction of SuperMolecule should always rely on the components, there is therefore no need to initialize it
    //apart from name and m_ppPersist

    if (pp) {
      setPersistentPointer(pp);

      const char* id = pp->XmlReadValue("id");
      if (id) setName(id);
      else{
        stringstream errorMsg;
        errorMsg << "Molecular name is absent. Default name <source> is used.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        string tempName = "source"; setName(tempName);
        //setFlag(true);
      }
      return true;
    }
    else{
      stringstream errorMsg;
      errorMsg << "Invalid PersistPtr.\n";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
  }

  //
  // Get cell density of states.
  //
  void ModelledMolecule::getCellDensityOfStates(vector<double> &cellDOS) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates() ;
    cellDOS.assign(m_cellDOS.begin(), m_cellDOS.end());
  }

  //
  // Get cell energies.
  //
  void ModelledMolecule::getCellEnergies(vector<double> &CellEne) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates() ;
    CellEne.assign(m_cellEne.begin(), m_cellEne.end());
  }

  //
  // Get grain density of states.
  //
  void ModelledMolecule::grnDensityOfStates(vector<double> &grainDOS) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates() ;
    grainDOS.assign(m_grainDOS.begin(), m_grainDOS.end());
  }

  //
  // Get grain energies.
  //
  void ModelledMolecule::grnEnergies(vector<double> &grainEne) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates() ;
    grainEne.assign(m_grainEne.begin(), m_grainEne.end());
  }

  //
  // Get Grain Boltzmann distribution.
  //
  void ModelledMolecule::grnBoltzDist(vector<double> &grainBoltzDist)
  {
    // If density of states have not already been calcualted then do so.
     if(0){stringstream errorMsg;
     errorMsg << "Before calcDensityOfStates(), Molecular name: " << getName();
     meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);}

    if (!m_cellDOS.size())
      calcDensityOfStates() ;

     if(0){stringstream errorMsg;
     errorMsg << "After calcDensityOfStates(), Molecular name: " << getName();
     meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);}

    int MaximumGrain = getEnv().MaxGrn ;
    double beta = getEnv().beta;

    // Calculate the Boltzmann distribution.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.

    double prtfn(0.0) ;
    for (int i = 0; i < MaximumGrain; ++i) {
        double tmp = log(m_grainDOS[i]) - beta*m_grainEne[i] + 10.0 ;
        tmp = exp(tmp) ;
        prtfn += tmp ;
        grainBoltzDist[i] = tmp ;
    }

    // Normalize the Boltzmann distribution.

    for (int i = 0; i < MaximumGrain; ++i) {
        grainBoltzDist[i] /= prtfn ;
    }
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
  // Get Grain canonical partition function.
  //
  double ModelledMolecule::grnCanPrtnFn() {
    // If density of states have not already been calculated then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates() ;

    double CanPrtnFn(0.0) ;

    // Calculate the ro-vibrational partition function based on the grain
    // densities of states, and not the molecular properties, for consistency.

    const int MaximumGrain = getEnv().MaxGrn ;
    const double beta = getEnv().beta;

    for (int i = 0; i < MaximumGrain; ++i) {
		if (m_grainDOS[i] > 0.0) 
			CanPrtnFn += exp( log(m_grainDOS[i]) - beta*m_grainEne[i] ) ;
    }

	// The following catches the case where the molecule is a single atom

	CanPrtnFn = max(CanPrtnFn, 1.0) ;

    // Electronic partition function.
    CanPrtnFn *= double(getSpinMultiplicity()) ;

    {stringstream errorMsg;
    errorMsg << "CanPrtnFn = " << CanPrtnFn << ", molecular name: " << getName();
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

    // Translational partition function.
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
      {stringstream errorMsg;
      errorMsg << "Rotational constants were not defined but requested in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);}
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

  //-------------------------------------------------------------------------------------------
  //
  // Methods related to the Master Equation.

  //
  // Initialize the Collision Operator.
  //
  bool CollidingMolecule::initCollisionOperator(double beta, Molecule *pBathGasMolecule)
  {
    // If density of states have not already been calcualted then do so.

    if (!m_cellDOS.size())
      calcDensityOfStates() ;

    // Calculate the collision frequency.

    m_collisionFrequency = collisionFrequency(beta, getEnv().conc, pBathGasMolecule) ;

    // Calculate the collision operator.

    if (!collisionOperator(beta)){
      std::stringstream errorMsg;
      errorMsg << "Failed building collision operator";
      mesmer::meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obError);
      return false;
    }
    return true;
  }

  //
  // Calculate collision operator
  //
  bool CollidingMolecule::collisionOperator(double beta)
  {
    if (test_rotConsts() < 0) return true;
    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //
    int i, j;
    int MaximumGrain = getEnv().MaxGrn;

    if(!m_DeltaEdown){
      stringstream errorMsg;
      errorMsg << "me:deltaEDown is necessary for " << getName() << ". Correct input file to remove this error.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    double alpha = 1.0/m_DeltaEdown ;
    //
    // Allocate memmory.
    //
    if (m_egme)                                 // Delete any existing matrix.
        delete m_egme ;

    m_egme = new dMatrix(MaximumGrain) ;              // Collision operator matrix.

    //
    // Initialisation and error checking.
    //
    for ( i = 0 ; i < MaximumGrain ; ++i ) {
      if (m_grainDOS[i] <= 0.0) {
        stringstream errorMsg;
        errorMsg << "Data indicates that grain " << i << " of the current colliding molecule has no states.";
        meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
    }
    if(0){stringstream errorMsg;
    errorMsg << "alpha = " << alpha << ", m_DeltaEdown = " << m_DeltaEdown;
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}
    //
    // The collision operator.
    //
    for ( i = 0 ; i < MaximumGrain ; ++i ) {
      double ei = m_grainEne[i] ;
      double ni = m_grainDOS[i] ;
      for ( j = i ; j < MaximumGrain ; ++j ) {
        //
        // Transfer to lower Energy -
        //
        (*m_egme)[i][j] = exp(-alpha*(m_grainEne[j] - ei)) ;
        //
        // Transfer to higher Energy (via detailed balance) -
        //
        (*m_egme)[j][i] = (*m_egme)[i][j] * (m_grainDOS[j]/ni)*exp(-beta*(m_grainEne[j] - ei)) ;
      }
    }
      //-----------------------------------------------------CHECKBLOCK
      if(1 && MaximumGrain < 11){stringstream errorMsg;
      errorMsg << "Collison operator:\n";
      for (int i = 0; i < MaximumGrain; ++i){
        for (int j = 0; j < MaximumGrain; ++j){
          formatFloat(errorMsg, to_double((*m_egme)[i][j]), 6,  14);
        }
        errorMsg << endl;
      }
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}
      //-----------------------------------------------------CHECKBLOCK

    //Normalisation
    m_egme->normalize();

    //account for collisional loss by subrtacting unity from the leading diagonal.
    for ( i = 0 ; i < MaximumGrain ; ++i ) (*m_egme)[i][i] -= 1.0 ;

      //-----------------------------------------------------CHECKBLOCK
      if(1 && MaximumGrain < 11){stringstream errorMsg;
      errorMsg << "Collison operator:\n";
      for (int i = 0; i < MaximumGrain; ++i){
        for (int j = 0; j < MaximumGrain; ++j){
          formatFloat(errorMsg, to_double((*m_egme)[i][j]), 6,  14);
        }
        errorMsg << endl;
      }
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}
      //-----------------------------------------------------CHECKBLOCK

    if (getEnv().collisionOCSEnabled){
      ctest << endl << "Collision operator column Sums" << endl << "{" << endl ;
      for ( i = 0 ; i < MaximumGrain ; ++i ) {
        double columnSum(0.0) ;
        for ( j = 0 ; j < MaximumGrain ; ++j ){
          columnSum += to_double((*m_egme)[j][i]) ;
        }
        ctest << columnSum << endl ;
      }
      ctest << "}" << endl;
    }

    //
    // Determine the equilibrium vector for symmetrization. The work
    // contains square root of the Boltzman distribution.
    //

    vector<double> work(MaximumGrain,0.0) ; // Work space.
    for ( i = 0 ; i < MaximumGrain ; ++i ) {
      double ei = log(m_grainDOS[i]) - beta*m_grainEne[i] + 10.0 ;
      ei = exp(ei) ;
      work[i] = sqrt(ei) ;
    }
    //
    // Symmetrization of the collision matrix.
    //
    for ( i = 1 ; i < MaximumGrain ; ++i ) {
      for ( j = 0 ; j < i ; ++j ) {
        (*m_egme)[j][i] *= (work[i]/work[j]) ;
        (*m_egme)[i][j]  = (*m_egme)[j][i] ;
      }
    }
    return true;
  }

  //
  // Calculate collision frequency.
  //
  double CollidingMolecule::collisionFrequency(double beta, const double conc, Molecule *pBathGasMolecule)
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
    double bthMass    = pBathGasMolecule->getMass();

    double bthSigma   = pBathGasMolecule->getSigma();
    if (!bthSigma){
      stringstream errorMsg;
      errorMsg << "me:sigma is necessary for " << pBathGasMolecule->getName() << ". Correct input file to remove this error.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    double bthEpsilon = pBathGasMolecule->getEpsilon();
    if (!bthEpsilon){
      stringstream errorMsg;
      errorMsg << "me:epsilon is necessary for " << pBathGasMolecule->getName() << ". Correct input file to remove this error.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

    double mu   = amu * getMass() * bthMass/(getMass() + bthMass) ;
    double eam  = sqrt(getEpsilon() * bthEpsilon) ;
    double sam  = (getSigma() + bthSigma) * 0.5;
    double tstr = 1. / (eam * beta);

    if(0){stringstream errorMsg;
    errorMsg << "mu = " << mu << ", eam = " << eam << ", sam = " << sam << ", tstr = " << tstr << ", beta = " << beta;
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

    // Calculate collision integral.
    double collFrq = A * exp(-log(tstr) * B) + C * exp(-D * tstr) + E * exp(-F * tstr) ;

    // Calculate molecular collision frequency.
    collFrq *= M_PI * sam * sam * 1.0e-20 * sqrt(8. * 1.381e-23 * temp/(M_PI * mu)) ;

    // Calculate overall collision frequency.
    collFrq *= conc * 1.0e+06 ;

    if(1){stringstream errorMsg;
    errorMsg << "Collision frequency of " << getName() << " is " << collFrq;
    meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);}

    return collFrq;
  }

  //
  // Calculate a reaction matrix element.
  //
  double CollidingMolecule::matrixElement(int eigveci, int eigvecj, vector<double> &k, int ndim)
  {
    double sum = 0.0 ;
    for (int i = 0 ; i < ndim ; ++i){
      sum +=  k[i]* to_double((*m_egme)[i][eigveci]*(*m_egme)[i][eigvecj]) ;

    }
    return sum ;
  }

  //
  // Copy collision operator to diagonal block of system matrix.
  //
  void CollidingMolecule::copyCollisionOperator(dMatrix *CollOptr,
                                                const int size,
                                                const int locate,
                                                const double RducdOmega) const
  {
    // Find size of system matrix.

    int smsize = static_cast<int>(CollOptr->size()) ;
    int MaximumGrain = getEnv().MaxGrn;

    // Check there is enough space in system matrix.

    if (locate + size > smsize) {
      stringstream errorMsg;
      errorMsg << "Error in the size of the system matrix.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
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

  //-------------------------------------------------------------------------------------------
  //
  // Private methods.

  //
  // Calculate the rovibrational density of states.
  //
  bool ModelledMolecule::calcDensityOfStates()
  {
    if (!m_pDensityOfStatesCalculator->countCellDOS(this)){
      return false;
    }
    calcGrainAverages(); testDensityOfStates() ;
    return true;
  }

  //
  // Test the rovibrational density of states for ModelledMolecule.
  //
  void ModelledMolecule::testDensityOfStates()
  {
    ctest << endl << "Test density of states for: " << getName() << endl << endl ;
    int MaximumGrain = getEnv().MaxGrn;
    int MaximumCell  = getEnv().MaxCell;

    string comment("Partition function calculation at various temperatures.\n qtot : partition function as a product of quantum mechanical partition functions for vibrations (1-D harmonic oscillator) and classifical partition functions for rotations\n sumc : (user calculated) cell based partition function \n sumg : (user calculated) grain based partition function ");

    PersistPtr ppList = getPersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment );

    ctest << "      T           qtot           sumc           sumg\n";

    //loop through predefined test temperatures
    for ( int n = 0 ; n < 29 ; ++n ) {
      double temp = 100.0*static_cast<double>(n + 2) ;
      double beta = 1.0/(boltzmann_RCpK*temp) ;

      // Calculate partition functions based on cells.

      double sumc  = 0.0 ;
      for ( int i = 0 ; i < MaximumCell ; ++i ) {
        sumc += m_cellDOS[i]*exp(-beta*m_cellEne[i]) ;
      }

      // Calculate partition functions based on grains.

      double sumg  = 0.0 ;
      for ( int i = 0 ; i < MaximumGrain ; ++i ) {
        sumg += m_grainDOS[i]*exp(-beta*m_grainEne[i]) ;
      }

      // Calculate partition functions using analytical formula (treat vibrations classically).

      double qtot = 1.0 ;

      vector<double> rotConst; int rotorType;

      if (!dynamic_cast<SuperMolecule *>(this)) {
        rotorType = get_rotConsts(rotConst);
      }
      else{rotorType = -4;}

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
          qtot *= (rotConst[0] / (get_Sym()*beta)) ;
          break;
        default:
          qtot = 0.;
      }
      formatFloat(ctest, temp,  6,  7) ;
      formatFloat(ctest, qtot,  6, 15) ;
      formatFloat(ctest, sumc,  6, 15) ;
      formatFloat(ctest, sumg,  6, 15) ;
      ctest << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
      ppItem->XmlWriteValueElement("me:T",    temp, 6);
      ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
      ppItem->XmlWriteValueElement("me:sumc", sumc, 6);
      ppItem->XmlWriteValueElement("me:sumg", sumg, 6);
    }
    if (getEnv().cellDOSEnabled){
      ctest << endl << "Cell density of states of " << getName() << endl << "{" << endl;
      for (int i = 0; i < MaximumCell; ++i){
        formatFloat(ctest, m_cellEne[i],  6,  15) ;
        formatFloat(ctest, m_cellDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }

    if (getEnv().grainDOSEnabled){
      ctest << endl << "Grain density of states of " << getName() << endl << "{" << endl;
      for (int i = 0; i < MaximumGrain; ++i){
        formatFloat(ctest, m_grainEne[i],  6,  15) ;
        formatFloat(ctest, m_grainDOS[i],  6,  15) ;
        ctest << endl ;
      }
      ctest << "}" << endl;
    }
  }

  //
  // Calculate the average grain energy and then number of states per grain.
  //
  void ModelledMolecule::calcGrainAverages()
  {
    int MaximumGrain = getEnv().MaxGrn;
    m_grainEne.resize(MaximumGrain, 0.) ;
    m_grainDOS.resize(MaximumGrain, 0.) ;

    //    int igsz = MAXCELL/MAXGRN ;

    // Check that there are enough cells.

    if (getEnv().GrainSize < 1) {
      stringstream errorMsg;
      errorMsg << "The number of Cells is insufficient to produce requested number of Grains.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;
    int grainsize = static_cast<int>(getEnv().GrainSize);
    for (int i = 0 ; i < MaximumGrain ; ++i ) {

      int idx3 = idx1 ;

      // Calculate the number of states in a grain.
      double gNOS = 0.0 ;
      for (int j = 0 ; j < grainsize ; ++j, ++idx1 ){
        gNOS += m_cellDOS[idx1] ;
      }

      // Calculate average energy of the grain if it contains sum states.
      if ( gNOS > 0.0 ) {
        double gSE = 0.0 ; // grain sum of state energy
        for (int j = 0 ; j < grainsize ; ++j, ++idx3 ){
          gSE += m_cellEne[idx3] * m_cellDOS[idx3] ;
        }
        m_grainDOS[idx2] = gNOS ;
        m_grainEne[idx2] = gSE/gNOS ;
      }

      idx2++ ;
    }

    // Issue warning if number of grains produced is less that requested.

    if ( idx2 != MaximumGrain ) {
      stringstream errorMsg;
      errorMsg << "Number of grains produced is not equal to that requested" << endl
               << "Number of grains requested: " << MaximumGrain << endl
               << "Number of grains produced : " << idx2 << " in " << getName();
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    else{
//      stringstream errorMsg;
//      errorMsg << "Number of grains requested: " << MaximumGrain << endl
//               << "Number of grains produced : " << idx2 << " in " << getName();
//      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

  }

  double ModelledMolecule::get_zpe() {
    if (m_ZPE_chk == -1){
      stringstream errorMsg;
      errorMsg << "m_ZPE was not defined but requested in " << getName() << ". Default value " << m_ZPE << " is given.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      --m_ZPE_chk;
      return m_ZPE;
    }
    else if (m_ZPE_chk < -1){
      --m_ZPE_chk;
      return m_ZPE;
    }
    ++m_ZPE_chk;
    return m_ZPE ;
  }

  void ModelledMolecule::set_zpe(double value) {
    m_ZPE = value;
    ++m_ZPE_chk;
  }

  double ModelledMolecule::get_Sym(void){
    if (m_Sym_chk == -1){
      stringstream errorMsg;
      errorMsg << "m_Sym was not defined but requested in " << getName() << ". Default value " << m_Sym << " is given.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
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
    vibFreq.clear();
    if (m_VibFreq_chk >=0){
      vibFreq.assign(m_VibFreq.begin(), m_VibFreq.end());
      ++m_VibFreq_chk;
    }
  }

  void ModelledMolecule::set_VibFreq(std::vector<double>& vibFreq, bool addition){
    if (addition){
      for (std::vector<double>::size_type i = 0; i < vibFreq.size(); ++i) m_VibFreq.push_back(vibFreq[i]);
    }
    else{
      m_VibFreq.clear();
      m_VibFreq.assign(vibFreq.begin(), vibFreq.end());
    }
  }

  int ModelledMolecule::getSpinMultiplicity()              {
    if (m_SpinMultiplicity_chk >= 0){
      ++m_SpinMultiplicity_chk;
      return m_SpinMultiplicity ;
    }
    else{
      stringstream errorMsg;
      errorMsg << "m_SpinMultiplicity was not defined but requested in " << getName() << ". Default value " << m_SpinMultiplicity << " is given.";
      meErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return m_SpinMultiplicity;
    }
  }

  void   ModelledMolecule::setSpinMultiplicity(int value){
    m_SpinMultiplicity = value;
  }

  // set composing member of the SuperMolecule, also copy necessary properties
  void SuperMolecule::setMembers(ModelledMolecule* mol1p, ModelledMolecule* mol2p){
    m_mol1 = mol1p; std::vector<double> vfMol1; m_mol1->get_VibFreq(vfMol1);
    m_mol2 = mol2p; std::vector<double> vfMol2; m_mol2->get_VibFreq(vfMol2);
    set_VibFreq(vfMol1);
    set_VibFreq(vfMol2, true);
    set_zpe(m_mol1->get_zpe() + m_mol2->get_zpe());
    setSpinMultiplicity(m_mol1->getSpinMultiplicity() + m_mol2->getSpinMultiplicity());
    set_DensityOfStatesCalculator(m_mol1->get_DensityOfStatesCalculator());
  }

}//namespace
