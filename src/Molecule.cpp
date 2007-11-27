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
  Molecule::Molecule():
    m_flag(0),
    m_ppPersist(NULL),
    m_Name(),
    m_Mass(0.0),
    m_Sigma(0.0),
    m_Epsilon(0.0)
  {}

  ModelledMolecule::ModelledMolecule():
    m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(0.0),
    m_ZPE(0.0),
    m_SpinMultiplicity(1),
    m_grnZpe(0),
    m_pDensityOfStatesCalculator(NULL),
    m_VibFreq(),
    m_cellEne(),
    m_cellDOS(),
    m_grainEne(),
    m_grainDOS()
  {}

  ModelledMolecule::~ModelledMolecule()
  {
    // Free any memory assigned for calculating densities of states.
    if (m_cellDOS.size()) m_grainDOS.clear();
    if (m_cellEne.size()) m_grainEne.clear();
    if (m_cellDOS.size()) m_cellDOS.clear();
    if (m_cellEne.size()) m_cellEne.clear();
    if (m_cellDOS.size()) m_VibFreq.clear();
  }

  CollidingMolecule::CollidingMolecule():
    m_DeltaEdown(0.0),
    m_collisionFrequency(0.0),
    m_ncolloptrsize(0),
    m_egme(NULL)
  {}

  CollidingMolecule::~CollidingMolecule()
  {
    if (m_egme != NULL) delete m_egme ;
  }

  TransitionState::TransitionState()
  {}

  SuperMolecule::SuperMolecule():
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
    stringstream errorMsg;
    m_ppPersist = pp;
    const char* id= m_ppPersist->XmlReadValue("id");
    if (id) m_Name = id;
    if (m_Name.empty()) { errorMsg << "Molecular name is absent.\n"; m_Name = "unknown"; setFlag(true); }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      errorMsg << "Molecule::Cannot find argument me:MW\n";
      //setFlag(true);
    }
    else { istringstream idata(txt); double mass(0.); idata >> mass; setMass(mass);}


    if (getFlag()){
      errorMsg << "Error(s) while initializing molecule: " << m_Name;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }
    return true;
  }

  bool BathGasMolecule::InitializeMolecule(PersistPtr pp)
  {
    stringstream errorMsg;

    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for Molecule " << getName() << " before constructing BathGasMolecule.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      errorMsg << "Molecule::Cannot find argument me:sigma.\n";
      setSigma(sigmaDefault);
      setFlag(true);
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt= ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      errorMsg << "Molecule::Cannot find argument me:epsilon.\n";
      setEpsilon(epsilonDefault);
      setFlag(true);
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    if (getFlag()){
      errorMsg << "Error(s) while initializing molecule: " << getName();
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }

    return true;
  }


  bool ModelledMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    stringstream errorMsg;
    PersistPtr oldpp = pp;

    if(!Molecule::InitializeMolecule(pp)){
      errorMsg << "InitializeMolecule for Molecule " << getName() << " before constructing ModelledMolecule with errors.";
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:vibFreqs");
    if(!txt){
      errorMsg << "ModelledMolecule::Cannot find argument me:vibFreqs";
      //setFlag(true); // it maybe an atom so not necessary to set this flag. Just produce warning.
    }
    else { istringstream idata(txt); double x; while (idata >> x) m_VibFreq.push_back(x);}

    txt= ppPropList->XmlReadProperty("me:rotConsts");
    if(!txt){
      errorMsg << "ModelledMolecule::Cannot find argument me:rotConsts";
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
    }

    txt= ppPropList->XmlReadProperty("me:symmetryNumber");
    if(!txt){
      errorMsg << "ModelledMolecule::Cannot find argument me:symmetryNumber. Set to default value 1.0";
      m_Sym = 1.0;
      //setFlag(true);
    }
    else { istringstream idata(txt); idata >> m_Sym; }

    txt= ppPropList->XmlReadProperty("me:ZPE");
    if(!txt){
      errorMsg << "ModelledMolecule::Cannot find argument me:ZPE";
      setFlag(true); // there has to have zero point energy in all molecules.
    }
    else { istringstream idata(txt); idata >> m_ZPE ; }

    // Determine the method of DOS calculation.
    const char* pDOSCMethodtxt = pp->XmlReadValue("me:DOSCMethod", false) ;
    if(pDOSCMethodtxt)
    {
      m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
      if(!m_pDensityOfStatesCalculator)
      {
        stringstream errorMsg;
        errorMsg << "Unknown method " << pDOSCMethodtxt
          << " for the calculation of Density Of States in ModelledMolecule "
          << getName();
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        m_pDensityOfStatesCalculator = NULL;
      }
    }
    else{
      pDOSCMethodtxt = "Classical rotors"; // must exist
      m_pDensityOfStatesCalculator = DensityOfStatesCalculator::Find(pDOSCMethodtxt);
    }

    txt= ppPropList->XmlReadProperty("me:spinMultiplicity");
    if(!txt){
      errorMsg << "ModelledMolecule::Cannot find argument me:spinMultiplicity. Set to default value 1.0";
      m_SpinMultiplicity = 1;
    }
    else
    { 
      istringstream idata(txt);
      idata >> m_SpinMultiplicity;
    }

    if (getFlag()){
      errorMsg << "Error(s) while initializing ModelledMolecule: " << getName();
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
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
      errorMsg << "InitializeMolecule failed for ModelledMolecule " << getName() << " before constructing CollidingMolecule.";
      setFlag(true);
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      errorMsg << "Molecule::Cannot find argument me:sigma. Default value " << sigmaDefault << " used.\n";
      setSigma(sigmaDefault);
      //setFlag(true);
      // sigma and epsilon are not always necessary.
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt= ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      errorMsg << "Molecule::Cannot find argument me:epsilon. Default value " << epsilonDefault << " used.\n";
      setEpsilon(epsilonDefault);
      //setFlag(true);
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    txt= ppPropList->XmlReadProperty("me:deltaEDown");
    if(!txt){
      errorMsg << "CollidingMolecule::Cannot find argument me:deltaEDown";
      //setFlag(true);
    }
    else { istringstream idata(txt); idata >> m_DeltaEdown; }

    if (getFlag()){
      errorMsg << "Error(s) while initializing CollidingMolecule: " << getName();
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

//    {stringstream errorMsg;
//     errorMsg << "Constructed CollidingMolecule " << getName() << " successfully.";
//     obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);}

    return true;
  }

  bool SuperMolecule::InitializeMolecule(PersistPtr pp)
  {
    //the construction of SuperMolecule should always rely on the components, there is therefore no need to initialize it
    //apart from name and m_ppPersist
    stringstream errorMsg;

    if (pp) {
      setPersistentPointer(pp);

      const char* id = pp->XmlReadValue("id");
      if (id) setName(id);
      else{
        errorMsg << "Molecular name is absent.\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        string tempName = "source"; setName(tempName);
        setFlag(true);
      }
      return true;
    }
    else{
      errorMsg << "Invalid PersistPtr.\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    if (getFlag()){
      errorMsg << "Error(s) while initializing SuperMolecule: " << getName();
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return false;
    }
  }

  //
  // Get cell density of states.
  //
  void ModelledMolecule::getCellDensityOfStates(vector<double> &cellDOS, const MesmerEnv &mEnv) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;
    cellDOS.assign(m_cellDOS.begin(), m_cellDOS.end());
  }

  //
  // Get cell energies.
  //
  void ModelledMolecule::getCellEnergies(vector<double> &CellEne, const MesmerEnv &mEnv) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;
    CellEne.assign(m_cellEne.begin(), m_cellEne.end());
  }

  //
  // Get grain density of states.
  //
  void ModelledMolecule::grnDensityOfStates(vector<double> &grainDOS, const MesmerEnv &mEnv) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;
    grainDOS.assign(m_grainDOS.begin(), m_grainDOS.end());
  }

  //
  // Get grain energies.
  //
  void ModelledMolecule::grnEnergies(vector<double> &grainEne, const MesmerEnv &mEnv) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;
    grainEne.assign(m_grainEne.begin(), m_grainEne.end());
  }

  //
  // Get Grain Boltzmann distribution.
  //
  void ModelledMolecule::grnBoltzDist(vector<double> &grainBoltzDist, const MesmerEnv &mEnv)
  {
    // If density of states have not already been calcualted then do so.
    {stringstream errorMsg;
     errorMsg << "Before calcDensityOfStates(), Molecular name: " << getName();
     obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);}

    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;

    {stringstream errorMsg;
     errorMsg << "After calcDensityOfStates(), Molecular name: " << getName();
     obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);}

    int MaximumGrain = mEnv.MaxGrn ;
    double beta = mEnv.beta;

    // Calculate the Boltzmann distribution.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.

    int i ;
    double prtfn(0.0) ;
    for (i = 0; i < MaximumGrain; i++) {
        double tmp = log(m_grainDOS[i]) - beta*m_grainEne[i] + 10.0 ;
        tmp = exp(tmp) ;
        prtfn += tmp ;
        grainBoltzDist[i] = tmp ;
    }

    // Normalize the Boltzmann distribution.

    for (i = 0; i < MaximumGrain; i++) {
        grainBoltzDist[i] /= prtfn ;
    }
  }

  //
  // Get Grain canonical partition function.
  //
  double ModelledMolecule::grnCanPrtnFn(const MesmerEnv &mEnv) {
    // If density of states have not already been calcualted then do so.
    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;

    double CanPrtnFn(0.0) ;

    // Calculate the ro-vibrational partition function based on the grain
    // densities of states, and not the molecular properties, for consistency.

    int MaximumGrain = mEnv.MaxGrn ;
    double beta = mEnv.beta;

    for (int i = 0; i < MaximumGrain; i++) {
      CanPrtnFn += exp( log(m_grainDOS[i]) - beta*m_grainEne[i] ) ;
    }

    // Electronic partition function.
    CanPrtnFn *= double(getSpinMultiplicity()) ;

    // Translational partition function.
    return CanPrtnFn ;
  }

  int ModelledMolecule::get_rotConsts(std::vector<double> &mmtsInt)
  {
    mmtsInt.clear();
    mmtsInt.push_back(m_RotCstA);
    mmtsInt.push_back(m_RotCstB);
    mmtsInt.push_back(m_RotCstC);
    /* now the classification of rotors is simplified to only three following types. 3-D rotors may have other
    attributes different from one another but in ILT they are treated as the same type. The function return values
    are temporary shorthand representations. */
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
  bool CollidingMolecule::initCollisionOperator(double beta, Molecule *pBathGasMolecule, const MesmerEnv &mEnv)
  {
    // If density of states have not already been calcualted then do so.

    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;

    // Calculate the collision frequency.

    m_collisionFrequency = collisionFrequency(beta, mEnv.conc, pBathGasMolecule) ;

    // Calculate the collision operator.

    if (!collisionOperator(beta, mEnv)){
      std::stringstream errorMsg;
      errorMsg << "Failed building collision operator";
      mesmer::obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obError);
      return false;
    }
    return true;
  }

  //
  // Calculate collision operator
  //
  bool CollidingMolecule::collisionOperator(double beta, const MesmerEnv &mEnv)
  {
    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //
    int i, j;
    int MaximumGrain = mEnv.MaxGrn;
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
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
    }
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

    //Normalisation
    m_egme->normalize();

    //account for collisional loss by subrtacting unity from the leading diagonal.
    for ( i = 0 ; i < MaximumGrain ; ++i ) (*m_egme)[i][i] -= 1.0 ;


    cout << endl << "Column Sums" << endl << endl ;
    for ( i = 0 ; i < MaximumGrain ; ++i ) {
      double columnSum(0.0) ;
      for ( j = 0 ; j < MaximumGrain ; ++j ){
        columnSum += to_double((*m_egme)[j][i]) ;
      }
      cout << columnSum << endl ;
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

  /*//
  // Diagonalize the Collision Operator. See ReactionManager::diagCollisionOperator()
  //
  void CollidingMolecule::  diagCollisionOperator() {

    // Allocate space for eigenvalues.

    int msize = m_egme->size() ;
    vector<double> rr(msize, 0.0) ;

    m_egme->diagonalize(&rr[0]) ;

    cout << endl ;
    for (int i(0) ; i < msize ; ++i)
        cout << rr[i] << endl;

    cout << endl ;
    for (int i(0) ; i < msize ; ++i)  {
      formatFloat(cout, m_grainEne[i],              6, 15) ;
      formatFloat(cout, (*m_egme)[i][MaximumGrain-1], 6, 15) ;
      formatFloat(cout, (*m_egme)[i][MaximumGrain-2], 6, 15) ;
      cout << endl ;
    }

  }
  */

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
    double bthEpsilon = pBathGasMolecule->getEpsilon();

    double mu   = amu * getMass() * bthMass/(getMass() + bthMass) ;
    double eam  = sqrt(getEpsilon() * bthEpsilon) ;
    double sam  = (getSigma() + bthSigma) * 0.5;
    double tstr = 1. / (eam * beta);

    // Calculate collision integral.
    double collFrq = A * exp(-log(tstr) * B) + C * exp(-D * tstr) + E * exp(-F * tstr) ;

    // Calculate molecular collision frequency.
    collFrq *= M_PI * sam * sam * 1.e-20 * sqrt(8. * 1.381e-23 * temp/(M_PI * mu)) ;

    // Calculate overall collision frequency.
    collFrq *= conc * 1.e+06 ;

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
                                                const double RducdOmega,
                                                const MesmerEnv &mEnv) const
  {
    // Find size of system matrix.

    int smsize = static_cast<int>(CollOptr->size()) ;
    int MaximumGrain = mEnv.MaxGrn;

    // Check there is enough space in system matrix.

    if (locate + size > smsize) {
      stringstream errorMsg;
      errorMsg << "Error in the size of the system matrix.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1) ;
    }

    // Copy collision operator to the diagonal block indicated by "locate"
    // and multiply by the reduced collision frequencey.
    int start = MaximumGrain - size;
    //CHL: I guess the well is concerned starting from the bottom (MaximumGrain - size) to the top (MaximumGrain)
    for (int i(0) ; i < size ; ++i) {
      int ii(locate + i) ; int ipos = i + start;
      for (int j(0) ; j < size ; ++j) {
        int jj(locate + j) ;
        (*CollOptr)[ii][jj] = RducdOmega * (*m_egme)[ipos][j + start] ;
      }
    }
  }

  //-------------------------------------------------------------------------------------------
  //
  // Private methods.

  //
  // Calculate the rovibrational density of states.
  //
  bool ModelledMolecule::calcDensityOfStates(const MesmerEnv &mEnv)
  {
    if (!m_pDensityOfStatesCalculator->countCellDOS(this, mEnv)){
      return false;
    }
    calcGrainAverages(mEnv); testDensityOfStates(mEnv) ;
    return true;
  }

  //
  // Test the rovibrational density of states for ModelledMolecule.
  //
  void ModelledMolecule::testDensityOfStates(const MesmerEnv &mEnv)
  {
    cout << endl << "Test density of states for ModelledMolecule: " << getName() << endl << endl ;
    int MaximumGrain = mEnv.MaxGrn;
    int MaximumCell  = mEnv.MaxCell;

    string comment("Partition function calculation at various temperatures.\n qtot : partition function as a product of quantum mechanical partition functions for vibrations (1-D harmonic oscillator) and classifical partition functions for rotations\n sumc : (user calculated) cell based partition function \n sumg : (user calculated) grain based partition function ");

    PersistPtr ppList = getPersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment );

    if (0)
    {
      cout << "Cell density of states:\n";
      for ( int i = 0 ; i < MaximumCell ; ++i ) {
        cout << m_cellDOS[i] << endl;
      }
      cout << "Grain density of states:\n";
      for ( int i = 0 ; i < MaximumGrain ; ++i ) {
        cout << m_grainDOS[i] << endl;
      }
    }
    
    cout << "      T           qtot           sumc           sumg\n";

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

      vector<double> rotConst; int rotorType = get_rotConsts(rotConst);

      switch(rotorType){
        case 2://3-D symmetric/asymmetric/spherical top
          for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) {
            qtot /= (1.0 - exp(-beta*m_VibFreq[j])) ;
          }
          qtot *= (sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/m_Sym) ;
          break;
        case 0://2-D linear
          for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) {
            qtot /= (1.0 - exp(-beta*m_VibFreq[j])) ;
          }
          qtot *= (rotConst[0] / (m_Sym*beta)) ;
          break;
        default:
          qtot = 0.;
      }
      formatFloat(cout, temp,  6,  7) ;
      formatFloat(cout, qtot,  6, 15) ;
      formatFloat(cout, sumc,  6, 15) ;
      formatFloat(cout, sumg,  6, 15) ;
      cout << endl ;

      //Add to XML document
      PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
      ppItem->XmlWriteValueElement("me:T",    temp, 6);
      ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
      ppItem->XmlWriteValueElement("me:sumc", sumc, 6);
      ppItem->XmlWriteValueElement("me:sumg", sumg, 6);
    }
  }

  //
  // Calculate the average grain energy and then number of states per grain.
  //
  void ModelledMolecule::calcGrainAverages(const MesmerEnv &mEnv)
  {
    int MaximumGrain = mEnv.MaxGrn;
    m_grainEne.resize(MaximumGrain) ;
    m_grainDOS.resize(MaximumGrain, 0.) ;

    //    int igsz = MAXCELL/MAXGRN ;

    // Check that there are enough cells.

    if (mEnv.GrainSize < 1) {
      stringstream errorMsg;
      errorMsg << "The number of Cells is insufficient to produce requested number of Grains.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i = 0 ; i < MaximumGrain ; ++i ) {

      int idx3 = idx1 ;

      // Calculate the number of states in a grain.

      double gNOS = 0.0 ; 
      for (int j = 0 ; j < mEnv.GrainSize ; ++j, ++idx1 )
        gNOS += m_cellDOS[idx1] ;

      // Calculate average energy of the grain if it contains sum states.

      if ( gNOS > 0.0 ) {

        double gSE = 0.0 ; // grain sum of state energy
        for (int j = 0 ; j < mEnv.GrainSize ; ++j, ++idx3 )
          gSE += m_cellEne[idx3] * m_cellDOS[idx3] ;

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
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    else{
//      stringstream errorMsg;
//      errorMsg << "Number of grains requested: " << MaximumGrain << endl
//               << "Number of grains produced : " << idx2 << " in " << getName();
//      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

  }
}//namespace
