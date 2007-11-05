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
    m_Mass(0.0),
    m_Sigma(0.0),
    m_Epsilon(0.0),
    m_SpinMultiplicity(1)
  {}

  ModelledMolecule::ModelledMolecule():
    m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(0.0),
    m_ZPE(0.0),
    m_cellEne(),
    m_cellDOS()
  {}

  CollidingMolecule::CollidingMolecule():
    m_egme(NULL),
    m_DeltaEdown(0.0),
    m_ncolloptrsize(0),
    m_collisionFrequency(0.0)
  {}

  ModelledMolecule::~ModelledMolecule()
  {
    // Free any memory assigned for calculating densities of states.
    if (m_cellDOS.size()) m_cellDOS.clear();
    if (m_cellEne.size()) m_cellEne.clear();
  }

  CollidingMolecule::~CollidingMolecule()
  {
    if (m_egme != NULL) delete m_egme ;
  }

  TransitionState::TransitionState()
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
  bool Molecule::InitializeMolecule(PersistPtr p)
  {
    m_ppPersist = p;
    const char* id= m_ppPersist->XmlReadValue("id");
    if(id)
        m_Name = id;
    return !m_Name.empty();
  }

  bool BathGasMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for Molecule before constructing BathGasMolecule.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:epsilon";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_Epsilon; } //extra block ensures idata is initiallised

    txt= ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:sigma";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_Sigma; }

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "BathGasMolecule::Cannot find argument me:MW";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_Mass; }

    return true;
    //TODO append name of molecule to error message
  }


  bool ModelledMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    if(!Molecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for Molecule before constructing ModelledMolecule.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:ZPE");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "ModelledMolecule::Cannot find argument me:ZPE";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_ZPE ; }

    txt= ppPropList->XmlReadProperty("me:rotConsts");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "ModelledMolecule::Cannot find argument me:rotConsts";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
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
      stringstream errorMsg;
      errorMsg << "ModelledMolecule::Cannot find argument me:symmetryNumber";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_Sym; }

    txt= ppPropList->XmlReadProperty("me:vibFreqs");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "ModelledMolecule::Cannot find argument me:vibFreqs";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); double x; while (idata >> x) m_VibFreq.push_back(x); }

//    txt= ppPropList->XmlReadProperty("me:spinMultiplicity");
//    if(!txt)
//        return false;
//    { istringstream idata(txt);
//    idata >> m_SpinMultiplicity;
//    }
    return true;
  }

  bool CollidingMolecule::InitializeMolecule(PersistPtr pp)
  {
    //Read base class parameters first
    if(!ModelledMolecule::InitializeMolecule(pp)){
      stringstream errorMsg;
      errorMsg << "InitializeMolecule failed for ModelledMolecule before constructing CollidingMolecule.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
        ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "CollidingMolecule::me:epsilon not assigned for this molecule";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    else { istringstream idata(txt); idata >> m_Epsilon; } //extra block ensures idata is initiallised

    txt= ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "CollidingMolecule::me:sigma not assigned for this molecule";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }
    else { istringstream idata(txt); idata >> m_Sigma; }

    txt= ppPropList->XmlReadProperty("me:MW");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "CollidingMolecule::Cannot find argument me:MW";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_Mass; }

    txt= ppPropList->XmlReadProperty("me:deltaEDown");
    if(!txt){
      stringstream errorMsg;
      errorMsg << "CollidingMolecule::Cannot find argument me:deltaEDown";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    else { istringstream idata(txt); idata >> m_DeltaEdown; }

    return true;
  }

  //
  // Get cell density of states.
  //
  void ModelledMolecule::cellDensityOfStates(vector<double> &cellDOS, const MesmerEnv &mEnv) {

    // If density of states have not already been calcualted then do so.

    if (!m_cellDOS.size())
        calcDensityOfStates(mEnv) ;

    for (int i = 0 ; i < mEnv.MaxCell ; ++i )
        cellDOS[i] = m_cellDOS[i] ;
  }

  //
  // Get cell energies.
  //
  void ModelledMolecule::cellEnergies(vector<double> &CellEne, const MesmerEnv &mEnv) {

    // If energies have not already been calcualted then do so.
    if (!m_cellDOS.size())
        calcDensityOfStates(mEnv) ;

    for (int i = 0 ; i < mEnv.MaxCell ; ++i )
        CellEne[i] = m_cellEne[i] ;
  }

  //
  // Get grain density of states.
  //
  void ModelledMolecule::grnDensityOfStates(vector<double> &grainDOS, const MesmerEnv &mEnv) {

    // If density of states have not already been calcualted then do so.

    if (!m_cellDOS.size())
        calcDensityOfStates(mEnv) ;

    grainDOS = m_grainDOS ;
  }

  //
  // Get grain energies.
  //
  void ModelledMolecule::grnEnergies(vector<double> &grainEne, const MesmerEnv &mEnv) {

    // If energies have not already been calculated then do so.

    if (!m_cellDOS.size())
        calcDensityOfStates(mEnv) ;

    grainEne = m_grainEne ;
  }

  //
  // Get Grain Boltzmann distribution.
  //
  void ModelledMolecule::grnBoltzDist(vector<double> &grainBoltzDist, const MesmerEnv &mEnv)
  {
      // If energies have not already been calculated then do so.

      if (!m_cellDOS.size())
          calcDensityOfStates(mEnv) ;

      int MaximumGrain = mEnv.MaxGrn ;
      double beta = 1.0/( boltzmann_RCpK * mEnv.temp ) ;

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

    if (!m_cellDOS.size())
      calcDensityOfStates(mEnv) ;

      double CanPrtnFn(0.0) ;

      // Calculate the ro-vibrational partition function based on the grain
      // densities of states, and not the molecular properties, for consistency.

      int MaximumGrain = mEnv.MaxGrn ;
      double beta = 1.0/( boltzmann_RCpK * mEnv.temp ) ;

      for (int i = 0; i < MaximumGrain; i++) {
        double tmp = log(m_grainDOS[i]) - beta*m_grainEne[i];
        tmp = exp(tmp) ;
        CanPrtnFn += tmp ;
      }

    // Electronic partition function.
    CanPrtnFn *= double(m_SpinMultiplicity) ;

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
  void CollidingMolecule::initCollisionOperator(double beta, Molecule *pBathGasMolecule, const MesmerEnv &mEnv)
  {
    // If density of states have not already been calcualted then do so.

    if (!m_cellDOS.size())
        calcDensityOfStates(mEnv) ;

    // Calculate the collision frequency.

    m_collisionFrequency = collisionFrequency(beta, mEnv.conc, pBathGasMolecule) ;

    // Calculate the collision operator.

    collisionOperator(beta, mEnv) ;
  }

  //
  // Calculate collision operator
  //
  void CollidingMolecule::collisionOperator(double beta, const MesmerEnv &mEnv)
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
        return ; //***THIS FUNCTION NEEDS AN ERROR RETURN****
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

    double mu   = amu * m_Mass * bthMass/(m_Mass + bthMass) ;
    double eam  = sqrt(m_Epsilon * bthEpsilon) ;
    double sam  = (m_Sigma + bthSigma) * 0.5;
    double tstr = temp/eam ;

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
                                                const double RducdOmega) const
  {
    // Find size of system matrix.

    int smsize = CollOptr->size() ;

    // Check there is enough space in system matrix.

    if (locate + size > smsize) {
      stringstream errorMsg;
      errorMsg << "Error in the size of the system matrix.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
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
  void ModelledMolecule::calcDensityOfStates(const MesmerEnv &mEnv)
  {
    cout << endl << "Number of frequencies: " << m_VibFreq.size() << endl ;

    m_cellDOS.clear(); m_cellEne.clear(); //make sure there is no residue left

    //
    // Initialize density of states array using calculated rotational
    // density of state.
    //

    //From inverse Laplace transform of rotors
    vector<double> rotConst; int rotorType = get_rotConsts(rotConst);
    double cnt = 0.;
    
    switch (rotorType){
      case 2: //3-D symmetric/asymmetric/spherical top
        cnt = sqrt(4./(rotConst[0] * rotConst[1] * rotConst[2]))/m_Sym ;
        for (int i = 0 ; i < mEnv.MaxCell ; ++i ) {
          m_cellEne.push_back(static_cast<double>(i) + 0.5);
          m_cellDOS.push_back(cnt*sqrt(m_cellEne[i]));
        }
        break;
      case 0: //2-D linear
        cnt = 1./ (rotConst[0] * m_Sym);
        for (int i = 0 ; i < mEnv.MaxCell ; ++i ){
          m_cellEne.push_back(static_cast<double>(i) + 0.5);
          m_cellDOS.push_back(cnt);
        }
        break;
      default:
        cnt = 0.;
        for (int i = 0 ; i < mEnv.MaxCell ; ++i ){
          m_cellEne.push_back(static_cast<double>(i) + 0.5);
          m_cellDOS.push_back(cnt);
        }
    }


    // Implementation of the Bayer-Swinehart algorithm.

    for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) {
      int iFreq = static_cast<int>(m_VibFreq[j] +.5) ;
      // +.5 makes sure it floors to the frequency that is closer to the integer
      // the original case has larger difference where if frequency = 392.95 it floors to 392
      for (int i = 0 ; i < mEnv.MaxCell - iFreq ; ++i ){
        m_cellDOS[i + iFreq] += m_cellDOS[i] ;
      }
    }

    // Calculate grain averages.

    calcGrainAverages(mEnv) ;

    //    if ( get_verbosity() )
    testDensityOfStates(mEnv) ;

    return ;
  }

  //
  // Test the rovibrational density of states.
  //
  void ModelledMolecule::testDensityOfStates(const MesmerEnv &mEnv)
  {
    cout << endl << "Test density of states: " << m_Name << endl << endl ;
    int MaximumGrain = mEnv.MaxGrn;
    double temp ;
    double beta ;

    string comment("Partition function calculation at various temperatures.\n qtot : using analytical formula\n sumc : based on cells\n sumg  : based on grains");

    PersistPtr ppList = m_ppPersist->XmlWriteMainElement("me:densityOfStatesList", comment );
    for ( int n = 0 ; n < 29 ; ++n ) {

      temp = 100.0*static_cast<double>(n + 2) ;
      beta = 1.0/(boltzmann_RCpK*temp) ;

      // Calculate partition functions based on cells.

      double sumc  = 0.0 ;
      for ( int i = 0 ; i < mEnv.MaxCell ; ++i ) {
        sumc += m_cellDOS[i]*exp(-beta*m_cellEne[i]) ;
      }

      // Calculate partition functions based on grains.

      double sumg  = 0.0 ;
      for ( int i = 0 ; i < MaximumGrain ; ++i ) {
        sumg += m_grainDOS[i]*exp(-beta*m_grainEne[i]) ;
      }

      // Calculate partition functions using analytical formula.

      double qtot = 1.0 ;
      
      vector<double> rotConst; int rotorType = get_rotConsts(rotConst);
      
      switch(rotorType){
        case 2:
          for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) {
            qtot /= (1.0 - exp(-beta*m_VibFreq[j])) ;
          }
          qtot *= sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/m_Sym ;
          break;
        case 0:
          for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) {
            qtot /= (1.0 - exp(-beta*m_VibFreq[j])) ;
          }
          qtot *= rotConst[0] / (m_Sym*beta) ;
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
    m_grainDOS.resize(MaximumGrain) ;

    //    int igsz = MAXCELL/MAXGRN ;

    // Check that there are enough cells.

    if (mEnv.GrainSize < 1) {
      stringstream errorMsg;
      errorMsg << "Not enought Cells to produce requested number of Grains.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i = 0 ; i < MaximumGrain ; ++i ) {

      int idx3 = idx1 ;

      // Calculate the number of states in a grain.

      double smt = 0.0 ;
      for (int j = 0 ; j < mEnv.GrainSize ; ++j, ++idx1 )
          smt += m_cellDOS[idx1] ;

      // Calculate average energy of the grain if it contains sum states.

      if ( smt > 0.0 ) {

        double smat = 0.0 ;
        for (int j = 0 ; j < mEnv.GrainSize ; ++j, ++idx3 )
          smat += m_cellEne[idx3] * m_cellDOS[idx3] ;

        m_grainDOS[idx2] = smt ;
        m_grainEne[idx2] = smat/smt ;
        idx2++ ;
      }
    }

    // Issue warning if number of grains produced is less that requested.

    if ( idx2 < MaximumGrain ) {
      stringstream errorMsg;
      errorMsg << "Number of grains produced is less than requested" << endl
               << "Number of grains requested: " << MaximumGrain << endl
               << "Number of grains produced : " << idx2 << ".";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
  }
}//namespace
