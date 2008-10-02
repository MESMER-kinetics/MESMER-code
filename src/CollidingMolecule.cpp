//
// CollidingMolecule.cpp
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
  CollidingMolecule::CollidingMolecule(const MesmerEnv& Env, MesmerFlags& Flags) : ModelledMolecule(Env, Flags),
    m_Sigma(sigmaDefault),
    m_Epsilon(epsilonDefault),
    m_DeltaEdownExponent(0.0),
    m_DeltaEdownRefTemp(298.0),
    m_DeltaEdown(0.0),
    m_collisionFrequency(0.0),
    m_ncolloptrsize(0),
    m_pDistributionCalculator(NULL),
    m_Sigma_chk(-1),
    m_Epsilon_chk(-1),
    m_DeltaEdown_chk(-1),
    m_grainFracBeta(0.),
    m_grainDist(0),
    m_egme(NULL)
  {}

  CollidingMolecule::~CollidingMolecule()
  {
    if (m_Sigma_chk == 0){
      cinfo << "m_Sigma is provided but not used in " << getName() << endl;
    }
    if (m_Epsilon_chk == 0){
      cinfo << "m_Epsilon is provided but not used in " << getName() << endl;
    }
    if (m_DeltaEdown_chk == 0){
      cinfo << "m_DeltaEdown is provided but not used in " << getName() << endl;
    }
    if (m_egme != NULL) delete m_egme ;
    if (m_grainDist.size()) m_grainDist.clear();
  }

  //
  //Initialization
  //
  bool CollidingMolecule::InitializeMolecule(PersistPtr pp)
  {
    PersistPtr oldpp = pp;

    //Read base class parameters first
    if(!ModelledMolecule::InitializeMolecule(pp)){
      cinfo << "InitializeMolecule failed for " << getName()
        << " before constructing CollidingMolecule." << endl;
      setFlag(true);
    }

    pp = oldpp;

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if(!ppPropList)
      ppPropList=pp; //Be forgiving; we can get by without a propertyList element

    const char* txt;

    txt= ppPropList->XmlReadProperty("me:sigma");
    if(!txt){
      cinfo << "Cannot find argument me:sigma in " << getName()
        << ". Default value " << sigmaDefault << " is used.\n" << endl;
      // sigma and epsilon are not always necessary.
    }
    else { istringstream idata(txt); double sigma(0.); idata >> sigma; setSigma(sigma);}

    txt= ppPropList->XmlReadProperty("me:epsilon");
    if(!txt){
      cinfo << "Cannot find argument me:epsilon in " << getName()
        << ". Default value " << epsilonDefault << " is used.\n" << endl;
      // sigma and epsilon are not always necessary.
    }
    else { istringstream idata(txt); double epsilon(0.); idata >> epsilon; setEpsilon(epsilon);} //extra block ensures idata is initiallised

    txt= ppPropList->XmlReadProperty("me:deltaEDown");
    if(!txt){
      cinfo << "Cannot find argument me:deltaEDown in " << getName() << endl;
      //setFlag(true);
      // deltaEDown is not always necessary. Hoever, it is not wise to provide a default value.
    }
    else {
      istringstream idata(txt);
      double value(0.0);
      idata >> value;
      const char* pLowertxt    = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "lower");
      const char* pUppertxt    = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "upper");
      const char* pStepStxt    = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "stepsize");
      const char* pRefTemptxt  = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "referenceTemperature");
      const char* pExponenttxt = ppPropList->XmlReadPropertyAttribute("me:deltaEDown", "exponent");
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
    const char* pDistCalcMethodtxt = pp->XmlReadValue("me:DistributionCalcMethod", false) ;
    if(pDistCalcMethodtxt)
    {
      m_pDistributionCalculator = DistributionCalculator::Find(pDistCalcMethodtxt);
      if(!m_pDistributionCalculator) // if the provided method cannot be found,
      {
        cinfo << "Unknown method " << pDistCalcMethodtxt
          << " for the calculation of distribution fraction in " << getName()
          << ". Please check spelling error. Default method <Boltzmann> is used." << endl;
        pDistCalcMethodtxt = "Boltzmann";
        m_pDistributionCalculator = DistributionCalculator::Find(pDistCalcMethodtxt);
      }
    }
    else{ // if no method is provided.
      cinfo << "No method for the calculation of distribution fraction in " << getName()
        << " is provided. Default method <Boltzmann> is used." << endl;
      pDistCalcMethodtxt = "Boltzmann"; // must exist
      m_pDistributionCalculator = DistributionCalculator::Find(pDistCalcMethodtxt);
    }

    if (getErrorFlag()){
      cerr << "Error(s) while initializing: " << getName();
      return false;
    }
    return true;
  }

  double CollidingMolecule::get_collisionFrequency() const {
    return m_collisionFrequency ;
  } ;

  void CollidingMolecule::set_colloptrsize(int ncolloptrsize) {
    m_ncolloptrsize = ncolloptrsize ;
  } ;

  int  CollidingMolecule::get_colloptrsize() const {
    return m_ncolloptrsize ;
  } ;

  void   CollidingMolecule::setSigma(double value)          {
    m_Sigma = value;
    m_Sigma_chk = 0;
  } ;

  double CollidingMolecule::getSigma()                      {
    if (m_Sigma_chk >= 0){
      ++m_Sigma_chk;
      return m_Sigma ;
    }
    else{
      cerr << "m_Sigma was not defined but requested in " << getName()
        << ". Default value " << sigmaDefault << " is used.\n";
      //exit(1);
      return m_Sigma ;
    }
  } ;

  const int CollidingMolecule::get_grnZPE(){
    double grnZpe = (get_zpe() - getEnv().EMin) / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cerr << "Grain zero point energy has to be a non-negative value.";

    return int(grnZpe);
  }

  void   CollidingMolecule::setEpsilon(double value)        {
    m_Epsilon = value;
    m_Epsilon_chk = 0;
  } ;

  double CollidingMolecule::getEpsilon()                    {
    if (m_Epsilon_chk >= 0){
      ++m_Epsilon_chk;
      return m_Epsilon ;
    }
    else{
      cerr << "m_Epsilon was not defined but requested in " << getName()
        << ". Default value " << epsilonDefault << " is used.\n";
      //exit(1);
      return m_Epsilon ;
    }
  } ;

  double CollidingMolecule::getDeltaEdown()                    {
    if (m_DeltaEdown_chk >= 0){
      ++m_DeltaEdown_chk;
      const double refTemp = getDeltaEdownRefTemp();
      const double dEdExp = getDeltaEdownExponent();
      const double dEdRef = m_DeltaEdown.get_value();
      const double temperature = 1. / (boltzmann_RCpK * getEnv().beta);
      return dEdRef * pow((temperature/refTemp),dEdExp);
    }
    else{
      cerr << "m_DeltaEdown was not defined but requested in " << getName();
      exit(1);
    }
  } ;

  //
  // Initialize the Collision Operator.
  //
  bool CollidingMolecule::initCollisionOperator(double beta, BathGasMolecule *pBathGasMolecule)
  {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates()){
      cerr << "Failed calculating DOS";
      return false;
    }

    // Calculate the collision frequency.
    m_collisionFrequency = collisionFrequency(beta, getEnv().conc, pBathGasMolecule) ;

    // Calculate the collision operator.
    if (!collisionOperator(beta)){
      cerr << "Failed building collision operator";
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

    double DEDown = getDeltaEdown();
    if(!DEDown){
      cerr << "me:deltaEDown is necessary for " << getName() << ". Correct input file to remove this error.";
      return false;
    }

    double alpha = 1.0/DEDown ;

    // Allocate memory.
    if (m_egme) delete m_egme ;                       // Delete any existing matrix.
    m_egme = new dMatrix(m_ncolloptrsize) ;           // Collision operator matrix.

    // Initialisation and error checking.
    for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
      if (m_grainDOS[i] <= 0.0) {
        cerr << "Data indicates that grain " << i << " of the current colliding molecule has no states.";
        return false;
      }
    }

    // The collision operator.
    for ( i = 0 ; i < m_ncolloptrsize ; ++i ) {
      double ei = m_grainEne[i] ;
      double ni = m_grainDOS[i] ;
      for ( j = i ; j < m_ncolloptrsize ; ++j ) {

        // Transfer to lower Energy -
        double transferDown = exp(-alpha*(m_grainEne[j] - ei)) ;
        (*m_egme)[i][j] = transferDown;

        // Transfer to higher Energy (via detailed balance) -
        double transferUp = (*m_egme)[i][j] * (m_grainDOS[j]/ni)*exp(-beta*(m_grainEne[j] - ei)) ;
        (*m_egme)[j][i] = transferUp;
      }
    }

    //Normalisation
    normalizeCollisionOperator();

    // print out of column sums to check normalization results
    if (getFlags().reactionOCSEnabled){
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
      popDist.push_back(sqrt(exp(log(m_grainDOS[idx]) - beta * m_grainEne[idx] + 10.0)));
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

    return true;
  }

  //
  // Normalize collision operator
  //
  void CollidingMolecule::normalizeCollisionOperator(){

    vector<double> work(m_ncolloptrsize) ;// Work space.
    //
    // Normalization of Probability matrix.
    // Normalising coefficients are found by using the fact that column sums
    // are unity. The procedure leads to a matrix that is of upper triangular
    // form and the normalisation constants are found by back substitution.
    //

    int i, j; //int makes sure the comparison to negative numbers meaningful (i >=0)

    double scaledRemain(0.0) ;
    for ( i = m_ncolloptrsize - 1 ; i >= 0 ; --i ) {

      double upperSum(0.0) ;
      for ( j = 0 ; j <= i ; ++j )
        upperSum += (*m_egme)[j][i] ;

      if (upperSum > 0.0){
        if (i < (int)m_ncolloptrsize - 1){
          scaledRemain = 0.0;
          for ( j = i + 1 ; j < (int)m_ncolloptrsize ; ++j )
            scaledRemain += (*m_egme)[j][i] * work[j] ;
        }
        work[i] = (1.0 - scaledRemain) / upperSum ;
      }
    }

    //
    // Apply normalization coefficients
    //
    for ( i = 0 ; i < (int)m_ncolloptrsize ; ++i ) {
      (*m_egme)[i][i] *= work[i] ;
      for ( j = i + 1 ; j < (int)m_ncolloptrsize ; ++j ) {
        (*m_egme)[j][i] *= work[j] ;
        (*m_egme)[i][j] *= work[j] ;
      }
    }

    //(*m_egme).showFinalBits(m_ncolloptrsize);
  }

  //
  // Calculate collision frequency.
  //
  double CollidingMolecule::collisionFrequency(double beta, const double conc, BathGasMolecule *pBathGasMolecule)
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
    bthMass = pBathGasMolecule->getMass();

    double bthSigma = 0.0;
    bthSigma = pBathGasMolecule->getSigma();

    if (!bthSigma)
      cerr << "me:sigma is necessary for " << pBathGasMolecule->getName()
      << ". Correct input file to remove this error." << endl;

    double bthEpsilon = 0.0;
    bthEpsilon = pBathGasMolecule->getEpsilon();

    if (!bthEpsilon)
      cerr << "me:epsilon is necessary for " << pBathGasMolecule->getName()
      << ". Correct input file to remove this error.";
    double mu   = amu * getMass() * bthMass/(getMass() + bthMass) ;
    double eam  = sqrt(getEpsilon() * bthEpsilon) ;
    double sam  = (getSigma() + bthSigma) * 0.5;
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
  void CollidingMolecule::copyCollisionOperator(qdMatrix *CollOptr,
    const int size,
    const int locate,
    const double RducdOmega) const
  {
    // Find size of system matrix.

    int smsize = static_cast<int>(CollOptr->size()) ;
    //int MaximumGrain = getEnv().MaxGrn;

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
  // calculates p(E)*exp(-EB)
  //
  void CollidingMolecule::grainDistribution(vector<double> &grainFrac, const int totalGrnNumber)
  {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";

    if (m_grainDist.size() != m_grainDOS.size() || getEnv().beta != m_grainFracBeta){
      m_pDistributionCalculator->calculateDistribution(m_grainDOS, m_grainEne, getEnv().beta, m_grainDist);
      m_grainFracBeta = getEnv().beta;
    }

    for (int i = 0; i < totalGrnNumber; ++i){
      grainFrac.push_back(m_grainDist[i]);
    }
  }

  //
  // Get normalized grain distribution.
  //
  void CollidingMolecule::normalizedInitialDistribution(vector<double> &grainFrac, const int totalGrnNumber)
  {
    grainDistribution(grainFrac, totalGrnNumber);

    double prtfn(0.);
    for (int i = 0; i < totalGrnNumber; ++i){
      prtfn += grainFrac[i];
    }

    for (int i = 0; i < totalGrnNumber; ++i){
      grainFrac[i] /= prtfn;
    }

    if (getFlags().grainBoltzmannEnabled){
      ctest << "\nGrain fraction for " << getName() << ":\n{\n";
      for (int i = 0; i < totalGrnNumber; ++i){
        ctest << grainFrac[i] << endl;
      }
      ctest << "}\n";
    }
  }

  //
  // Get normalized cell distribution.
  //
  void CollidingMolecule::normalizedCellBoltzmannDistribution(vector<double> &cellFrac, const int startingCell)
  {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";
    
    cellFrac.clear();

    const int MaximumCell = getEnv().MaxCell;
    vector<double> tempCellFrac, cellEne;
    getCellEnergies(MaximumCell, cellEne);
    
    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    for (int i = 0; i < MaximumCell; ++i) {
      tempCellFrac.push_back(exp(log(m_cellDOS[i]) - getEnv().beta * cellEne[i] + 10.0));
    }

    double prtfn(0.);
    for (int i = 0; i < MaximumCell; ++i)
      prtfn += tempCellFrac[i];

    for (int i = 0; i < MaximumCell; ++i)
      tempCellFrac[i] /= prtfn;

    for (int i = startingCell; i < MaximumCell; ++i)
      cellFrac.push_back(tempCellFrac[i]);

  }


  //
  // Get normalized grain distribution.
  //
  void CollidingMolecule::normalizedGrnBoltzmannDistribution(vector<double> &grainFrac, const int totalGrnNumber, const int startGrnIdx, const int ignoreCellNumber)
  {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";

    vector<double> tempGrnFrac;
    grainFrac.clear();

    // Calculate unnormalized Boltzmann dist.
    // Note the extra 10.0 is to prevent underflow, it is removed during normalization.
    for (int i = 0; i < totalGrnNumber; ++i) {
      tempGrnFrac.push_back(exp(log(m_grainDOS[i]) - getEnv().beta * m_grainEne[i] + 10.0));
    }

    double prtfn(0.);
    for (int i = 0; i < totalGrnNumber; ++i){
      prtfn += tempGrnFrac[i];
    }

    for (int i = 0; i < totalGrnNumber; ++i){
      tempGrnFrac[i] /= prtfn;
    }

    //---------------------------
    //
    if (ignoreCellNumber == 0){
      grainFrac = tempGrnFrac;
    }
    else{
      // As there are cells ignored, the population of the first grain participates the reaction will be changed.
      // deal with the partial grain.
      const int MaximumCell = getEnv().MaxCell;
      const int gsz = getEnv().GrainSize;
      const int cellOffset = get_cellOffset();
      const int grnStartCell = startGrnIdx * gsz - cellOffset;
      double partialDOS(0.0);
      for (int i(ignoreCellNumber); i < gsz; ++i){
        partialDOS += m_cellDOS[i + grnStartCell];
      }
      vector<double> cellEne;
      getCellEnergies(MaximumCell, cellEne);
      double partialAvgEne(0.0);
      for (int i(ignoreCellNumber); i < gsz; ++i){
        partialAvgEne += m_cellDOS[i + grnStartCell] * cellEne[i + grnStartCell];
      }
      partialAvgEne /= partialDOS;
      double partialFrac = exp(log(partialDOS) - getEnv().beta * partialAvgEne + 10.0);
      partialFrac /= prtfn;
      grainFrac.push_back(partialFrac);
      for (int i(startGrnIdx+1); i < int(tempGrnFrac.size()); ++i){
        grainFrac.push_back(tempGrnFrac[i]);
      }
    }
    //---------------------------

  }

  void CollidingMolecule::set_DistributionCalculator(DistributionCalculator* value){
    m_pDistributionCalculator = value;
  }

  DistributionCalculator* CollidingMolecule::get_DistributionCalculator(){
    return m_pDistributionCalculator;
  }
}//namespace
