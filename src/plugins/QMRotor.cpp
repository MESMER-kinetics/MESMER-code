#include "../DensityOfStates.h"
#include "../Molecule.h"
#include "../gDensityOfStates.h"

using namespace std;
namespace mesmer
{
  class QMRotor : public DensityOfStatesCalculator
  {
  public:

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, const MesmerEnv& env);

    // Function to calculate contribution to canonical partition function.
    virtual void canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    QMRotor(const char* id) : m_id(id) { Register(); }

    virtual const char* getID()  { return m_id; }
    virtual bool includesRotations(){return true;}

    virtual QMRotor* Clone() { return new QMRotor(*this); }

    virtual ~QMRotor() {}

  private:
    const char* m_id;

  private:

    void asymmetricRotor(double A, double B, double C, int J, double kpp, vector<double> &Er) ;

  } ;

  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  QMRotor theQMRotor("QMRotors");

  //************************************************************

  // Provide a function to define particular counts of the DOS of a molecule.
  bool QMRotor::countCellDOS(gDensityOfStates* pDOS, const MesmerEnv& env)
  {
    const size_t MaximumCell = env.MaxCell ;
    const double cellSize = env.CellSize ;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellSize, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0) ;

    //
    // Initialize density of states array using calculated quantum mechanical
    // rotational density of state.
    //
    vector<double> rotConst; 
    RotationalTop rotorType = pDOS->get_rotConsts(rotConst);
    double sym = pDOS->get_Sym();
    double qele = pDOS->getSpinMultiplicity();
    size_t i_e(0);

    // Note: rotConst[0] (A) >= rotConst[1] (B) >= rotConst[2] (C)
    double rcA(rotConst[0]), rcB(rotConst[1]), rcC(rotConst[2]);

    switch (rotorType) {

    case NONLINEAR: //3-D symmetric/asymmetric/spherical top

      // The following code tests for the type of top and, where possible, uses an analytic 
      // solution for the energy levels.

      if (rcA == rcC || ((rcA - rcC)/rcC < .01)) { // spherical top

        rcA = (rcA + rcB + rcC) / 3.0;
        for (int j(0);; ++j ){
          i_e = nint(rcA * double(j * (j + 1))/cellSize);
          if (i_e > MaximumCell) break;
          int sqrdg(2 * j + 1);
          cellDOS[i_e] = qele * double(sqrdg * sqrdg) / sym;
        }

      } else {

        // Asymmetry parameter Kappa varies from -1 for a prolate symmetric top to 1 for an oblate symmetric top.
        double Kappa = (2. * rcB - rcA - rcC)/(rcA - rcC);

        // if (0) { // Near symmetric top.
        if (abs(Kappa) > 0.95) { // Near symmetric top.

          double rcDiff(0.0) ;
          int maxJ(0) ;

          if (Kappa > 0.95) { // Near oblate symmetric top.

            // A true oblate symmetric top has rotational constants A = B > C.
            // Energy given by: E = B J (J + 1) + (C - B) K^2
            // The closer Kappa is to 1, the closer it is an oblate rotor.
            rcB    = (rcB + rcA) / 2.0 ;
            rcDiff = rcC - rcB ;
            // Determine the maximum J possible for MaximumCell.
            maxJ = int((-rcB + sqrt(rcB*rcB + 4.0*rcC * double(MaximumCell)))/(2.0*rcC)); 

          } else { // Near prolate symmetric top.

            // A true prolate symmetric top has rotational constants A > B = C.
            // Energy given by: E = B J (J + 1) + (A - B) K^2
            // The closer Kappa is to -1, the closer it is an prolate rotor.
            rcB    = (rcB + rcC) / 2.0 ;
            rcDiff = rcA - rcB ;
            // Determine the maximum J possible for MaximumCell.
            maxJ = int((-rcB + sqrt(rcB*rcB + 4.0*rcB * double(MaximumCell)))/(2.0*rcB)); 

          }

          for (int j(0); j <= maxJ; ++j ){
            double d_ei = rcB * double(j * (j + 1)); // B J (J + 1)
            for (int k(-j) ; k <= j; ++k ){
              i_e = nint((d_ei + rcDiff * double(k * k))/cellSize) ; 
              if (i_e < MaximumCell)
                cellDOS[i_e] += qele * double(2 * j + 1) / sym;
            }
          }          

        } else { // General asymmetric top.

          // This method can be expensive. Issue a warning if there is more than four non-Hydrogen atoms.
          if (pDOS->IsHeavyTop(4)) {
            string name(pDOS->getHost()->getName()) ;
            cinfo << "Warning: " << name << " is an asymmetric top containing five or more non-hydrogen" << endl ;
            cinfo << "atoms. The asymmetric top rotation states may be expensive to calculate. You may " << endl ;
            cinfo << "wish to consider using a classical treatment." << endl ;
          }

          bool withInRange(true) ;
          for (int j(0); withInRange ; ++j ) {
            vector<double> Er ;
            asymmetricRotor(rcA, rcB, rcC, j, Kappa, Er) ;
            withInRange = false ;
            for (size_t k(0); k < Er.size() ; ++k ){
              i_e = nint(Er[k]/cellSize) ;
              if (i_e < MaximumCell) {
                withInRange = true ;
                cellDOS[i_e] += qele * double(2 * j + 1) / sym;
              }
            }
          }
        }
      }
      break;
    case LINEAR: //2-D linear
      for (int j(0);; ++j ){
        i_e = nint(rcA * double(j * (j + 1))/cellSize);
        if (i_e > MaximumCell){
          break;
        }
        cellDOS[i_e] += qele * double(2 * j + 1) / sym;
      }
      break;
    default: // Assume atom.
      cellDOS[0] = qele ;
      break;
    }

    // Electronic excited states.
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
    vector<double> tmpCellDOS(cellDOS);
    for (size_t j(0) ; j < eleExc.size() ; ++j){
      size_t nr = nint(eleExc[j]/cellSize) ;
      if (nr < MaximumCell) {
        for (size_t i(0) ; i < MaximumCell - nr ; i++ ) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i] ;
        }
      }
    }

    pDOS->setCellDensityOfStates(tmpCellDOS) ;

    return true;
  }

  //
  // Energy levels of an asymmetric molecule. 
  //
  // The rotational constant B varies from A to C. The King, Hainer & Cross 
  // notation of energy levels is used (JCP, Vol. 11, p. 27 (1943)). The 
  // eigenfunctions are expanded in the prolate symmetric top basis set.
  // (See also Zare.)
  // 
  // ED(K)=E(K,K) are the diagonal elements of the energy matrix, 
  // ED0=E(0,0). EF(K)=E(K-2,K) are the off-diagonal elements of the 
  // energy matrix.
  //

  void QMRotor::asymmetricRotor(double A, double B, double C, int J, double kpp, vector<double> &Er) {

    int NMAX = max(3,2*J+1) ;
    vector<double> E(NMAX,0.0) ;
    vector<double> R(NMAX,0.0) ;
    vector<double> Ed(NMAX,0.0) ;
    vector<double> Ef(NMAX,0.0) ;

    double jsqd = double(J*(J + 1)) ;
    double f    =  (kpp - 1.0)/2.0 ;
    double h    = -(kpp + 1.0)/2.0 ;
    double Ed0  = f*jsqd ;

    for (int k = 0 ; k < J ; k++) {
      double kk = double(k + 1) ;
      Ed[k] = Ed0 + (1.0 - f)*kk*kk ;
      double Ee = (jsqd - (kk-2.0)*(kk-1.0))*(jsqd - (kk-1.0)*kk)/4.0 ;
      Ef[k] = h*sqrt(Ee) ;
    }
    Ef[1] *= sqrt(2.0) ;

    //
    // E+ Block.
    //
    int i(0) ;
    int N = J/2 + 1 ;
    R[0] = Ed0 ;
    E[0] = 0.0 ;
    for (i = 1 ; i < N ; i++) {
      E[i] = Ef[i*2-1] ;
      R[i] = Ed[i*2-1] ;
    }

    TMatrix<double>::tqlev(&R[0], &E[0], N) ;

    double ene = (A+C)*jsqd/2.0 ;
    double en2 = (A-C)/2.0 ;
    for (i = 0 ; i < N ; i++) {
      Er.push_back(ene + en2*R[i]) ;
    }

    //
    // E- Block.
    //
    N = J/2 ;
    for (i = 0 ; i < N ; i++) {
      E[i+1] = Ef[i*2+3] ;
      R[i]   = Ed[i*2+1] ;
    }

    TMatrix<double>::tqlev(&R[0], &E[0], N) ;

    for (i = 0 ; i < N ; i++) {
      Er.push_back(ene + en2*R[i]) ;
    }

    //
    // O+ Block.
    //
    N = (J+1)/2 ;
    R[0] = Ed[0] + Ef[0] ; // E(1,1) + E(-1,1)
    E[0] = 0.0 ;
    for (i = 1 ; i < N ; i++) {
      E[i] = Ef[i*2] ;
      R[i] = Ed[i*2] ;
    }

    TMatrix<double>::tqlev(&R[0], &E[0], N) ;

    for (i = 0 ; i < N ; i++) {
      Er.push_back(ene + en2*R[i]) ;
    }

    //
    // O- Block.
    //
    N = (J+1)/2 ;
    R[0] = Ed[0] - Ef[0] ; // E(1,1) - E(-1,1)
    E[0] = 0.0 ;
    for (i = 1 ; i < N ; i++) {
      E[i] = Ef[i*2] ;
      R[i] = Ed[i*2] ;
    }

    TMatrix<double>::tqlev(&R[0], &E[0], N) ;

    for (i = 0 ; i < N ; i++) {
      Er.push_back(ene + en2*R[i]) ;
    }

    return ;
  }

  // Calculate contribution to canonical partition function.
  void QMRotor::canPrtnFnCntrb(gDensityOfStates* gdos, double beta, double &PrtnFn, double &IntrlEne, double &varEne) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);
    double sym = gdos->get_Sym();

    double qtot(1.0), ene(0.0), var(0.0) ; 
    qtot *= double(gdos->getSpinMultiplicity());

    // Classical approximations are used here for convenience. These will
    // break down at low temperature, in this case the cell based solution
    // is probably more accurate, though the use of symmetry numbers will
    // also be incorrect under such condiditons.
    switch(rotorType){
    case NONLINEAR://3-D symmetric/asymmetric/spherical top
      qtot *= (sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/sym) ;
      ene   = 3.0/(2.0*beta) ;
      var   = 3.0*boltzmann_RCpK/2.0 ;
      break;
    case LINEAR://2-D linear
      qtot /= (rotConst[0]*sym*beta) ;
      ene   = 1.0/beta ;
      var   = boltzmann_RCpK ;
      break;
    default:
      break; // Assume atom.
    }

    // Electronic excited states.
    vector<double> eleExc;
    gdos->getEleExcitation(eleExc);
    double qelec(1.0), Eelec(0.0), varEelec(0.0)  ;
    for (size_t j(0) ; j < eleExc.size() ; ++j){
      qelec += exp(-beta*eleExc[j]) ;
      Eelec += eleExc[j]*exp(-beta*eleExc[j]) ;
      varEelec += eleExc[j]*eleExc[j]*exp(-beta*eleExc[j]) ;
    }
    Eelec /= qelec ;
    varEelec = varEelec/qelec - Eelec*Eelec ;

    PrtnFn   *= qtot*qelec ;
    IntrlEne += ene + Eelec ;
    varEne   += var + varEelec ;
  }  

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int QMRotor::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    unsigned int nDOF(0) ;
    switch(rotorType){
    case NONLINEAR:
      nDOF = 3 ;
      break;
    case LINEAR:
      nDOF = 2 ;
      break;
    default:
      // Assume atom.
      break; 
    }

    return nDOF ;
  }

}//namespace
