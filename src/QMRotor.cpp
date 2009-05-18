#include "QMRotor.h"

using namespace std;
namespace mesmer
{
  //************************************************************
  //Global instance, defining its id (usually the only instance) but here with an alternative name
  QMRotor theQMRotor("QMRotors");
  QMRotor oldQMRotor("QM rotors");
  //************************************************************

  // Provide a function to define particular counts of the DOS of a molecule.
  bool QMRotor::countCellDOS(gDensityOfStates* pDOS, int MaximumCell, PersistPtr ppDOSC)
  {
    vector<double> VibFreq ; 
    pDOS->get_VibFreq(VibFreq) ;

    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0) ;

    //
    // Initialize density of states array using calculated rotational
    // density of state.
    //

    //From inverse Laplace transform of rotors
    vector<double> rotConst; int rotorType = pDOS->get_rotConsts(rotConst);
    double sym = pDOS->get_Sym();
    double qele = pDOS->getSpinMultiplicity();
    int i_e(0);

    // Note: rotConst[0] (A) >= rotConst[1] (B) >= rotConst[2] (C)
    double rcA(rotConst[0]), rcB(rotConst[1]), rcC(rotConst[2]);

    switch (rotorType){
      // Although there is possibility that two or three rotational constants are very close but not equal, 
      // it is considered as they are different. When a user input rotational constants into the input file,
      // it is one's own decision to make slightly different rotational constants equal for alternative
      // interpretation; otherwise they will be treated as different even very close to each other.
    case 2: //3-D symmetric/asymmetric/spherical top
      if (rcA == rcC || ((rcA - rcC)/rcC < .01)){ // spherical top
        if ((rcA - rcC)/rcC < .01) rcA = (rcA + rcB + rcC) / 3.0;
        for (int j(0);; ++j ){
          i_e = int(rcA * (double)(j * (j + 1)));
          if (i_e > MaximumCell) break;
          int sqrdg(2 * j + 1);
          cellDOS[i_e] = qele * (double)(sqrdg * sqrdg) * sym;
        }
      }
      else{

        // Asymmetry parameter Kappa varies from -1 for a prolate symmetric top to 1 for an oblate symmetric top.
        double Kappa = (2. * rcB - rcA - rcC)/(rcA - rcC);

        if (Kappa > 0){ // near oblate symmetric top
          // A true oblate symmetric top has rotational constants A = B > C. The closer Kappa is to 1, the closer it 
          // is an oblate rotor.
          if (Kappa != 1.){
            if (Kappa < 0.9) cerr << "The rotor is not close to symmetric oblate, risking of incorrect interpretation.";
            rcB = .5 * (rcB + rcA);
          }
          int maxJ = int((-rcB + sqrt(rcB*rcB +4. * rcC * (double)(MaximumCell)))/
            (2. * rcC)); // A function to calculate the maximum J possible for MaximumCell.
          for (int j(0); j <= maxJ; ++j ){
            double d_ei = rcB * (double)(j * (j + 1)); // B J (J + 1)
            for (int k(-j) ; k <= j; ++k ){
              i_e = int (d_ei + (rcC - rcB) * (double)(k * k)); // B J (J + 1) + (C - B) K^2
              if (i_e >= MaximumCell){
                if (k < 0) k = -k;
                continue;
              }
              cellDOS[i_e] += qele * (double)(2 * j + 1) * sym;
              // use += operator because there can be duplicated energy JK sums
            }
          }
        }
        else if (Kappa < 0){ // near prolate symmetric top
          // A true prolate symmetric top has rotational constants A > B = C. The closer Kappa is to -1, the closer it 
          // is an prolate rotor.
          if (Kappa != -1.){
            if (Kappa > -0.9) cerr << "The rotor is not close to symmetric prolate, risking of incorrect interpretation.";
            rcB = .5 * (rcB + rcC);
          }
          for (int j(0);; ++j ){
            double d_ei = rcB * (double)(j * (j + 1)); // B J (J + 1)
            if (int(d_ei) >= MaximumCell) break; // if J exceeds the energy level, K contribution also ignored as A>B.
            for (int k(-j) ; k <= j; ++k ){
              i_e = int (d_ei + (rcA - rcB) * (double)(k * k)); // B J (J + 1) + (A - B) K^2
              if (i_e >= MaximumCell) break;
              cellDOS[i_e] += qele * (double)(2 * j + 1) * sym;
              // use += operator because there can be duplicated energy JK sums
            }
          }
        }
        else{ // complete asymmetric top (chances are rare)
          cinfo << "There is no appropriate way to use QM direct count for a complete asymmetric top rotor.";
          return false;
        }
      }
      break;
    case 0: //2-D linear
      for (int j(0);; ++j ){
        i_e = int(rcA * (double)(j * (j + 1)));
        if (i_e > MaximumCell){
          break;
        }
        cellDOS[i_e] += qele * (double)(2 * j + 1) * sym;
      }
      break;
    default:
      break;
      //DO nothing
      /*
      cnt = 0.;
      for (int i(0); i < MaximumCell ; ++i ) 
      cellDOS[i] = cnt ;
      */
    }

    // Implementation of the Beyer-Swinehart algorithm.
    Beyer_Swinehart(VibFreq, cellDOS);

    //electronic degeneracy
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
    if (!eleExc.empty()){
      for (int j(0); j < int(eleExc.size()) ; ++j){
        int iele = static_cast<int>(eleExc[j]);
        for (int i = (MaximumCell - 1); i >= (iele - 1); --i){
          cellDOS[i] += cellDOS[i - iele + 1];
        }
      }
    }

    pDOS->setCellDensityOfStates(cellDOS) ;

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
  // energy matrix.  ED1=E(1,1)+E(-1,1), ED2=E(1,1)-E(-1,1).
  //

  void QMRotor::asymmetricRotor(double A, double B, double C, int J, double *pkpp, double *Er, double *Ed, double *Ef) {

    const int NMAX = 200 ;
    vector<double> E(NMAX,0.0) ;
    vector<double> R(NMAX,0.0) ;
    
    double kpp = (2.0*B - A - C)/(A - C) ; // Eccentricity of top.

    double jj   = double(J) ;
    double jsqd = jj*(jj + 1.0) ;
    double f    =  (kpp - 1.0)/2.0 ;
    double h    = -(kpp + 1.0)/2.0 ;

    double Ed0 = f*jsqd ;

    for (int k = 0 ; k < J ; k++) {
      double kk = double(k + 1) ;
      Ed[k] = Ed0 + (1.0 - f)*kk*kk ;
      double Ee = (jsqd - (kk-2.0)*(kk-1.0))*(jsqd - (kk-1.0)*kk)/4.0 ;
      Ef[k] = h*sqrt(Ee) ;
    }

    Ef[2] *= sqrt(2.0) ;
    double Ed1 = Ed[1] + Ef[1] ;
    double Ed2 = Ed[1] - Ef[1] ;

    //
    // E+ Block.
    //
    int i(0) ;
    int N = J/2 + 1 ;
    R[0] = Ed0 ;
    for (i = 0 ; i < N-1 ; i++) {
      E[i+1] = Ef[i*2] ;
      R[i+1] = Ed[i*2] ;
    }

    TMatrix<double>::tqlev(&R[0], &E[0], N) ;
    
    int m0(0) ;
    double ene = (A+C)*jsqd/2.0 ;
    double en2 = (A-C)/2.0 ;
    for (i = 0 ; i < N ; i++, m0++) {
      Er[m0] = ene + en2*R[i] ;
    }

    //
    // E- Block.
    //
    N = J/2 ;
    for (i = 0 ; i < N ; i++) {
      E[i+1] = Ef[i*2+2] ;
      R[i]   = Ed[i*2] ;
    }
    
    TMatrix<double>::tqlev(&R[0], &E[0], N) ;
        
    for (i = 0 ; i < N ; i++, m0++) {
      Er[m0] = ene + en2*R[i] ;
    }

    //
    // O+ Block.
    //
    N = (J+1)/2 ;
    R[1] = Ed1 ;
    for (i = 0 ; i < N-1 ; i++) {
      E[i+1] = Ef[i*2+1] ;
      R[i+1] = Ed[i*2+1] ;
    }
    
    TMatrix<double>::tqlev(&R[0], &E[0], N) ;
    
    for (i = 0 ; i < N ; i++, m0++) {
      Er[m0] = ene + en2*R[i] ;
    }
    
    //
    // O- Block.
    //
    N = (J+1)/2 ;
    R[1] = Ed2 ;
    for (i = 0 ; i < N-1 ; i++) {
        E[i+1] = Ef[i*2+1] ;
        R[i+1] = Ed[i*2+1] ;
    }
    
    TMatrix<double>::tqlev(&R[0], &E[0], N) ;
    
    for (i = 0 ; i < N ; i++, m0++) {
      Er[m0] = ene + en2*R[i] ;
    }
    
    //
    // Sort eigenvalues.
    //
    //      DO 110 N = 1,M0 - 1
    //         DO 100 I = N + 1,M0
    //              IF (ER(N).GT.ER(I)) THEN
    //                  EE    = ER(N)
    //                  ER(N) = ER(I)
    //                  ER(I) = EE
    //              END IF
    //  100     CONTINUE
    //  110 CONTINUE
    //

  }

}//namespace
