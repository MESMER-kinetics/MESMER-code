#include "MesmerMath.h"


//convolutes DOSs
void DOSconvolution(const std::vector<double> &f1,
                    const std::vector<double> &f2,
                    std::vector<double> &conv)
{
  for (unsigned int i = 0; i < f1.size(); ++i){
    conv[i] = 0.;
    for (unsigned int j = 0; j <= i; ++j)
      conv[i] += f1[j] * f2[i + 1 - j];
  }
}


//convolutes rovibrational DOSs
//bool convRV(const std::vector<double> &p_f1,
//            const std::vector<double> &p_f2,
//            const std::vector<double> &p_r1,
//            const std::vector<double> &p_r2,
//            const int &p_sym1,
//            const int &p_sym2,
//            const int &p_edg1,
//            const int &p_edg2)
//{
//  //populating vibrational frequencies
//
//  //declaring internal variables
//  std::vector<double> f1;
//  std::vector<double> f2;
//  std::vector<double> f3;
//  std::vector<double> r1;
//  std::vector<double> r2;
//  int sym1 = 1;
//  int sym2 = 1;
//  int edg1 = 1;
//  int edg2 = 1;
//
//  if (p_f1.size() < p_f2.size()){ // p_f1 has fewer vibrational modes -> revert their order
//    f1 = p_f2; f2 = p_f1;
//    r1 = p_r2; r2 = p_r1;
//    sym1 = p_sym2; sym2 = p_sym1;
//    edg1 = p_edg2; edg2 = p_edg1;
//  }
//  else{
//    f1 = p_f1; f2 = p_f2;
//    r1 = p_r1; r2 = p_r2;
//    sym1 = p_sym1; sym2 = p_sym2;
//    edg1 = p_edg1; edg2 = p_edg2;
//  }
//
//  //check vibrational density of states
//  if (!f1.size()){
//    std::stringstream errorMsg;
//    errorMsg << "No vibrational frequencies for the major reactant";
//    mesmer::obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obError);
//    return false;
//  }
//  if (static_cast<double>(sym2) > 1.0 && !f2.size()){
//    std::stringstream errorMsg;
//    errorMsg << "Symmetry number of the secondary reactant > 1 but no vibrational frequencies exists.";
//    mesmer::obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obError);
//    return false;
//  }
//
//  double cnt = 0.0;
//
//  for (int i = 0; i < f1.size(); ++i) f3.push_back(f1[i]);
//  if (f2.size()){
//    for (int i = 0; i < f2.size(); ++i) f3.push_back(f2[i]);
//    if (r1.size() == 1 && r2.size() == 1){ //both are linear rotors
//      cnt = edg1 * edg2 /(r1[0] * r2[0] * sym1 * sym2);
//    }
//    else if (r1.size() == 1){ // the first molecule is a linear rotor, but the second is not
//      cnt = edg1 * edg2 * sqrt(M_PI/(r2[0]*r2[1]*r2[2]))/sym;
//      cnt /= r1[0];
//    }
//    else if (r2.size() == 1){ // the second molecule is a linear rotor, but the first is not
//      cnt = edg1 * edg2 * sqrt(M_PI/(r1[0]*r1[1]*r1[2]))/sym;
//      cnt /= r2[0];
//    }
//    else{ //both are not linear rotors
//
//    }
//  }
//  else{ //the second reactant has only one nucleus, needs to treat the first reactant only.
//  }
//
//
//  mfq = f3.size();
//
//  for (int i = 0; i < mfq; ++i) {
//    n = frq[i];
//    for (j=1;j<=ncll-n;j++) tn[n+j] += tn[j];
//  }
//
//  return true;
//}
