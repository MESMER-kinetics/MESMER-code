//-------------------------------------------------------------------------------------------
//
// logExp.cpp
//
// Author: Chi-Hsiu Liang
// Date:   01/Jul/2009
//
//-------------------------------------------------------------------------------------------
#include "logExp.h"
#include <stdexcept>

using namespace std;

namespace logExpGroup
{
  const double logExp::Tinfinity = std::numeric_limits<double>::infinity();

  const int logExp::log10Limit = int(-log10(numeric_limits<double>::epsilon())+1.0);

  void logExp::LPAddition      (const logExp& b, const bool& isAdd){

    /*
    Addition of two numbers of different magnitudes is not a small task consider if we need to do this very
    many times.
    In this class, we will always do addition in the magnitude of the larger number.
    Looks a bit lengthy, but in most of the cases it should run quick... I hope.
    */

    static unsigned int countAdd(0);
    ++countAdd;


#ifdef LogExpDebugMode
    const double tempDoubleValue(isAdd ? doubleValue + b.doubleValue : doubleValue - b.doubleValue);
#endif

    if (b.logValue == -Tinfinity){ // b equals to zero
      // do nothing
    }
    else if (logValue == -Tinfinity){ // *this equals to zero
      logValue  = b.logValue;
      rank      = b.rank;
      boolSign  = (b.boolSign == isAdd);
    }
    else{ // *this and b are non-zero
      int diffLog = int(rank - b.rank);
      if (abs(diffLog) >= log10Limit){ // addition of very different magnitudes numbers (keep the greater magnitude one)
        if (diffLog>0){ // the magnitude of the left number is much greater than the right
          // do nothing, as the magnitudes are too different to be accounted by the current numerical type double
        }
        else{
          // copy the value of b into *this
          logValue  = b.logValue;
          rank      = b.rank;
          boolSign  = (b.boolSign == isAdd);
        }
      }
      else if (diffLog > 0){ // a's rank is greater than b
        double fvb = std::pow(10.,b.logValue);
        double fva = std::pow(10.,logValue);
        double fvc(0.0);

        // Use the boolSign of a, so no need to change a's boolSign.
        // Whether it is deduction or addition
        // would depend on whether boolSign == (boolSign == isAdd)
        if (boolSign == (b.boolSign == isAdd)){
          fvc = fva + fvb * std::pow(10.0, -diffLog);
        }
        else{
          fvc = fva - fvb * std::pow(10.0, -diffLog);
        }
        // check the magnitude of the number after addition or deduction
        checkAndAssign(fvc);
      }
      else if (diffLog == 0){ // addition or deduction of two numbers with equal magnitude
        double fvc(0.0);
        bool tempBool(b.boolSign == isAdd);
        if (boolSign == tempBool){ // addition
          checkAndAssign(std::pow(10.0, logValue) + std::pow(10.0, b.logValue));
        }
        else{ // deduction: need to compare the magnitude
          fvc = std::pow(10.0, logValue) - std::pow(10.0, b.logValue);
          if (fvc == 0.0){
            logValue = -Tinfinity; rank = 0; boolSign = true;
          }
          else{
            if (fvc < 0.0){
              boolSign = (b.boolSign == isAdd);
            }
            // check the magnitude of the number after addition or deduction
            checkAndAssign(std::abs(fvc));
          }
        }
      }
      else{  // b is greater than a (use the sign (boolSign && isAdd))
        rank = b.rank;
        double fvb = std::pow(10.0, b.logValue);
        double fva = std::pow(10.0, logValue);
        double fvc(0.0);

        // Use the boolSign of (b.boolSign == isAdd).
        bool tempBool(b.boolSign == isAdd);
        // Whether it is addition or deduction
        // would depend on whether boolSign == (b.boolSign == isAdd)
        if (boolSign == tempBool){ // ADDING VALUES
            fvc =   fva * std::pow(10.0, diffLog) + fvb;
        }
        else{                                   // DEDUCTING VALUES
            fvc =  -fva * std::pow(10.0, diffLog) + fvb;
        }

        // assign the boolean value in the end.
        boolSign = tempBool;

        // check the magnitude of the number after addition or deduction
        checkAndAssign(fvc);
      }
    }


#ifdef LogExpDebugMode
    convValue = to_Type(*this);
    const double tempD((tempDoubleValue - convValue) / tempDoubleValue);
    if (std::abs(tempD) > 1.0E-2){
      const double tempE(tempD);
    }
    doubleValue = tempDoubleValue;
#endif

  };


  void logExp::LPMultiplication(const logExp& b, const bool& isMul){
    // A multiplication process involves additions of the log values of both ranks, whereas the division
    // is simply the opposite.

    static unsigned int countMul(0);

    ++countMul;

#ifdef LogExpDebugMode
    doubleValue = isMul ? doubleValue * b.doubleValue : doubleValue / b.doubleValue;
#endif

    if (logValue != -Tinfinity){ // if a == 0, there is no point to do this multiplication.
      if (b.logValue == -Tinfinity){ // b == 0
        if (isMul == false){
          throw (std::runtime_error("Error dividing a number with a zero."));
        }
        else{ // a multiply with zero
          logValue = -Tinfinity;
          rank = 0;
          boolSign = 0;
        }
      }
      else{ // b != 0
        if (logValue != -Tinfinity){
          boolSign = (boolSign == b.boolSign);
          if (isMul){
            rank += b.rank;
            logValue += b.logValue;
          }
          else{
            rank -= b.rank;
            logValue -= b.logValue;
          }

          if (logValue < 0.0 || logValue >= 1.0){
            const double flrLogV(floor(logValue));
            logValue -= flrLogV;
            rank += int(flrLogV);
          }
        }
      }
    }


#ifdef LogExpDebugMode
    convValue = to_Type(*this);

    double tempD(0.0);
    if ((convValue != doubleValue) && doubleValue != 0.0){
      tempD = (doubleValue - convValue) / doubleValue;
    }
    if (std::abs(tempD) > 1.0E-2){
      const double tempE(tempD);
    }
#endif

  };

  // Conversion to double
  double to_Type(const logExp& a){
    const double temp(a.getBoolSign() ?  std::pow(10.,a.getLogValue()+double(a.getRank()))
                                      : -std::pow(10.,a.getLogValue()+double(a.getRank())));
    return temp;
  }

  logExp pow(const logExp& a, int p){
    logExp b(a);
    double rmd(b.getLogValue() * p);
    int quo(b.getRank() * p);
    double flrRmd(floor(rmd));
    b.setRank(quo + int(flrRmd));
    b.setLogValue(rmd - flrRmd);
    return b;
  }

  // EXP() definition
  logExp exp(const logExp& a){
    /*Stratagy:
    Use the power series to calculate the exp. This method is good for small
    numbers, as for big numbers the accuracy is achieved by expanding the power
    series to infinity. However, in our application we aim to make approximation
    to small numbers.
    e^x = 1 + x + x^2 / 2! + x^3 / 3! + x^4 / 4! + ...
    */
    logExp b(a+1.0);

    /*
    Need to find a way to make sure the value is converged while not exceeding
    the numerical limit. For now we only do to 10 to serve the purpose;
    */
    logExp XDenominator(1.0);
    for (int i (2); i < 10; ++i){
      XDenominator *= double(i);
      b += (pow(a, i) / XDenominator);
    }

    return b;
  }

  // ABS() definitions
  logExp abs(const logExp& a){
    logExp b(a);
    b.setBoolSign(true);

#ifdef LogExpDebugMode
    b.doubleValue = std::abs(a.doubleValue);
#endif

    return b;
  };

  // SQRT() definition
  logExp sqrt(const logExp& a){
    logExp b(a);

#ifdef LogExpDebugMode
    b.doubleValue = std::sqrt(a.doubleValue);
#endif

    if (b.getBoolSign()){
      const double hLogTV((b.getLogValue() + double(b.getRank())) / 2.0);
      const double fhLogTV(floor(hLogTV));
      b.setLogValue(hLogTV-fhLogTV);
      b.setRank(int(fhLogTV));
    }
    else{
      throw (std::runtime_error("Error taking square root of a negative number."));
    }

#ifdef LogExpDebugMode
    b.convValue = to_Type(b);
#endif

    return b;
  };

  void append_expn(std::string &str, int expn) {
    int k;

    str += (expn < 0 ? '-' : '+');
    expn = std::abs(expn);

    if (expn >= 100) {
      k = (expn / 100);
      str += '0' + k;
      expn -= 100*k;
    }

    k = (expn / 10);
    str += '0' + k;
    expn -= 10*k;

    str += '0' + expn;
  }


  //ostream operator
  std::ostream& operator<<(std::ostream &os, const logExp &a){
    bool showpos = (os.flags() & std::ios_base::showpos) != 0;
    bool uppercase =  (os.flags() & std::ios_base::uppercase) != 0;
    std::ios_base::fmtflags float_field = os.flags() & std::ios_base::floatfield;
    std::ios_base::fmtflags adjust_field = os.flags() & std::ios_base::adjustfield;

  return os << a.write(os.precision(), os.width(), float_field,
      adjust_field, showpos, uppercase, os.fill());
  }

  //Binary operators
  logExp operator/(double a, const logExp &b) { return logExp(a) / b; }
  logExp operator*(double a, const logExp &b) { return logExp(a) * b; }

  //Comparison operators: compare with type double on the right hand side
  bool operator==(const logExp& a, const double& b){ return (to_Type(a) == b); }
  bool operator!=(const logExp& a, const double& b){ return (to_Type(a) != b); }
  bool operator< (const logExp& a, const double& b){ return (to_Type(a) <  b); }
  bool operator> (const logExp& a, const double& b){ return (to_Type(a) >  b); }
  bool operator>=(const logExp& a, const double& b){ return (to_Type(a) >= b); }
  bool operator<=(const logExp& a, const double& b){ return (to_Type(a) <= b); }

  //Comparison operators: compare with type double on the left hand side
  bool operator==(const double& a, const logExp& b){ return (a == to_Type(b)); }
  bool operator!=(const double& a, const logExp& b){ return (a != to_Type(b)); }
  bool operator< (const double& a, const logExp& b){ return (a <  to_Type(b)); }
  bool operator> (const double& a, const logExp& b){ return (a >  to_Type(b)); }
  bool operator>=(const double& a, const logExp& b){ return (a >= to_Type(b)); }
  bool operator<=(const double& a, const logExp& b){ return (a <= to_Type(b)); }

    //Comparison operators: compare with type double on the right hand side
  bool operator==(const logExp& a, const logExp& b){
    return (to_Type(a) == to_Type(b));
  }
  bool operator!=(const logExp& a, const logExp& b){ return (to_Type(a) != to_Type(b)); }
  bool operator< (const logExp& a, const logExp& b){ return (to_Type(a) <  to_Type(b)); }
  bool operator> (const logExp& a, const logExp& b){ return (to_Type(a) >  to_Type(b)); }
  bool operator>=(const logExp& a, const logExp& b){ return (to_Type(a) >= to_Type(b)); }
  bool operator<=(const logExp& a, const logExp& b){ return (to_Type(a) <= to_Type(b)); }


}//namespace

