//-------------------------------------------------------------------------------------------
//
// logExp.h
//
// Author: Chi-Hsiu Liang
// Date:   01/Jul/2009
//
//-------------------------------------------------------------------------------------------
#ifndef GUARD_logExp_h
#define GUARD_logExp_h

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include "MesmerPrecision.h"

#include <map>
#include "oberror.h"


//#define LogExpDebugMode true

namespace logExpGroup
{

/*
This class is here mainly to resolve the numerical issues triggerd by extreme conditions (P/double) in
master equation calculations. The class was built on the belief that numbers are more important on
their magnitudes than the precision in most of the master equation simulations.

This class is capable of doing numerical range practically from 10^(-IntRange) to 10^(IntRange),
where IntRange is the maximum number of int and is thus platform dependent. The numerical precision
(number of digits) of this class depends on the numerical type used. In MESMER, the type options
are double, double-double and quad-double.
*/

  void append_expn(std::string &str, int expn);


  class logExp
  {
    /*
    Class logExp is a class getting values stored in log scale.

    Some static variables has to be put in the diagonalization routine to capture the number of
    deductions and additions that are out of numerical precision. For example, an addition
    of a big number and a small number resulting an unchanged big number. And, an deduction
    of two very close non-zero numbers.
    */

  private:

    static const double Tinfinity;
    static const int log10Limit;

    double    logValue;
    int  rank;
    bool boolSign;


  public:


#ifdef LogExpDebugMode
    double convValue; // value converted by this method
    double doubleValue; // double precision calculated value
#endif

    //Initializer operator from double
    logExp(double val=0.0)
    {
      *this = val;
    }

    //Assignment operator
    logExp& operator=(const double& val){
      if (val == 0.0){
        // basically the sum of rank and logValue is negative infinity, and the exponential of it
        // is zero.
        logValue = -Tinfinity;
        rank     = 0;
        boolSign = true;
      }
      else{
        //first need to know the magnitude of the value
        const double logV(log10(val >= 0.0 ? val : -val));
        const double flrLogV(floor(logV));
        rank = int(flrLogV);
        boolSign=(val > 0.0);
        logValue=logV-flrLogV;
      }

#ifdef LogExpDebugMode
      doubleValue = val;
      convValue = (boolSign ?  pow(10.,logValue+double(rank))
                            : -pow(10.,logValue+double(rank)));
#endif

      return *this;
    }

    //Unary minus
    logExp operator-() const {
       logExp a(*this);
       a.setBoolSign(a.getBoolSign() ? false : true);

#ifdef LogExpDebugMode
      a.doubleValue = -(this->doubleValue);
      a.convValue = (boolSign ?  pow(10.,logValue+double(rank))
                              : -pow(10.,logValue+double(rank)));
#endif
       return a;
    }

    //Compound Assignment Operators
    logExp& operator= (const logExp& r) {
      logValue = r.logValue; rank = r.rank; boolSign = r.boolSign;
#ifdef LogExpDebugMode
      doubleValue = r.doubleValue;
      convValue = r.convValue;
#endif
      return *this;
    }
    logExp& operator+=(const logExp& r) {
      LPAddition(r, true);         return *this;
    }
    logExp& operator-=(const logExp& r) { LPAddition(r, false);        return *this; }
    logExp& operator*=(const logExp& r) { LPMultiplication(r, true);   return *this; }
    logExp& operator/=(const logExp& r) {
      LPMultiplication(r, false);  return *this;
    }

    //Binary arithmetic operators
    const logExp operator+(const logExp& other) const { return logExp(*this) += other;}
    const logExp operator-(const logExp& other) const { return logExp(*this) -= other;}
    const logExp operator*(const logExp& other) const {
      return logExp(*this) *= other;
    }
    const logExp operator/(const logExp& other) const { return logExp(*this) /= other;}

    //Compound Assignment Operators with type double
    logExp& operator+=(const double& r) { LPAddition(logExp(r), true);         return *this; }
    logExp& operator-=(const double& r) { LPAddition(logExp(r), false);        return *this; }
    logExp& operator*=(const double& r) { LPMultiplication(logExp(r), true);   return *this; }
    logExp& operator/=(const double& r) { LPMultiplication(logExp(r), false);  return *this; }

    //Binary arithmetic operators with type double
    const logExp operator+(const double& other) const { return logExp(*this) += other;}
    const logExp operator-(const double& other) const { return logExp(*this) -= other;}
    const logExp operator*(const double& other) const { return logExp(*this) *= other;}
    const logExp operator/(const double& other) const { return logExp(*this) /= other;}

    std::string write(int precision, int width,
      std::ios_base::fmtflags float_field, std::ios_base::fmtflags adjust_field,
      bool showpos, bool uppercase, char fill) const {
      std::string s;

      //-------- This section is copied and modified following similar definition in the QD package
      bool fixed = (float_field & std::ios_base::fixed) != 0;
      bool sgn = true;
      int i, e(rank);

      if (!boolSign)
        s += '-';
      else if (showpos)
        s += '+';
      else
        sgn = false;

      if (logValue == -Tinfinity) {
        /* Zero case */
        s += '0';
        if (precision > 0) {
          s += '.';
          s.append(precision, '0');
        }
      } else {
        /* Non-zero case */
        int off = (fixed ? (1 + rank) : 1);
        int d = precision + off;

        if (fixed && d <= 0) {
          s += '0';
          if (precision > 0) {
            s += '.';
            s.append(precision, '0');
          }
        } else {
          char *t = new char[d+1];
          int j;

          to_digits(t, e, d);

          if (fixed) {
            if (off > 0) {
              for (i = 0; i < off; i++) s += t[i];
              if (precision > 0) {
                s += '.';
                for (j = 0; j < precision; j++, i++) s += t[i];
              }
            } else {
              s += "0.";
              if (off < 0) s.append(-off, '0');
              for (i = 0; i < d; i++) s += t[i];
            }
          } else {
            s += t[0];
            if (precision > 0) s += '.';

            for (i = 1; i <= precision; i++)
              s += t[i];

            delete [] t;
          }
        }
      }

      if (!fixed) {
        /* Fill in exponent part */
        s += uppercase ? 'E' : 'e';
        append_expn(s, e);
      }

      /* Fill in the blanks */
      int len = int(s.length());
      if (len < width) {
        int delta = width - len;
        if (adjust_field & std::ios_base::internal) {
          if (sgn)
            s.insert(static_cast<std::string::size_type>(1), delta, fill);
          else
            s.insert(static_cast<std::string::size_type>(0), delta, fill);
        } else if (adjust_field & std::ios_base::left) {
          s.append(delta, fill);
        } else {
          s.insert(static_cast<std::string::size_type>(0), delta, fill);
        }
      }
      //-------- This section is copied and modified following similar definition in the QD package


      return s;
    }

    void to_digits(char *s, int &expn, int precision) const {
      int D = precision + 1;  /* number of digits to compute */

      int i, d;

      /* First determine the (approximate) exponent. */
      double r(pow(10.,logValue));

      /* Extract the digits */
      for (i = 0; i < D; i++) {
        d = static_cast<int>(r);
        r -= d;
        r *= 10.0;

        s[i] = static_cast<char>(d + '0');
      }


      ///* Round, handle carry */
      //if (s[D-1] >= '5') {
      //  s[D-2]++;

      //  i = D-2;
      //  while (i > 0 && s[i] > '9') {
      //    s[i] -= 10;
      //    s[--i]++;
      //  }
      //}

      s[precision] = 0;
    }

    //Modifiers
    void setRank(const int& r){rank = r;}
    void setLogValue(const double& r){logValue = r;}
    void setBoolSign(const bool& r){boolSign = r;}

    //Accessors
    const int getRank(void) const {return rank;}
    const double getLogValue(void) const {return logValue;}
    const bool getBoolSign(void) const {return boolSign;}

  private:

    void LPAddition      (const logExp& b, const bool& isAdd);
    void LPMultiplication(const logExp& b, const bool& isMul);

    void checkAndAssign(const double& floatValue){
      const double tempLogV(log10(floatValue));
      if (floatValue >= 10.0 || floatValue < 1.0){
        const double flrLogV(floor(tempLogV));
        logValue = tempLogV-flrLogV;
        rank += int(flrLogV);
      }
      else{ // assign value directly
        logValue = tempLogV;
      }
    }

  };

  // Conversion to double
  double to_Type(const logExp& a);

  // POW() definition
  logExp pow(const logExp& a, int p);

  // EXP() definition
  logExp exp(const logExp& a);

  // ABS() definition
  logExp abs(const logExp& a);

  // SQRT() definition
  logExp sqrt(const logExp& a);



  //ostream operator
  std::ostream& operator<<(std::ostream &os, const logExp &a);

  //Binary operators
  logExp operator/(double a, const logExp &b);
  logExp operator*(double a, const logExp &b);

  //Comparison operators: compare with type double on the right hand side
  bool operator==(const logExp& a, const double& b);
  bool operator!=(const logExp& a, const double& b);
  bool operator< (const logExp& a, const double& b);
  bool operator> (const logExp& a, const double& b);
  bool operator>=(const logExp& a, const double& b);
  bool operator<=(const logExp& a, const double& b);

  //Comparison operators: compare with type double on the left hand side
  bool operator==(const double& a, const logExp& b);
  bool operator!=(const double& a, const logExp& b);
  bool operator< (const double& a, const logExp& b);
  bool operator> (const double& a, const logExp& b);
  bool operator>=(const double& a, const logExp& b);
  bool operator<=(const double& a, const logExp& b);

    //Comparison operators: compare with type double on the right hand side
  bool operator==(const logExp& a, const logExp& b);
  bool operator!=(const logExp& a, const logExp& b);
  bool operator< (const logExp& a, const logExp& b);
  bool operator> (const logExp& a, const logExp& b);
  bool operator>=(const logExp& a, const logExp& b);
  bool operator<=(const logExp& a, const logExp& b);


}//namespace
#endif

