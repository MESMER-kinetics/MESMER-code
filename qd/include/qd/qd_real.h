/*
 * include/qd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Quad-double precision (>= 212-bit significand) floating point arithmetic
 * package, written in ANSI C++, taking full advantage of operator overloading.
 * Uses similar techniques as that of David Bailey's double-double package 
 * and that of Jonathan Shewchuk's adaptive precision floating point 
 * arithmetic package.  See
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   http://www.cs.cmu.edu/~quake/robust.html
 *
 * for more details.
 *
 * Yozo Hida
 */
#ifndef _QD_QD_REAL_H
#define _QD_QD_REAL_H

#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/dd_real.h>
#include <Eigen/Dense>

class QD_API qd_real {
protected:

  /* Eliminates any zeros in the middle component(s). */
  void zero_elim();
  void zero_elim(double &e);

  void renorm();
  void renorm(double &e);

  void quick_accum(double d, double &e);
  void quick_prod_accum(double a, double b, double &e);

  /* Computes  qd * d  where d is known to be a power of 2.
     This can be done component wise.                      */
  friend QD_API qd_real mul_pwr2(const qd_real &qd, double d);

  /* Sin / Cos Tables */
  static const qd_real _pi1024;  /* Used in sin / cos routine. */
  static const qd_real sin_table[];
  static const qd_real cos_table[];

public:
  double x[4];    /* The Components. */

  /* Protected Constructor. */
  qd_real(double x0, double x1, double x2, double x3);
  qd_real(const double *xx);

  /* Useful constants. */
  static const qd_real _2pi;
  static const qd_real _pi;
  static const qd_real _3pi4;
  static const qd_real _pi2;
  static const qd_real _pi4;
  static const qd_real _pi8;
  static const qd_real _e;
  static const qd_real _log2;
  static const qd_real _log10;
  static const qd_real _nan;
  static const qd_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const qd_real _max;
  static const qd_real _safe_max;
  static const int _ndigits;

  /* Constructors */
  qd_real();
  qd_real(const char *s);
  qd_real(const dd_real &dd);
  qd_real(double d);
  qd_real(int i);

  /* Member Access */
  double operator[](int i) const;
  double &operator[](int i);

  static void abort(const char *msg);

  bool isnan() const;
  bool isfinite() const { return QD_ISFINITE(x[0])!=0; }
  bool isinf() const { return QD_ISINF(x[0]); }

  /* Addition */
  friend QD_API qd_real operator+(const qd_real &a); //CM
  friend QD_API qd_real operator+(const qd_real &a, const qd_real &b);
  friend QD_API qd_real operator+(const dd_real &a, const qd_real &b);
  friend QD_API qd_real operator+(const qd_real &a, const dd_real &b);
  friend QD_API qd_real operator+(const qd_real &a, double b);
  friend QD_API qd_real operator+(double a, const qd_real &b);
  static qd_real ieee_add(const qd_real &a, const qd_real &b);
  static qd_real sloppy_add(const qd_real &a, const qd_real &b);

  /* Self-Addition */
  qd_real &operator+=(double a);
  qd_real &operator+=(const dd_real &a);
  qd_real &operator+=(const qd_real &a);

  /* Subtraction */
  friend QD_API qd_real operator-(const qd_real &a, const qd_real &b);
  friend QD_API qd_real operator-(const dd_real &a, const qd_real &b);
  friend QD_API qd_real operator-(const qd_real &a, const dd_real &b);
  friend QD_API qd_real operator-(const qd_real &a, double b);
  friend QD_API qd_real operator-(double a, const qd_real &b);

  /* Self-Subtraction */
  qd_real &operator-=(double a);
  qd_real &operator-=(const dd_real &a);
  qd_real &operator-=(const qd_real &a);

  /* Multiplication */
  friend QD_API qd_real operator*(const qd_real &a, const qd_real &b);
  friend QD_API qd_real operator*(const dd_real &a, const qd_real &b);
  friend QD_API qd_real operator*(const qd_real &a, const dd_real &b);
  friend QD_API qd_real operator*(const qd_real &a, double b);
  friend QD_API qd_real operator*(double a, const qd_real &b);
  static qd_real sloppy_mul(const qd_real &a, const qd_real &b);
  static qd_real accurate_mul(const qd_real &a, const qd_real &b);

  /* Self-Multiplication */
  qd_real &operator*=(double a);
  qd_real &operator*=(const dd_real &a);
  qd_real &operator*=(const qd_real &a);

  /* Division */
  friend QD_API qd_real operator/(const qd_real &a, const qd_real &b);
  friend QD_API qd_real operator/(const dd_real &a, const qd_real &b);
  friend QD_API qd_real operator/(const qd_real &a, const dd_real &b);
  friend QD_API qd_real operator/(const qd_real &a, double b);
  friend QD_API qd_real operator/(double a, const qd_real &b);
  static qd_real sloppy_div(const qd_real &a, const dd_real &b);
  static qd_real accurate_div(const qd_real &a, const dd_real &b);
  static qd_real sloppy_div(const qd_real &a, const qd_real &b);
  static qd_real accurate_div(const qd_real &a, const qd_real &b);

  /* Self-Division */
  qd_real &operator/=(double a);
  qd_real &operator/=(const dd_real &a);
  qd_real &operator/=(const qd_real &a);

  /* Square, Square Root, Power, N-th Root. */
  friend QD_API qd_real sqr(const qd_real &a);
  friend QD_API qd_real sqrt(const qd_real &a);
  friend QD_API qd_real pow(const qd_real &a, int n);
  friend QD_API qd_real pow(const qd_real &a, const qd_real &b);
  friend QD_API qd_real npwr(const qd_real &a, int n);
  qd_real operator^(int n) const;
  friend QD_API qd_real nroot(const qd_real &a, int n);

  /* Unary Minus */
  qd_real operator-() const;

  /* Remainder */
  friend QD_API qd_real rem(const qd_real &a, const qd_real &b);
  friend QD_API qd_real drem(const qd_real &a, const qd_real &b);
  friend QD_API qd_real divrem(const qd_real &a, const qd_real &b, qd_real &r);

  /* Assignment */
  qd_real &operator=(double a);
  qd_real &operator=(int i);
  qd_real &operator=(const dd_real &a);
  qd_real &operator=(const char *s);

  /* Cast */
  friend dd_real to_dd_real(const qd_real &a);
  friend double  to_double(const qd_real &a);
  friend int     to_int(const qd_real &a);
  explicit operator const double() const; //CM from Eigen forum

  /* Equality Comparison */
  friend QD_API bool operator==(const qd_real &a, const qd_real &b);
  friend QD_API bool operator==(const qd_real &a, const dd_real &b);
  friend QD_API bool operator==(const dd_real &a, const qd_real &b);
  friend QD_API bool operator==(double a, const qd_real &b);
  friend QD_API bool operator==(const qd_real &a, double b);

  /* Less-Than Comparison */
  friend QD_API bool operator<(const qd_real &a, const qd_real &b);
  friend QD_API bool operator<(const qd_real &a, const dd_real &b);
  friend QD_API bool operator<(const dd_real &a, const qd_real &b);
  friend QD_API bool operator<(double a, const qd_real &b);
  friend QD_API bool operator<(const qd_real &a, double b);

  /* Greater-Than Comparison */
  friend QD_API bool operator>(const qd_real &a, const qd_real &b);
  friend QD_API bool operator>(const qd_real &a, const dd_real &b);
  friend QD_API bool operator>(const dd_real &a, const qd_real &b);
  friend QD_API bool operator>(double a, const qd_real &b);
  friend QD_API bool operator>(const qd_real &a, double b);

  /* Less-Than-Or-Equal-To Comparison */
  friend QD_API bool operator<=(const qd_real &a, const qd_real &b);
  friend QD_API bool operator<=(const qd_real &a, const dd_real &b);
  friend QD_API bool operator<=(const dd_real &a, const qd_real &b);
  friend QD_API bool operator<=(double a, const qd_real &b);
  friend QD_API bool operator<=(const qd_real &a, double b);

  /* Greater-Than-Or-Equal-To Comparison */
  friend QD_API bool operator>=(const qd_real &a, const qd_real &b);
  friend QD_API bool operator>=(const qd_real &a, const dd_real &b);
  friend QD_API bool operator>=(const dd_real &a, const qd_real &b);
  friend QD_API bool operator>=(double a, const qd_real &b);
  friend QD_API bool operator>=(const qd_real &a, double b);

  /* Not-Equal-To Comparison */
  friend QD_API bool operator!=(const qd_real &a, const qd_real &b);
  friend QD_API bool operator!=(const qd_real &a, const dd_real &b);
  friend QD_API bool operator!=(const dd_real &a, const qd_real &b);
  friend QD_API bool operator!=(double a, const qd_real &b);
  friend QD_API bool operator!=(const qd_real &a, double b);

  /* Other Comparisons.  These are faster than directly
     comparing to 0 or 1.                               */
  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  /* Micellaneous algebraic operations */
  friend QD_API qd_real fabs(const qd_real &a);
  friend QD_API qd_real abs(const qd_real &a);    /* same as fabs */

  /* Computes  a * 2^n. */
  friend QD_API qd_real ldexp(const qd_real &a, int n);

  /* Rounding */
  friend QD_API qd_real nint(const qd_real &a);
  friend QD_API qd_real quick_nint(const qd_real &a);
  friend QD_API qd_real floor(const qd_real &a);
  friend QD_API qd_real ceil(const qd_real &a);
  friend QD_API qd_real aint(const qd_real &a);

  /* Trigonometric Functions */
  friend QD_API qd_real sin(const qd_real &a);
  friend QD_API qd_real cos(const qd_real &a);
  friend QD_API qd_real tan(const qd_real &a);
  friend QD_API void sincos(const qd_real &a, qd_real &s, qd_real &c);

  /* Inverse Trigonometric Functions */
  friend QD_API qd_real asin(const qd_real &a);
  friend QD_API qd_real acos(const qd_real &a);
  friend QD_API qd_real atan(const qd_real &a);
  friend QD_API qd_real atan2(const qd_real &y, const qd_real &x);

  /* Exponential / Logarithm */
  friend QD_API qd_real exp(const qd_real &a);
  friend QD_API qd_real log(const qd_real &a);
  friend QD_API qd_real log10(const qd_real &a);

  /* Hyperbolic Functions */
  friend QD_API qd_real sinh(const qd_real &a);
  friend QD_API qd_real cosh(const qd_real &a);
  friend QD_API qd_real tanh(const qd_real &a);
  friend QD_API void sincosh(const qd_real &a, qd_real &sin_qd, qd_real &cos_qd);

  /* Inverse Hyperbolic Functions */
  friend QD_API qd_real asinh(const qd_real &a);
  friend QD_API qd_real acosh(const qd_real &a);
  friend QD_API qd_real atanh(const qd_real &a);

  /* Random Number */
  static qd_real rand(void);
  friend QD_API qd_real qdrand(void);

  /* Min / Max */
  friend QD_API qd_real max(const qd_real &a, const qd_real &b);
  friend QD_API qd_real max(const qd_real &a, const qd_real &b, const qd_real &c);
  friend QD_API qd_real min(const qd_real &a, const qd_real &b);
  friend QD_API qd_real min(const qd_real &a, const qd_real &b, const qd_real &c);

  /* Input / Output */
  friend QD_API std::ostream &operator<<(std::ostream &s, const qd_real &a);
  friend QD_API std::istream &operator>>(std::istream &s, qd_real &a);
  
  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;  
      /* Note: s must hold d+8 chars. */
  std::string write(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags floatfield = static_cast<std::ios_base::fmtflags>(0), 
      std::ios_base::fmtflags adjustfield = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  static int read(const char *s, qd_real &a);

  /* Debugging methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static qd_real debug_rand();

};

namespace std {
  template <>
  class numeric_limits<qd_real> : public numeric_limits<double> {
  public:
    inline static double epsilon() { return qd_real::_eps; }
    inline static double min() { return qd_real::_min_normalized; }
    inline static qd_real max() { return qd_real::_max; }
    inline static qd_real safe_max() { return qd_real::_safe_max; }
    static const int digits = 209;
    static const int digits10 = 62;
  };
}

QD_API qd_real polyeval(const qd_real *c, int n, const qd_real &x);
QD_API qd_real polyroot(const qd_real *c, int n, 
    const qd_real &x0, int max_iter = 64, double thresh = 0.0);

QD_API qd_real qdrand(void);
//QD_API qd_real sqrt(const qd_real &a); CM duplicate?

QD_API inline bool isnan(const qd_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const qd_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const qd_real &a) { return a.isinf(); }


#ifdef QD_INLINE
#include <qd/qd_inline.h>
#endif

namespace Eigen {
  //http://eigen.tuxfamily.org/dox/NumTraits_8h_source.html
  template<> struct NumTraits<qd_real> : GenericNumTraits<qd_real>
  {
    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 10,
      AddCost = 10,
      MulCost = 40
    };

    /*
    inline static Real epsilon(long Precision = mpfr::mpreal::get_default_prec()) { return mpfr::machine_epsilon(Precision); }
    inline static Real epsilon(const Real& x) { return mpfr::machine_epsilon(x); }
    */

    static inline qd_real dummy_precision()
    {
      // make sure to override this for floating-point types
      return dd_real(1.21543267145725e-63);
    }

  };
}
#endif /* _QD_QD_REAL_H */

