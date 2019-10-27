/*
 * tests/qd_test.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This contains some simple tests to sanity check the double-double
 * and quad-double library.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <qd/qd_real.h>
#include "../error.h"

using std::cerr;
using std::endl;

using std::abs;
using std::sqrt;
using std::strcmp;
using std::exit;

using mesmer::ctest;

// Global flags passed to the main program.
static bool flag_test_dd = false;
static bool flag_test_qd = false;
bool flag_verbose = true;

bool print_result(bool result) {
  if (result)
    ctest << "Test passed." << endl;
  else
    ctest << "Test FAILED." << endl;
  return result;
}

template <class T>
class TestSuite {
  static const int double_digits;
public:
  bool test1();
  bool test2();
  bool test3();
  bool test4();
  bool test5();
  bool test6();
  bool testall();
};

template <class T>
const int TestSuite<T>::double_digits = 6;

/* Test 1.   Polynomial Evaluation / Polynomial Solving */
template <class T>
bool TestSuite<T>::test1() {
  ctest << endl;
  ctest << "Test 1.  (Polynomial)." << endl;

  static const int n = 8;
  T *c = new T[n];
  T x, y;

  for (int i = 0; i < n; i++)
    c[i] = static_cast<double>(i+1);

  x = polyroot(c, n-1, T(0.0));
  y = polyeval(c, n-1, x);

  if (flag_verbose) {
    ctest.precision(T::_ndigits);
    ctest << "Root Found:  x  = " << x << endl;
    ctest << "           p(x) = " << y << endl;
  }

  delete [] c;
  return (to_double(y) < 4.0 * T::_eps);
}

/* Test 2.  Machin's Formula for Pi. */
template <class T>
bool TestSuite<T>::test2() {

  ctest << endl;
  ctest << "Test 2.  (Machin's Formula for Pi)." << endl;
  
  /* Use the Machin's arctangent formula:

       pi / 4  =  4 arctan(1/5) - arctan(1/239)

     The arctangent is computed based on the Taylor series expansion

       arctan(x) = x - x^3 / 3 + x^5 / 5 - x^7 / 7 + ...
  */

  T s1, s2, t, r;
  int k;
  int sign;
  double d;
  double err;

  /* Compute arctan(1/5) */
  d = 1.0;
  t = T(1.0) / 5.0;
  r = sqr(t);
  s1 = 0.0;
  k = 0;

  sign = 1;
  while (t > T::_eps) {
    k++;
    if (sign < 0)
      s1 -= (t / d);
    else
      s1 += (t / d);

    d += 2.0;
    t *= r;
    sign = -sign;
  }

  if (flag_verbose)
    ctest << k << " Iterations" << endl;

  /* Compute arctan(1/239) */
  d = 1.0;
  t = T(1.0) / 239.0;
  r = sqr(t);
  s2 = 0.0;
  k = 0;

  sign = 1;
  while (t > T::_eps) {
    k++;
    if (sign < 0)
      s2 -= (t / d);
    else
      s2 += (t / d);

    d += 2.0;
    t *= r;
    sign = -sign;
  }

  if (flag_verbose)
    ctest << k << " Iterations" << endl;

  T p = 4.0 * s1 - s2;

  p *= 4.0;
  err = abs(to_double(p - T::_pi));

  if (flag_verbose) {
    ctest.precision(T::_ndigits);
    ctest << "   pi = " << p << endl;
    ctest << "  _pi = " << T::_pi << endl;

    ctest.precision(double_digits);
    ctest << "error = " << err << " = " << err / T::_eps << " eps" << endl;
  }

  return (err < 8.0 * T::_eps);
}

/* Test 3.  Salamin-Brent Quadratic Formula for Pi. */
template <class T>
bool TestSuite<T>::test3() {
  ctest << endl;
  ctest << "Test 3.  (Salamin-Brent Quadratic Formula for Pi)." << endl;
  ctest.precision(T::_ndigits);

  T a, b, s, p;
  T a_new, b_new, p_old;
  double m;
  double err;
  const int max_iter = 20;

  a = 1.0;
  b = sqrt(T(0.5));
  s = 0.5;
  m = 1.0;

  p = 2.0 * sqr(a) / s;
  if (flag_verbose)
    ctest << "Iteration  0: " << p << endl;
  for (int i = 1; i <= max_iter; i++) {
    m *= 2.0;
    a_new = 0.5 * (a + b);
    b_new = a * b;
    s -= m * (sqr(a_new) - b_new);
    a = a_new;
    b = sqrt(b_new);
    p_old = p;
    p = 2.0 * sqr(a) / s;
    if (flag_verbose)
      ctest << "Iteration " << std::setw(2) << i << ": " << p << endl;
    if (abs(to_double(p - p_old)) < 64 * T::_eps)
      break;
  }

  err = abs(to_double(p - T::_pi));

  if (flag_verbose) {
    ctest << "         _pi: " << T::_pi << endl;
    ctest.precision(double_digits);
    ctest << "       error: " << err << " = " << err / T::_eps << " eps" << endl;
  }

  // for some reason, this test gives relatively large error compared
  // to other tests.  May need to be looked at more closely.
  return (err < 1024.0 * T::_eps);
}

/* Test 4.  Borwein Quartic Formula for Pi. */
template <class T>
bool TestSuite<T>::test4() {
  ctest << endl;
  ctest << "Test 4.  (Borwein Quartic Formula for Pi)." << endl;
  ctest.precision(T::_ndigits);

  T a, y, p, r, p_old;
  double m;
  double err;
  const int max_iter = 20;

  a = 6.0 - 4.0 * sqrt(T(2.0));
  y = sqrt(T(2.0)) - 1.0;
  m = 2.0;

  p = 1.0 / a;
  if (flag_verbose)
    ctest << "Iteration  0: " << p << endl;

  for (int i = 1; i <= max_iter; i++) {
    m *= 4.0;
    r = nroot(1.0 - sqr(sqr(y)), 4);
    y = (1.0 - r) / (1.0 + r);
    a = a * sqr(sqr(1.0 + y)) - m * y * (1.0 + y + sqr(y));
    
    p_old = p;
    p = 1.0 / a;
    if (flag_verbose)
      ctest << "Iteration " << std::setw(2) << i << ": " << p << endl;
    if (abs(to_double(p - p_old)) < 16 * T::_eps)
      break;
  }

  err = abs(to_double(p - T::_pi));
  if (flag_verbose) {
    ctest << "         _pi: " << T::_pi << endl;
    ctest.precision(double_digits);
    ctest << "       error: " << err << " = " << err / T::_eps << " eps" << endl;
  }  

  return (err < 256.0 * T::_eps);
}

/* Test 5.  Taylor Series Formula for E. */
template <class T>
bool TestSuite<T>::test5() {

  ctest << endl;
  ctest << "Test 5.  (Taylor Series Formula for E)." << endl;
  ctest.precision(T::_ndigits);

  /* Use Taylor series

       e = 1 + 1 + 1/2! + 1/3! + 1/4! + ...

     To compute e.
  */

  T s = 2.0, t = 1.0;
  double n = 1.0;
  double delta;
  int i = 0;

  while (t > T::_eps) {
    i++;
    n += 1.0;
    t /= n;
    s += t;
  }

  delta = abs(to_double(s - T::_e));

  if (flag_verbose) {
    ctest << "    e = " << s << endl;
    ctest << "   _e = " << T::_e << endl;

    ctest.precision(double_digits);
    ctest << "error = " << delta << " = " << delta / T::_eps << " eps" << endl;
    ctest << i << " iterations." << endl;
  }

  return (delta < 64.0 * T::_eps);
}

/* Test 6.  Taylor Series Formula for log 2.*/
template <class T>
bool TestSuite<T>::test6() {
  ctest << endl;
  ctest << "Test 6.  (Taylor Series Formula for Log 2)." << endl;
  ctest.precision(T::_ndigits);

  /* Use the Taylor series

      -log(1-x) = x + x^2/2 + x^3/3 + x^4/4 + ...

     with x = 1/2 to get  log(1/2) = -log 2.
  */

  T s = 0.5;
  T t = 0.5;
  double delta;
  double n = 1.0;
  double i = 0;

  while (abs(t) > T::_eps) {
    i++;
    n += 1.0;
    t *= 0.5;
    s += (t/n);
  }

  delta = abs(to_double(s - T::_log2));

  if (flag_verbose) {
    ctest << " log2 = " << s << endl;
    ctest << "_log2 = " << T::_log2 << endl;

    ctest.precision(double_digits);
    ctest << "error = " << delta << " = " << (delta / T::_eps) 
         << " eps" << endl;
    ctest << i << " iterations." << endl;
  }

  return (delta < 4.0 * T::_eps);
}

template <class T>
bool TestSuite<T>::testall() {
  bool pass = true;
  pass &= print_result(test1());
  pass &= print_result(test2());
  pass &= print_result(test3());
  pass &= print_result(test4());
  pass &= print_result(test5());
  pass &= print_result(test6());
  return pass;
}

