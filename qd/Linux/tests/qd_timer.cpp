/*
 * tests/qd_timer.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2004
 *
 * Contains code to time basic operations.
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <qd/qd_real.h>
#include <qd/fpu.h>
#include "tictoc.h"

using std::cout;
using std::cerr;
using std::endl;
using std::sqrt;
using std::strcmp;
using std::setw;
using std::setprecision;
using std::fixed;

// Global flags passed to the main program.
static bool flag_test_double = false;
static bool flag_test_dd = false;
static bool flag_test_qd = false;
static bool flag_verbose = false;
static bool flag_long = false;

template <class T>
class TestSuite {
public:
  void test1();
  void test2();
  void test3();
  void test4();
  void test5();
  void test6();
  void testall();
  T pi();
};

template <class T>
T TestSuite<T>::pi() { return T::_pi; }

template <>
double TestSuite<double>::pi() { return 3.141592653589793116; }

template <class T>
void TestSuite<T>::test1() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing addition..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  std::ios_base::fmtflags fmt = cout.flags();
  if (flag_long) n *= 10;

  T a = 1.0 / T(7.0);
  T b = 0.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b += a;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "b = " << b << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   add: ";
  }

  cout << fixed << setprecision(6) << setw(10) << t/n*1.0e6 << " us" << endl;
  cout.flags(fmt);
}

template <class T>
void TestSuite<T>::test2() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing multiplication ..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  if (flag_long) n *= 10;

  T a = 1.0 + 1.0 / T(static_cast<double>(n));
  T b = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b *= a;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "b = " << b << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   mul: ";
  }

  cout << fixed << setprecision(6) << setw(10) << t/n*1.0e6 << " us" << endl;
}

template <class T>
void TestSuite<T>::test3() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing division ..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  if (flag_long) n *= 10;

  T a = 1.0 + 1.0 / T(static_cast<double>(n));
  T b = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b /= a;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "b = " << b << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   div: ";
  }

  cout << fixed << setprecision(6) << setw(10) << t/n*1.0e6 << " us" << endl;
}

template <class T>
void TestSuite<T>::test4() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing square root ..." << endl;
  }

  int n = 10000, i;
  tictoc tv;
  double t;
  if (flag_long) n *= 10;

  T a = 0.0;
  T b = 2.0 + pi();

  tic(&tv);
  for (i = 0; i < n; i++) {
    a = sqrt(a + b);
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "a = " << a << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "  sqrt: ";
  }

  cout << fixed << setprecision(6) << setw(10) << t/n*1.0e6 << " us" << endl;
}

template <class T>
void TestSuite<T>::test5() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing sin ..." << endl;
  }

  int n = 2000, i;
  tictoc tv;
  double t;
  if (flag_long) n *= 10;

  T a = 0.0;
  T b = pi() / static_cast<double>(n);
  T c = 0.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    a += b;
    c += sin(a);
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "a = " << a << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   sin: ";
  }

  cout << fixed << setprecision(6) << setw(10) << t/n*1.0e6 << " us" << endl;
}

template <class T>
void TestSuite<T>::test6() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing log ..." << endl;
  }

  int n = 1000, i;
  tictoc tv;
  double t;
  if (flag_long) n *= 10;

  T a = 0.0;
  T c = exp(T(-50.1));
  T d = exp(T(100.2) / double(n));

  tic(&tv);
  for (i = 0; i < n; i++) {
    a = a + log(c);
    c *= d;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "a = " << a << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   log: ";
  }

  cout << fixed << setprecision(6) << setw(10) << t/n*1.0e6 << " us" << endl;
}

template <class T>
void TestSuite<T>::testall() {
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
}

void print_usage() {
  cout << "qd_test [-h] [-dd] [-qd] [-all]" << endl;
  cout << "  Performs timing tests of the quad-double library." << endl;
  cout << "  By default, double-double and quad-double arithmetics" << endl;
  cout << "  are timed." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -double   Time arithmetic of double." << endl;
  cout << "  -dd       Time arithmetic of double-double." << endl;
  cout << "  -qd       Time arithmetic of quad-double." << endl;
  cout << "  -all      Perform both double-double and quad-double tests." << endl;
  cout << "  -v        Verbose output." << endl;
  cout << "  -long     Perform a longer timing loop." << endl;
}

int main(int argc, char *argv[]) {
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  /* Parse the arguments. */
  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      std::exit(0);
    } else if (strcmp(arg, "-double") == 0) {
      flag_test_double = true;
    } else if (strcmp(arg, "-dd") == 0) {
      flag_test_dd = true;
    } else if (strcmp(arg, "-qd") == 0) {
      flag_test_qd = true;
    } else if (strcmp(arg, "-all") == 0) {
      flag_test_double = flag_test_dd = flag_test_qd = true;
    } else if (strcmp(arg, "-v") == 0) {
      flag_verbose = true;
    } else if (strcmp(arg, "-long") == 0) {
      flag_long = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  /* If no flag, test both double-double and quad-double. */
  if (!flag_test_double && !flag_test_dd && !flag_test_qd) {
    flag_test_dd = true;
    flag_test_qd = true;
  }

  if (flag_test_double) {
    TestSuite<double> test;

    cout << endl;
    cout << "Timing double" << endl;
    cout << "-------------" << endl;
    test.testall();
  }

  if (flag_test_dd) {
    TestSuite<dd_real> test;

    cout << endl;
    cout << "Timing dd_real" << endl;
    cout << "--------------" << endl;
    test.testall();
  }

  if (flag_test_qd) {
    TestSuite<qd_real> test;

    cout << endl;
    cout << "Timing qd_real" << endl;
    cout << "--------------" << endl;
    test.testall();
  }
  
  fpu_fix_end(&old_cw);
  return 0;
}

