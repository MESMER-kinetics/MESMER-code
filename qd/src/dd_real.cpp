/*
 * src/dd_real.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Contains implementation of non-inlined functions of double-double
 * package.  Inlined functions are found in dd_inline.h (in include directory).
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "config.h"
#include <qd/dd_real.h>
#include "util.h"

#include <qd/bits.h>

#ifndef QD_INLINE
#include <qd/dd_inline.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;
using std::ios_base;
using std::string;
using std::setw;

/* This routine is called whenever a fatal error occurs. */
void dd_real::abort(const char *msg) { 
  if (msg) { cerr << "ERROR " << msg << endl; }
}

/* Computes the square root of the double-double number dd.
   NOTE: dd must be a non-negative number.                   */
QD_API dd_real sqrt(const dd_real &a) {
  /* Strategy:  Use Karp's trick:  if x is an approximation
     to sqrt(a), then

        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision.
  */

  if (a.is_zero())
    return 0.0;

  if (a.is_negative()) {
    dd_real::abort("(dd_real::sqrt): Negative argument.");
    return dd_real::_nan;
  }

  double x = 1.0 / std::sqrt(a.hi);
  double ax = a.hi * x;
  return dd_real::add(ax, (a - dd_real::sqr(ax)).hi * (x * 0.5));
}

/* Computes the square root of a double in double-double precision. 
   NOTE: d must not be negative.                                   */
dd_real dd_real::sqrt(double d) {
  return ::sqrt(dd_real(d));
}

/* Computes the n-th root of the double-double number a.
   NOTE: n must be a positive integer.  
   NOTE: If n is even, then a must not be negative.       */
dd_real nroot(const dd_real &a, int n) {
  /* Strategy:  Use Newton iteration for the function

          f(x) = x^(-n) - a

     to find its root a^{-1/n}.  The iteration is thus

          x' = x + x * (1 - a * x^n) / n

     which converges quadratically.  We can then find 
    a^{1/n} by taking the reciprocal.
  */

  if (n <= 0) {
    dd_real::abort("(dd_real::nroot): N must be positive.");
    return dd_real::_nan;
  }

  if (n%2 == 0 && a.is_negative()) {
    dd_real::abort("(dd_real::nroot): Negative argument.");
    return dd_real::_nan;
  }

  if (n == 1) {
    return a;
  } 
  if (n == 2) {
    return sqrt(a);
  }

  if (a.is_zero())
    return 0.0;

  /* Note  a^{-1/n} = exp(-log(a)/n) */
  dd_real r = abs(a);
  dd_real x = std::exp(-std::log(r.hi) / n);

  /* Perform Newton's iteration. */
  x += x * (1.0 - r * npwr(x, n)) / static_cast<double>(n);
  if (a.hi < 0.0)
    x = -x;
  return 1.0/x;
}

/* Computes the n-th power of a double-double number. 
   NOTE:  0^0 causes an error.                         */
dd_real npwr(const dd_real &a, int n) {
  
  if (n == 0) {
    if (a.is_zero()) {
      dd_real::abort("(dd_real::npwr): Invalid argument.");
      return dd_real::_nan;
    }
    return 1.0;
  }

  dd_real r = a;
  dd_real s = 1.0;
  int N = std::abs(n);

  if (N > 1) {
    /* Use binary exponentiation */
    while (N > 0) {
      if (N % 2 == 1) {
        s *= r;
      }
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }
  } else {
    s = r;
  }

  /* Compute the reciprocal if n is negative. */
  if (n < 0)
    return (1.0 / s);
  
  return s;
}

dd_real pow(const dd_real &a, int n) {
  return npwr(a, n);
}

dd_real pow(const dd_real &a, const dd_real &b) {
  return exp(b * log(a));
}

/* Exponential.  Computes exp(x) in double-double precision. */
dd_real exp(const dd_real &a) {
  /* Strategy:  We first reduce the size of x by noting that
     
          exp(kr + m) = exp(m) * exp(r)^k

     Thus by choosing m to be a multiple of log(2) closest
     to x, we can make |kr| <= log(2) / 2 = 0.3466.  Now
     we can set k = 64, so that |r| <= 0.000542.  Then

          exp(x) = exp(kr + s log 2) = (2^s) * [exp(r)]^64

     Then exp(r) is evaluated using the familiar Taylor series.
     Reducing the argument substantially speeds up the convergence.
  */  

  const int k = 64;

  if (a.hi <= -709.0)
    return 0.0;

  if (a.hi >=  709.0)
    return dd_real::_inf;

  if (a.is_zero()) {
    return 1.0;
  }

  if (a.is_one()) {
    return dd_real::_e;
  }

  int z = to_int(nint(a / dd_real::_log2));
  dd_real r = (a - dd_real::_log2 * static_cast<double>(z)) / static_cast<double>(k);
  dd_real s, t, f, p;
  double m;

  s = 1.0 + r;
  p = sqr(r);
  m = 2.0;
  f = 2.0;
  t = p / f;
  do {
    s += t;
    p *= r;
    m += 1.0;
    f *= m;
    t = p / f;
  } while (std::abs(to_double(t)) > 1.0e-35);

  s += t;
  r = pow(s, k);
  r = mul_pwr2(r, std::ldexp(1.0, z));

  return r;
}

/* Logarithm.  Computes log(x) in double-double precision.
   This is a natural logarithm (i.e., base e).            */
dd_real log(const dd_real &a) {
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x) 
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.
           
     Only one iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

  if (a.is_one()) {
    return 0.0;
  }

  if (a.hi <= 0.0) {
    dd_real::abort("(dd_real::log): Non-positive argument.");
    return dd_real::_nan;
  }

  dd_real x = std::log(a.hi);   /* Initial approximation */

  x = x + a * exp(-x) - 1.0;
  return x;
}

dd_real log10(const dd_real &a) {
  return log(a) / dd_real::_log10;
}

/* Computes sin(a) and cos(a) using Taylor series.
   Assumes |a| <= pi/32.                           */
static void sincos_taylor(const dd_real &a, 
                          dd_real &sin_a, dd_real &cos_a) {
  const double thresh = 1.0e-35 * std::abs(to_double(a));
  dd_real t;  /* Term being added. */
  dd_real s;  /* Current partial sum. */
  dd_real x;  /* = -sqr(a) */
  double m;

  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  x = -sqr(a);
  s = a;
  t = a;
  m = 1.0;
  do {
    m += 2.0;
    t *= x;
    t /= (m*(m-1.0));
    s += t;
  } while (std::abs(to_double(t)) > thresh);

  sin_a = s;
  cos_a = sqrt(1.0 - sqr(s));
}

static dd_real sin_taylor(const dd_real &a) {
  const double thresh = 1.0e-35 * std::abs(to_double(a));
  dd_real t;  /* Term being added. */
  dd_real s;  /* Current partial sum. */
  dd_real x;  /* = -sqr(a) */
  double m;

  if (a.is_zero()) {
    return 0.0;
  }

  x = -sqr(a);
  s = a;
  t = a;
  m = 1.0;
  do {
    m += 2.0;
    t *= x;
    t /= (m*(m-1.0));
    s += t;
  } while (std::abs(to_double(t)) > thresh);

  return s;
}

static dd_real cos_taylor(const dd_real &a) {
  const double thresh = 1.0e-35 * std::abs(to_double(a));
  dd_real t;  /* Term being added. */
  dd_real s;  /* Current partial sum. */
  dd_real x;  /* = -sqr(a) */
  double m;

  if (a.is_zero()) {
    return 1.0;
  }

  x = -sqr(a);
  t = 0.5 * x;
  s = 1.0 + t;
  m = 2.0;
  do {
    m += 2.0;
    t *= x;
    t /= (m*(m-1.0));
    s += t;
  } while (std::abs(to_double(t)) > thresh);

  return s;
}

dd_real sin(const dd_real &a) {  

  /* Strategy.  To compute sin(x), we choose integers a, b so that

       x = s + a * (pi/2) + b * (pi/16)

     and |s| <= pi/32.  Using the fact that 

       sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))

     we can compute sin(x) from sin(s), cos(s).  This greatly 
     increases the convergence of the sine Taylor series. */

  if (a.is_zero()) {
    return 0.0;
  }

  /* First reduce modulo 2*pi so that |r| <= pi. */
  dd_real r = drem(a, dd_real::_2pi);

  /* Now reduce by modulo pi/2 and then by pi/16 so that
     we obtain numbers a, b, and t. */
  dd_real t;
  dd_real sin_t, cos_t;
  dd_real s, c;
  int j = to_int(divrem(r, dd_real::_pi2, t));
  int abs_j = std::abs(j);
  int k = to_int(divrem(t, dd_real::_pi16, t));
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    dd_real::abort("(dd_real::sin): Cannot reduce modulo pi/2.");
    return dd_real::_nan;
  }

  if (abs_k > 4) {
    dd_real::abort("(dd_real::sin): Cannot reduce modulo pi/16.");
    return dd_real::_nan;
  }

  if (abs_j == 0) {
    if (k == 0) {
      r = sin_taylor(t);
    } else if (k > 0) {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = u * sin_t + v * cos_t;
    } else {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = u * sin_t - v * cos_t;
    }
  } else if (j == 1) {
    if (k == 0) {
      r = cos_taylor(t);
    } else if (k > 0) {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = u * cos_t - v * sin_t;
    } else {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k == 0) {
      r = -cos_taylor(t);
    } else if (k > 0) {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = v * sin_t - u * cos_t;
    } else if (k < 0) {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = -u * cos_t - v * sin_t;
    }
  } else {
    if (k == 0) {
      r = -sin_taylor(t);
    } else if (k > 0) {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = -u * sin_t - v * cos_t;
    } else {
      dd_real u = dd_real::cos_table[abs_k-1];
      dd_real v = dd_real::sin_table[abs_k-1];
      dd_real sin_t, cos_t;
      sincos_taylor(t, sin_t, cos_t);
      r = v * cos_t - u * sin_t;
    }
  }

  return r;
}

dd_real cos(const dd_real &a) {

  if (a.is_zero()) {
    return 1.0;
  }

  /* First reduce modulo 2*pi so that |r| <= pi. */
  dd_real r = drem(a, dd_real::_2pi);

  /* Now reduce by modulo pi/2 and then by pi/16 so that
     we obtain numbers a, b, and t. */
  dd_real t;
  dd_real sin_t, cos_t;
  dd_real s, c;
  int j = to_int(divrem(r, dd_real::_pi2, t));
  int abs_j = std::abs(j);
  int k = to_int(divrem(t, dd_real::_pi16, t));
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    dd_real::abort("(dd_real::cos): Cannot reduce modulo pi/2.");
    return dd_real::_nan;
  }

  if (abs_k > 4) {
    dd_real::abort("(dd_real::cos): Cannot reduce modulo pi/16.");
    return dd_real::_nan;
  }

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    dd_real u = dd_real::cos_table[abs_k-1];
    dd_real v = dd_real::sin_table[abs_k-1];

    if (k > 0) {
      s = u * sin_t + v * cos_t;
      c = u * cos_t - v * sin_t;
    } else {
      s = u * sin_t - v * cos_t;
      c = u * cos_t + v * sin_t;
    }
  }

  if (abs_j == 0) {
    r = c;
  } else if (j == 1) {
    r = -s;
  } else if (j == -1) {
    r = s;
  } else {
    r = -c;
  }

  return r;
}

void sincos(const dd_real &a, dd_real &sin_a, dd_real &cos_a) {

  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  /* First reduce modulo 2*pi so that |r| <= pi. */
  dd_real r = drem(a, dd_real::_2pi);

  /* Now reduce by modulo pi/2 and then by pi/16 so that
     we obtain numbers a, b, and t. */
  dd_real t;
  dd_real sin_t, cos_t;
  dd_real s, c;
  int j = to_int(divrem(r, dd_real::_pi2, t));
  int abs_j = std::abs(j);
  int k = to_int(divrem(t, dd_real::_pi16, t));
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    dd_real::abort("(dd_real::sincos): Cannot reduce modulo pi/2.");
    cos_a = sin_a = dd_real::_nan;
    return;
  }

  if (abs_k > 4) {
    dd_real::abort("(dd_real::sincos): Cannot reduce modulo pi/16.");
    cos_a = sin_a = dd_real::_nan;
    return;
  }

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    dd_real u = dd_real::cos_table[abs_k-1];
    dd_real v = dd_real::sin_table[abs_k-1];

    if (k > 0) {
      s = u * sin_t + v * cos_t;
      c = u * cos_t - v * sin_t;
    } else {
      s = u * sin_t - v * cos_t;
      c = u * cos_t + v * sin_t;
    }
  }

  if (abs_j == 0) {
    sin_a = s;
    cos_a = c;
  } else if (j == 1) {
    sin_a = c;
    cos_a = -s;
  } else if (j == -1) {
    sin_a = -c;
    cos_a = s;
  } else {
    sin_a = -s;
    cos_a = -c;
  }
  
}

dd_real atan(const dd_real &a) {
  return atan2(a, dd_real(1.0));
}

dd_real atan2(const dd_real &y, const dd_real &x) {
  /* Strategy: Instead of using Taylor series to compute 
     arctan, we instead use Newton's iteration to solve
     the equation

        sin(z) = y/r    or    cos(z) = x/r

     where r = sqrt(x^2 + y^2).
     The iteration is given by

        z' = z + (y - sin(z)) / cos(z)          (for equation 1)
        z' = z - (x - cos(z)) / sin(z)          (for equation 2)

     Here, x and y are normalized so that x^2 + y^2 = 1.
     If |x| > |y|, then first iteration is used since the 
     denominator is larger.  Otherwise, the second is used.
  */

  if (x.is_zero()) {
    
    if (y.is_zero()) {
      /* Both x and y is zero. */
      dd_real::abort("(dd_real::atan2): Both arguments zero.");
      return dd_real::_nan;
    }

    return (y.is_positive()) ? dd_real::_pi2 : -dd_real::_pi2;
  } else if (y.is_zero()) {
    return (x.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }

  if (x == y) {
    return (y.is_positive()) ? dd_real::_pi4 : -dd_real::_3pi4;
  }

  if (x == -y) {
    return (y.is_positive()) ? dd_real::_3pi4 : -dd_real::_pi4;
  }

  dd_real r = sqrt(sqr(x) + sqr(y));
  dd_real xx = x / r;
  dd_real yy = y / r;

  /* Compute double precision approximation to atan. */
  dd_real z = std::atan2(to_double(y), to_double(x));
  dd_real sin_z, cos_z;

  if (xx > yy) {
    /* Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)  */
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
  } else {
    /* Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)  */
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
  }

  return z;
}

dd_real tan(const dd_real &a) {
  dd_real s, c;
  sincos(a, s, c);
  return s/c;
}

dd_real asin(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    dd_real::abort("(dd_real::asin): Argument out of domain.");
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real::_pi2 : -dd_real::_pi2;
  }

  return atan2(a, sqrt(1.0 - sqr(a)));
}

dd_real acos(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    dd_real::abort("(dd_real::acos): Argument out of domain.");
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }

  return atan2(sqrt(1.0 - sqr(a)), a);
}
 
dd_real sinh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (abs(a) > 0.05) {
    dd_real ea = exp(a);
    return mul_pwr2(ea - inv(ea), 0.5);
  }

  /* since a is small, using the above formula gives
     a lot of cancellation.  So use Taylor series.   */
  dd_real s = a;
  dd_real t = a;
  dd_real r = sqr(t);
  double m = 1.0;
  double thresh = std::abs((to_double(a)) * dd_real::_eps);

  do {
    m += 2.0;
    t *= r;
    t /= (m-1) * m;

    s += t;
  } while (abs(t) > thresh);

  return s;

}

dd_real cosh(const dd_real &a) {
  if (a.is_zero()) {
    return 1.0;
  }

  dd_real ea = exp(a);
  return mul_pwr2(ea + inv(ea), 0.5);
}

dd_real tanh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  dd_real ea = exp(a);
  dd_real inv_ea = inv(ea);
  return (ea - inv_ea) / (ea + inv_ea);
}

void sincosh(const dd_real &a, dd_real &sinh_a, dd_real &cosh_a) {
  sinh_a = sinh(a);
  cosh_a = cosh(a);
}

dd_real asinh(const dd_real &a) {
  return log(a + sqrt(sqr(a) + 1.0));
}

dd_real acosh(const dd_real &a) {
  if (a < 1.0) {
    dd_real::abort("(dd_real::acosh): Argument out of domain.");
    return dd_real::_nan;
  }

  return log(a + sqrt(sqr(a) - 1.0));
}

dd_real atanh(const dd_real &a) {
  if (abs(a) >= 1.0) {
    dd_real::abort("(dd_real::atanh): Argument out of domain.");
    return dd_real::_nan;
  }

  return mul_pwr2(log((1.0 + a) / (1.0 - a)), 0.5);
}

QD_API dd_real ddrand() {
  static const double m_const = 4.6566128730773926e-10;  /* = 2^{-31} */
  double m = m_const;
  dd_real r = 0.0;
  double d;

  /* Strategy:  Generate 31 bits at a time, using lrand48 
     random number generator.  Shift the bits, and reapeat
     4 times. */

  for (int i = 0; i < 4; i++, m *= m_const) {
//    d = lrand48() * m;
    d = std::rand() * m;
    r += d;
  }

  return r;
}

/* polyeval(c, n, x)
   Evaluates the given n-th degree polynomial at x.
   The polynomial is given by the array of (n+1) coefficients. */
dd_real polyeval(const dd_real *c, int n, const dd_real &x) {
  /* Just use Horner's method of polynomial evaluation. */
  dd_real r = c[n];
  
  for (int i = n-1; i >= 0; i--) {
    r *= x;
    r += c[i];
  }

  return r;
}

/* polyroot(c, n, x0)
   Given an n-th degree polynomial, finds a root close to 
   the given guess x0.  Note that this uses simple Newton
   iteration scheme, and does not work for multiple roots.  */
QD_API dd_real polyroot(const dd_real *c, int n, 
    const dd_real &x0, int max_iter, double thresh) {
  dd_real x = x0;
  dd_real f;
  dd_real *d = new dd_real[n];
  bool conv = false;
  int i;
  double max_c = std::abs(to_double(c[0]));
  double v;

  if (thresh == 0.0) thresh = dd_real::_eps;

  /* Compute the coefficients of the derivatives. */
  for (i = 1; i <= n; i++) {
    v = std::abs(to_double(c[i]));
    if (v > max_c) max_c = v;
    d[i-1] = c[i] * static_cast<double>(i);
  }
  thresh *= max_c;

  /* Newton iteration. */
  for (i = 0; i < max_iter; i++) {
    f = polyeval(c, n, x);

    if (abs(f) < thresh) {
      conv = true;
      break;
    }
    x -= (f / polyeval(d, n-1, x));
  }
  delete [] d;

  if (!conv) {
    dd_real::abort("(dd_real::polyroot): Failed to converge.");
    return dd_real::_nan;
  }

  return x;
}


/* Constructor.  Reads a double-double number from the string s
   and constructs a double-double number.                         */
dd_real::dd_real(const char *s) {
  if (dd_real::read(s, *this)) {
    dd_real::abort("(dd_real::dd_real): INPUT ERROR.");
    *this = dd_real::_nan;
  }
}

dd_real &dd_real::operator=(const char *s) {
  if (dd_real::read(s, *this)) {
    dd_real::abort("(dd_real::operator=): INPUT ERROR.");
    *this = dd_real::_nan;
  }
  return *this;
}

/* Outputs the double-double number dd. */
ostream &operator<<(ostream &os, const dd_real &dd) {
  bool showpos = (os.flags() & ios_base::showpos) != 0;
  bool uppercase =  (os.flags() & ios_base::uppercase) != 0;
  ios_base::fmtflags float_field = os.flags() & ios_base::floatfield;
  ios_base::fmtflags adjust_field = os.flags() & ios_base::adjustfield;

  return os << dd.write(os.precision(), os.width(), float_field, 
      adjust_field, showpos, uppercase, os.fill());
}

/* Reads in the double-double number a. */
istream &operator>>(istream &s, dd_real &a) {
  char str[255];
  s >> str;
  a = dd_real(str);
  return s;
}

void dd_real::to_digits(char *s, int &expn, int precision) const {
  int D = precision + 1;  /* number of digits to compute */

  dd_real r = abs(*this);
  int e;  /* exponent */
  int i, d;

  if (hi == 0.0) {
    /* this == 0.0 */
    expn = 0;
    for (i = 0; i < precision; i++) s[i] = '0';
    return;
  }

  /* First determine the (approximate) exponent. */
  e = to_int(std::floor(std::log10(std::abs(hi))));

  if (e < -300) {
    r *= dd_real(10.0) ^ 300;
    r /= dd_real(10.0) ^ (e + 300);
  } else if (e > 300) {
    r = ldexp(r, -53);
    r /= dd_real(10.0) ^ e;
    r = ldexp(r, 53);
  } else {
    r /= dd_real(10.0) ^ e;
  }

  /* Fix exponent if we are off by one */
  if (r >= 10.0) {
    r /= 10.0;
    e++;
  } else if (r < 1.0) {
    r *= 10.0;
    e--;
  }

  if (r >= 10.0 || r < 1.0) {
    dd_real::abort("(dd_real::to_digits): can't compute exponent.");
    return;
  }

  /* Extract the digits */
  for (i = 0; i < D; i++) {
    d = static_cast<int>(r.hi);
    r -= d;
    r *= 10.0;

    s[i] = static_cast<char>(d + '0');
  }

  /* Fix negative digits. */
  for (i = D-1; i > 0; i--) {
    if (s[i] < '0') {
      s[i-1]--;
      s[i] += 10;
    }
  }

  if (s[0] <= '0') {
    dd_real::abort("(dd_real::to_digits): non-positive leading digit.");
    return;
  }

  /* Round, handle carry */
  if (s[D-1] >= '5') {
    s[D-2]++;

    i = D-2;
    while (i > 0 && s[i] > '9') {
      s[i] -= 10;
      s[--i]++;
    }
  }

  /* If first digit is 10, shift everything. */
  if (s[0] > '9') { 
    e++; 
    for (i = precision; i >= 2; i--) s[i] = s[i-1]; 
    s[0] = '1';
    s[1] = '0';
  }

  s[precision] = 0;
  expn = e;
}

/* Writes the double-double number into the string s.
   The integer d specifies how many significant digits to write.
   The string s must be able to hold at least (d+8) characters.  
   showpos indicates whether to use the + sign, and uppercase indicates
   whether the E or e is to be used for the exponent. */
void dd_real::write(char *s, int precision, bool showpos, bool uppercase) const {
  char *t = new char[precision + 1];
  int e, i, j;

  if (isnan()) {
    if (uppercase) {
      s[0] = 'N'; s[1] = 'A'; s[2] = 'N';
    } else {
      s[0] = 'n'; s[1] = 'a'; s[2] = 'n';
    }
    s[3] = 0;
    return;
  }

  i = 0;
  if (hi < 0.0)
    s[i++] = '-';
  else if (hi >= 0.0 && showpos)
    s[i++] = '+';

  if (isinf()) {
    if (uppercase) {
      s[i++] = 'i'; s[i++] = 'n'; s[i++] = 'f';
    } else {
      s[i++] = 'I'; s[i++] = 'N'; s[i++] = 'F';
    }
    s[i++] = 0;
    return;
  }

  to_digits(t, e, precision);

  s[i++] = t[0];
  s[i++] = '.';

  for (j = 1; j < precision; j++, i++)
    s[i] = t[j];

  /* Fill in exponent part */
  s[i++] = uppercase ? 'E' : 'e';
  std::sprintf(&s[i], "%d", e);
  delete [] t;
}

string dd_real::write(int precision, int width, 
    ios_base::fmtflags float_field, ios_base::fmtflags adjust_field, 
    bool showpos, bool uppercase, char fill) const {
  string s;
  bool fixed = (float_field & ios_base::fixed) != 0;
  bool sgn = true;
  int i, e = 0;

  if (isnan()) {
    s = uppercase ? "NAN" : "nan";
    sgn = false;
  } else {
    if (*this < 0.0)
      s += '-';
    else if (showpos)
      s += '+';
    else
      sgn = false;

    if (isinf()) {
      s += uppercase ? "INF" : "inf";
    } else if (*this == 0.0) {
      /* Zero case */
      s += '0';
      if (precision > 0) {
        s += '.';
        s.append(precision, '0');
      }
    } else {
      /* Non-zero case */
      int off = (fixed ? (1 + to_int(floor(log10(abs(*this))))) : 1);
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

    if (!fixed && !isinf()) {
      /* Fill in exponent part */
      s += uppercase ? 'E' : 'e';
      append_expn(s, e);
    }
  }

  /* Fill in the blanks */
  int len = s.length();
  if (len < width) {
    int delta = width - len;
    if (adjust_field & ios_base::internal) {
      if (sgn)
        s.insert(static_cast<string::size_type>(1), delta, fill);
      else
        s.insert(static_cast<string::size_type>(0), delta, fill);
    } else if (adjust_field & ios_base::left) {
      s.append(delta, fill);
    } else {
      s.insert(static_cast<string::size_type>(0), delta, fill);
    }
  }

  return s;
}

/* Reads in a double-double number from the string s. */
int dd_real::read(const char *s, dd_real &a) {
  const char *p = s;
  char ch;
  int sign = 0;
  int point = -1;
  int nd = 0;
  int e = 0;
  bool done = false;
  dd_real r = 0.0;
  int nread;
  
  /* Skip any leading spaces */
  while (*p == ' ')
    p++;

  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') {
      int d = ch - '0';
      r *= 10.0;
      r += static_cast<double>(d);
      nd++;
    } else {

      switch (ch) {

      case '.':
        if (point >= 0)
          return -1;        
        point = nd;
        break;

      case '-':
      case '+':
        if (sign != 0 || nd > 0)
          return -1;
        sign = (ch == '-') ? -1 : 1;
        break;

      case 'E':
      case 'e':
        nread = std::sscanf(p+1, "%d", &e);
        done = true;
        if (nread != 1)
          return -1;
        break;

      default:
        return -1;
      }
    }
    
    p++;
  }

  if (point >= 0) {
    e -= (nd - point);
  }

  if (e != 0) {
    r *= (dd_real(10.0) ^ e);
  }

  a = (sign == -1) ? -r : r;
  return 0;
}

/* Debugging routines */
void dd_real::dump(const string &name, std::ostream &os) const {
  std::ios_base::fmtflags old_flags = os.flags();
  std::streamsize old_prec = os.precision(19);
  os << std::scientific;

  if (name.length() > 0) os << name << " = ";
  os << "[ " << setw(27) << hi << ", " << setw(27) << lo << " ]" << endl;

  os.precision(old_prec);
  os.flags(old_flags);
}

void dd_real::dump_bits(const string &name, std::ostream &os) const {
  string::size_type len = name.length();
  if (len > 0) {
    os << name << " = ";
    len +=3;
  }
  os << "[ ";
  len += 2;
  print_double_info(os, hi);
  os << endl;
  for (string::size_type i = 0; i < len; i++) os << ' ';
  print_double_info(os, lo);
  os << " ]" << endl;
}

dd_real dd_real::debug_rand() { 

  if (std::rand() % 2 == 0)
    return ddrand();

  int expn = 0;
  dd_real a = 0.0;
  double d;
  for (int i = 0; i < 2; i++) {
    d = std::ldexp(static_cast<double>(std::rand()) / RAND_MAX, -expn);
    a += d;
    expn = expn + 54 + std::rand() % 200;
  }
  return a;
}
