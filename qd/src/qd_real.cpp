/*
 * src/qd_real.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Contains implementation of non-inlined functions of quad-double
 * package.  Inlined functions are found in qd_inline.h (in include directory).
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "config.h"
#include <qd/qd_real.h>
#include "util.h"

#include <qd/bits.h>

#ifndef QD_INLINE
#include <qd/qd_inline.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::ios_base;
using std::string;
using std::setw;

using namespace qd;

void qd_real::abort(const char *msg) {
  if (msg) { cerr << "ERROR " << msg << endl; }
}

/********** Multiplications **********/

qd_real nint(const qd_real &a) {
  double x0, x1, x2, x3;

  x0 = nint(a[0]);
  x1 = x2 = x3 = 0.0;

  if (x0 == a[0]) {
    /* First double is already an integer. */
    x1 = nint(a[1]);

    if (x1 == a[1]) {
      /* Second double is already an integer. */
      x2 = nint(a[2]);
      
      if (x2 == a[2]) {
        /* Third double is already an integer. */
        x3 = nint(a[3]);
      } else {
        if (std::abs(x2 - a[2]) == 0.5 && a[3] < 0.0) {
          x2 -= 1.0;
        }
      }

    } else {
      if (std::abs(x1 - a[1]) == 0.5 && a[2] < 0.0) {
          x1 -= 1.0;
      }
    }

  } else {
    /* First double is not an integer. */
      if (std::abs(x0 - a[0]) == 0.5 && a[1] < 0.0) {
          x0 -= 1.0;
      }
  }
  
  renorm(x0, x1, x2, x3);
  return qd_real(x0, x1, x2, x3);
}

qd_real floor(const qd_real &a) {
  double x0, x1, x2, x3;
  x1 = x2 = x3 = 0.0;
  x0 = std::floor(a[0]);

  if (x0 == a[0]) {
    x1 = std::floor(a[1]);
    
    if (x1 == a[1]) {
      x2 = std::floor(a[2]);

      if (x2 == a[2]) {
        x3 = std::floor(a[3]);
      }
    }

    renorm(x0, x1, x2, x3);
    return qd_real(x0, x1, x2, x3);
  }

  return qd_real(x0, x1, x2, x3);
}

qd_real ceil(const qd_real &a) {
  double x0, x1, x2, x3;
  x1 = x2 = x3 = 0.0;
  x0 = std::ceil(a[0]);

  if (x0 == a[0]) {
    x1 = std::ceil(a[1]);
    
    if (x1 == a[1]) {
      x2 = std::ceil(a[2]);

      if (x2 == a[2]) {
        x3 = std::ceil(a[3]);
      }
    }

    renorm(x0, x1, x2, x3);
    return qd_real(x0, x1, x2, x3);
  }

  return qd_real(x0, x1, x2, x3);
}



/********** Divisions **********/
/* quad-double / double */
qd_real operator/(const qd_real &a, double b) {
  /* Strategy:  compute approximate quotient using high order
     doubles, and then correct it 3 times using the remainder.
     (Analogous to long division.)                             */
  double t0, t1;
  double q0, q1, q2, q3;
  qd_real r;

  q0 = a[0] / b;  /* approximate quotient */

  /* Compute the remainder  a - q0 * b */
  t0 = two_prod(q0, b, t1);
  r = a - dd_real(t0, t1);

  /* Compute the first correction */
  q1 = r[0] / b;
  t0 = two_prod(q1, b, t1);
  r -= dd_real(t0, t1);

  /* Second correction to the quotient. */
  q2 = r[0] / b;
  t0 = two_prod(q2, b, t1);
  r -= dd_real(t0, t1);

  /* Final correction to the quotient. */
  q3 = r[0] / b;

  renorm(q0, q1, q2, q3);
  return qd_real(q0, q1, q2, q3);
}

qd_real::qd_real(const char *s) {
  if (qd_real::read(s, *this)) {
    qd_real::abort("(qd_real::qd_real): INPUT ERROR.");
    *this = qd_real::_nan;
  }
}

qd_real &qd_real::operator=(const char *s) {
  if (qd_real::read(s, *this)) {
    qd_real::abort("(qd_real::operator=): INPUT ERROR.");
    *this = qd_real::_nan;
  }
  return *this;
}

istream &operator>>(istream &s, qd_real &qd) {
  char str[255];
  s >> str;
  qd = qd_real(str);
  return s;
}

ostream &operator<<(ostream &os, const qd_real &qd) {
  bool showpos = (os.flags() & ios_base::showpos) != 0;
  bool uppercase = (os.flags() & ios_base::uppercase) != 0;
  ios_base::fmtflags float_field = os.flags() & ios_base::floatfield;
  ios_base::fmtflags adjust_field = os.flags() & ios_base::adjustfield;
  return os << qd.write(os.precision(), os.width(), float_field, 
      adjust_field, showpos, uppercase, os.fill());
}

/* Read a quad-double from s. */
int qd_real::read(const char *s, qd_real &qd) {
  const char *p = s;
  char ch;
  int sign = 0;
  int point = -1;  /* location of decimal point */
  int nd = 0;      /* number of digits read */
  int e = 0;       /* exponent. */
  bool done = false;
  qd_real r = 0.0;  /* number being read */

  /* Skip any leading spaces */
  while (*p == ' ') p++;

  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') {
      /* It's a digit */
      int d = ch - '0';
      r *= 10.0;
      r += static_cast<double>(d);
      nd++;
    } else {
      /* Non-digit */
      switch (ch) {
      case '.':
        if (point >= 0)
          return -1;   /* we've already encountered a decimal point. */
        point = nd;
        break;
      case '-':
      case '+':
        if (sign != 0 || nd > 0)
          return -1;  /* we've already encountered a sign, or if its
                            not at first position. */
        sign = (ch == '-') ? -1 : 1;
        break;
      case 'E':
      case 'e':
        int nread;
        nread = std::sscanf(p+1, "%d", &e);
        done = true;
        if (nread != 1)
          return -1;  /* read of exponent failed. */
        break;
      case ' ':
        done = true;
        break;
      default:
        return -1;
      
      }
    }

    p++;
  }



  /* Adjust exponent to account for decimal point */
  if (point >= 0) {
    e -= (nd - point);
  }

  /* Multiply the the exponent */
  if (e != 0) {
    r *= (qd_real(10.0) ^ e);
  }

  qd = (sign < 0) ? -r : r;
  return 0;
}

void qd_real::to_digits(char *s, int &expn, int precision) const {
  int D = precision + 1;  /* number of digits to compute */

  qd_real r = abs(*this);
  int e;  /* exponent */
  int i, d;

  if (x[0] == 0.0) {
    /* this == 0.0 */
    expn = 0;
    for (i = 0; i < precision; i++) s[i] = '0';
    return;
  }

  /* First determine the (approximate) exponent. */
  e = static_cast<int>(std::floor(std::log10(std::abs(x[0]))));

  if (e < -300) {
    r *= qd_real(10.0) ^ 300;
    r /= qd_real(10.0) ^ (e + 300);
  } else if (e > 300) {
    r = ldexp(r, -53);
    r /= qd_real(10.0) ^ e;
    r = ldexp(r, 53);
  } else {
    r /= qd_real(10.0) ^ e;
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
    qd_real::abort("(qd_real::to_digits): can't compute exponent.");
    return;
  }

  /* Extract the digits */
  for (i = 0; i < D; i++) {
    d = static_cast<int>(r[0]);
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
    qd_real::abort("(qd_real::to_digits): non-positive leading digit.");
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

/* Writes the quad-double number into the string s.
   The integer d specifies how many significant digits to write.
   The string s must be able to hold at least (d+8) characters.  
   showpos indicates whether to use the + sign, and uppercase indicates
   whether the E or e is to be used for the exponent. */
void qd_real::write(char *s, int precision, bool showpos, bool uppercase) const {
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
  if (x[0] < 0.0)
    s[i++] = '-';
  else if (x[0] >= 0.0 && showpos)
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

string qd_real::write(int precision, int width, 
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

/* Computes  qd^n, where n is an integer. */
qd_real pow(const qd_real &a, int n) {
  if (n == 0)
    return 1.0;

  qd_real r = a;   /* odd-case multiplier */
  qd_real s = 1.0;  /* current answer */
  int N = std::abs(n);

  if (N > 1) {

    /* Use binary exponentiation. */
    while (N > 0) {
      if (N % 2 == 1) {
        /* If odd, multiply by r. Note eventually N = 1, so this
         eventually executes. */
        s *= r;
      }
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }

  } else {
    s = r;
  }

  if (n < 0)
    return (1.0 / s);

  return s;
}

qd_real pow(const qd_real &a, const qd_real &b) {
  return exp(b * log(a));
}

qd_real npwr(const qd_real &a, int n) {
  return pow(a, n);
}

/* Debugging routines */
void qd_real::dump_bits(const string &name, std::ostream &os) const {
  string::size_type len = name.length();
  if (len > 0) {
    os << name << " = ";
    len += 3;
  }
  os << "[ ";
  len += 2;
  for (int j = 0; j < 4; j++) {
    if (j > 0) for (string::size_type i = 0; i < len; i++) os << ' ';
    print_double_info(os, x[j]);
    if (j < 3)
      os << endl;
    else
      os << " ]" << endl;
  }
}

void qd_real::dump(const string &name, std::ostream &os) const {
  std::ios_base::fmtflags old_flags = os.flags();
  std::streamsize old_prec = os.precision(19);
  os << std::scientific;

  string::size_type len = name.length();
  if (len > 0) {
    os << name << " = ";
    len += 3;
  }
  os << "[ ";
  len += 2;
  os << setw(27) << x[0] << ", " << setw(26) << x[1] << "," << endl;
  for (string::size_type i = 0; i < len; i++) os << ' ';
  os << setw(27) << x[2] << ", " << setw(26) << x[3] << "  ]" << endl;

  os.precision(old_prec);
  os.flags(old_flags);
}

/* Divisions */
/* quad-double / double-double */
qd_real qd_real::sloppy_div(const qd_real &a, const dd_real &b) {
  double q0, q1, q2, q3;
  qd_real r;

  q0 = a[0] / b._hi();
  r = a - q0 * b;

  q1 = r[0] / b._hi();
  r -= (q1 * b);

  q2 = r[0] / b._hi();
  r -= (q2 * b);

  q3 = r[0] / b._hi();

  ::renorm(q0, q1, q2, q3);
  return qd_real(q0, q1, q2, q3);
}

qd_real qd_real::accurate_div(const qd_real &a, const dd_real &b) {
  double q0, q1, q2, q3;
  qd_real r;

  q0 = a[0] / b._hi();
  r = a - q0 * b;

  q1 = r[0] / b._hi();
  r -= (q1 * b);

  q2 = r[0] / b._hi();
  r -= (q2 * b);

  q3 = r[0] / b._hi();

  r -= (q3 * b);
  double q4 = r[0] / b._hi();

  ::renorm(q0, q1, q2, q3, q4);
  return qd_real(q0, q1, q2, q3);
}

/* quad-double / quad-double */
qd_real qd_real::sloppy_div(const qd_real &a, const qd_real &b) {
  double q0, q1, q2, q3;

  qd_real r;

  q0 = a[0] / b[0];
  r = a - (b * q0);

  q1 = r[0] / b[0];
  r -= (b * q1);

  q2 = r[0] / b[0];
  r -= (b * q2);

  q3 = r[0] / b[0];

  ::renorm(q0, q1, q2, q3);

  return qd_real(q0, q1, q2, q3);
}

qd_real qd_real::accurate_div(const qd_real &a, const qd_real &b) {
  double q0, q1, q2, q3;

  qd_real r;

  q0 = a[0] / b[0];
  r = a - (b * q0);

  q1 = r[0] / b[0];
  r -= (b * q1);

  q2 = r[0] / b[0];
  r -= (b * q2);

  q3 = r[0] / b[0];

  r -= (b * q3);
  double q4 = r[0] / b[0];

  ::renorm(q0, q1, q2, q3, q4);

  return qd_real(q0, q1, q2, q3);
}

QD_API qd_real sqrt(const qd_real &a) {
  /* Strategy:  

     Perform the following Newton iteration:

       x' = x + (1 - a * x^2) * x / 2;
       
     which converges to 1/sqrt(a), starting with the
     double precision approximation to 1/sqrt(a).
     Since Newton's iteration more or less doubles the
     number of correct digits, we only need to perform it 
     twice.
  */

  if (a.is_zero())
    return 0.0;

  if (a.is_negative()) {
    qd_real::abort("(qd_real::sqrt): Negative argument.");
    return qd_real::_nan;
  }

  qd_real r = (1.0 / std::sqrt(a[0]));
  qd_real h = a * 0.5;

  r += ((0.5 - h * sqr(r)) * r);
  r += ((0.5 - h * sqr(r)) * r);
  r += ((0.5 - h * sqr(r)) * r);

  r *= a;
  return r;
}


/* Computes the n-th root of a */
qd_real nroot(const qd_real &a, int n) {
  /* Strategy:  Use Newton's iteration to solve
     
        1/(x^n) - a = 0

     Newton iteration becomes

        x' = x + x * (1 - a * x^n) / n

     Since Newton's iteration converges quadratically, 
     we only need to perform it twice.

   */

  if (a == 0.0) {
    return qd_real(0.0);
  }

  qd_real r = std::pow(a[0], -1.0/n);

  double dbl_n = static_cast<double>(n);
  r += r * (1.0 - a * (r ^ n)) / dbl_n;
  r += r * (1.0 - a * (r ^ n)) / dbl_n;
  r += r * (1.0 - a * (r ^ n)) / dbl_n;

  return 1.0 / r;
}

qd_real exp(const qd_real &a) {
  /* Strategy:  We first reduce the size of x by noting that
     
          exp(kr + m) = exp(m) * exp(r)^k

     Thus by choosing m to be a multiple of log(2) closest
     to x, we can make |kr| <= log(2) / 2 = 0.3466.  Now
     we can set k = 256, so that |r| <= 0.00136.  Then

          exp(x) = exp(kr + s log 2) = (2^s) * [exp(r)]^256

     Then exp(r) is evaluated using the familiar Taylor series.
     Reducing the argument substantially speeds up the convergence.
  */  

  const int k = 256;

  if (a[0] <= -709.0)
    return 0.0;

  if (a[0] >=  709.0)
    return qd_real::_inf;

  if (a.is_zero()) {
    return 1.0;
  }

  if (a.is_one()) {
    return qd_real::_e;
  }

  int z = to_int(nint(a / qd_real::_log2));
  qd_real r = (a - qd_real::_log2 * static_cast<double>(z)) / static_cast<double>(k);
  qd_real s, p;
  double m;
  double thresh = qd_real::_eps;

  p = sqr(r) / 2.0;
  s = 1.0 + r + p;
  m = 2.0;
  do {
    m += 1.0;
    p *= r;
    p /= m;
    s += p;
  } while (std::abs(to_double(p)) > thresh);

  r = pow(s, k);
  r = mul_pwr2(r, std::ldexp(1.0, z));

  return r;  
}

/* Logarithm.  Computes log(x) in quad-double precision.
   This is a natural logarithm (i.e., base e).            */
qd_real log(const qd_real &a) {
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x) 
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.
           
     Two iteration is needed, since Newton's iteration 
     approximately doubles the number of digits per iteration. */

  if (a.is_one()) {
    return 0.0;
  }

  if (a[0] <= 0.0) {
    qd_real::abort("(qd_real::log): Non-positive argument.");
    return qd_real::_nan;
  }

  qd_real x = std::log(a[0]);   /* Initial approximation */

  x = x + a * exp(-x) - 1.0;
  x = x + a * exp(-x) - 1.0;
  x = x + a * exp(-x) - 1.0;

  return x;
}

qd_real log10(const qd_real &a) {
  return log(a) / qd_real::_log10;
}

/* Computes sin(a) and cos(a) using Taylor series.
   Assumes |a| <= pi/2048.                           */
static void sincos_taylor(const qd_real &a, 
                          qd_real &sin_a, qd_real &cos_a) {
  const double thresh = qd_real::_eps * std::abs(to_double(a));
  qd_real p;  /* Current power of a. */
  qd_real s;  /* Current partial sum. */
  qd_real x;  /* = -sqr(a) */
  double m;

  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  x = -sqr(a);
  s = a;
  p = a;
  m = 1.0;
  do {
    p *= x;
    m += 2.0;
    p /= (m*(m-1));
    s += p;
  } while (std::abs(to_double(p)) > thresh);

  sin_a = s;
  cos_a = sqrt(1.0 - sqr(s));
}

qd_real sin(const qd_real &a) {

  /* Strategy.  To compute sin(x), we choose integers a, b so that

       x = s + a * (pi/2) + b * (pi/1024)

     and |s| <= pi/2048.  Using a precomputed table of
     sin(k pi / 1024) and cos(k pi / 1024), we can compute
     sin(x) from sin(s) and cos(s).  This greatly increases the
     convergence of the sine Taylor series.                          */

  if (a.is_zero()) {
    return 0.0;
  }

  /* First reduce modulo 2*pi so that |r| <= pi. */
  qd_real r = drem(a, qd_real::_2pi);

  /* Now reduce by modulo pi/2 and then by pi/1024 so that
     we obtain numbers a, b, and t. */
  qd_real t;
  qd_real sin_t, cos_t;
  qd_real s, c;
  int j = to_int(divrem(r, qd_real::_pi2, t));
  int abs_j = std::abs(j);
  int k = to_int(divrem(t, qd_real::_pi1024, t));
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    qd_real::abort("(qd_real::sin): Cannot reduce modulo pi/2.");
    return qd_real::_nan;
  }

  if (abs_k > 256) {
    qd_real::abort("(qd_real::sin): Cannot reduce modulo pi/1024.");
    return qd_real::_nan;
  }

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    qd_real u = qd_real::cos_table[abs_k-1];
    qd_real v = qd_real::sin_table[abs_k-1];

    if (k > 0) {
      s = u * sin_t + v * cos_t;
      c = u * cos_t - v * sin_t;
    } else {
      s = u * sin_t - v * cos_t;
      c = u * cos_t + v * sin_t;
    }
  }

  if (abs_j == 0) {
    r = s;
  } else if (j == 1) {
    r = c;
  } else if (j == -1) {
    r = -c;
  } else {
    r = -s;
  }

  return r;
}

qd_real cos(const qd_real &a) {

  if (a.is_zero()) {
    return 1.0;
  }

  /* First reduce modulo 2*pi so that |r| <= pi. */
  qd_real r = drem(a, qd_real::_2pi);

  /* Now reduce by modulo pi/2 and then by pi/1024 so that
     we obtain numbers a, b, and t. */
  qd_real t;
  qd_real sin_t, cos_t;
  qd_real s, c;
  int j = to_int(divrem(r, qd_real::_pi2, t));
  int abs_j = std::abs(j);
  int k = to_int(divrem(t, qd_real::_pi1024, t));
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    qd_real::abort("(qd_real::cos): Cannot reduce modulo pi/2.");
    return qd_real::_nan;
  }

  if (abs_k > 256) {
    qd_real::abort("(qd_real::cos): Cannot reduce modulo pi/1024.");
    return qd_real::_nan;
  }

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    qd_real u = qd_real::cos_table[abs_k-1];
    qd_real v = qd_real::sin_table[abs_k-1];

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

void sincos(const qd_real &a, qd_real &sin_a, qd_real &cos_a) {

  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  /* First reduce modulo 2*pi so that |r| <= pi. */
  qd_real r = drem(a, qd_real::_2pi);

  /* Now reduce by modulo pi/2 and then by pi/1024 so that
     we obtain numbers a, b, and t. */
  qd_real t;
  qd_real sin_t, cos_t;
  qd_real s, c;
  int j = to_int(divrem(r, qd_real::_pi2, t));
  int abs_j = std::abs(j);
  int k = to_int(divrem(t, qd_real::_pi1024, t));
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    qd_real::abort("(qd_real::sincos): Cannot reduce modulo pi/2.");
    cos_a = sin_a = qd_real::_nan;
    return;
  }

  if (abs_k > 256) {
    qd_real::abort("(qd_real::sincos): Cannot reduce modulo pi/1024.");
    cos_a = sin_a = qd_real::_nan;
    return;
  }

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    qd_real u = qd_real::cos_table[abs_k-1];
    qd_real v = qd_real::sin_table[abs_k-1];

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

qd_real atan(const qd_real &a) {
  return atan2(a, qd_real(1.0));
}

qd_real atan2(const qd_real &y, const qd_real &x) {
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
      qd_real::abort("(qd_real::atan2): Both arguments zero.");
      return qd_real::_nan;
    }

    return (y.is_positive()) ? qd_real::_pi2 : -qd_real::_pi2;
  } else if (y.is_zero()) {
    return (x.is_positive()) ? qd_real(0.0) : qd_real::_pi;
  }

  if (x == y) {
    return (y.is_positive()) ? qd_real::_pi4 : -qd_real::_3pi4;
  }

  if (x == -y) {
    return (y.is_positive()) ? qd_real::_3pi4 : -qd_real::_pi4;
  }

  qd_real r = sqrt(sqr(x) + sqr(y));
  qd_real xx = x / r;
  qd_real yy = y / r;

  /* Compute double precision approximation to atan. */
  qd_real z = std::atan2(to_double(y), to_double(x));
  qd_real sin_z, cos_z;

  if (xx > yy) {
    /* Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)  */
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
  } else {
    /* Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)  */
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
  }

  return z;
}


qd_real drem(const qd_real &a, const qd_real &b) {
  qd_real n = nint(a/b);
  return (a - n * b);
}

qd_real divrem(const qd_real &a, const qd_real &b, qd_real &r) {
  qd_real n = nint(a/b);
  r = a - n * b;
  return n;
}

qd_real tan(const qd_real &a) {
  qd_real s, c;
  sincos(a, s, c);
  return s/c;
}

qd_real asin(const qd_real &a) {
  qd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    qd_real::abort("(qd_real::asin): Argument out of domain.");
    return qd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? qd_real::_pi2 : -qd_real::_pi2;
  }

  return atan2(a, sqrt(1.0 - sqr(a)));
}

qd_real acos(const qd_real &a) {
  qd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    qd_real::abort("(qd_real::acos): Argument out of domain.");
    return qd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? qd_real(0.0) : qd_real::_pi;
  }

  return atan2(sqrt(1.0 - sqr(a)), a);
}
 
qd_real sinh(const qd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (abs(a) > 0.05) {
    qd_real ea = exp(a);
    return mul_pwr2(ea - inv(ea), 0.5);
  }

  /* Since a is small, using the above formula gives
     a lot of cancellation.   So use Taylor series. */
  qd_real s = a;
  qd_real t = a;
  qd_real r = sqr(t);
  double m = 1.0;
  double thresh = std::abs(to_double(a) * qd_real::_eps);

  do {
    m += 2.0;
    t *= r;
    t /= (m-1) * m;

    s += t;
  } while (abs(t) > thresh);

  return s;
}

qd_real cosh(const qd_real &a) {
  if (a.is_zero()) {
    return 1.0;
  }

  qd_real ea = exp(a);
  return mul_pwr2(ea + inv(ea), 0.5);
}

qd_real tanh(const qd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (std::abs(to_double(a)) > 0.05) {
    qd_real ea = exp(a);
    qd_real inv_ea = inv(ea);
    return (ea - inv_ea) / (ea + inv_ea);
  } else {
    qd_real s, c;
    s = sinh(a);
    c = sqrt(1.0 + sqr(s));
    return s / c;
  }
}

void sincosh(const qd_real &a, qd_real &s, qd_real &c) {
  if (std::abs(to_double(a)) <= 0.05) {
    s = sinh(a);
    c = sqrt(1.0 + sqr(s));
  } else {
    qd_real ea = exp(a);
    qd_real inv_ea = inv(ea);
    s = mul_pwr2(ea - inv_ea, 0.5);
    c = mul_pwr2(ea + inv_ea, 0.5);
  }
}

qd_real asinh(const qd_real &a) {
  return log(a + sqrt(sqr(a) + 1.0));
}

qd_real acosh(const qd_real &a) {
  if (a < 1.0) {
    qd_real::abort("(qd_real::acosh): Argument out of domain.");
    return qd_real::_nan;
  }

  return log(a + sqrt(sqr(a) - 1.0));
}

qd_real atanh(const qd_real &a) {
  if (abs(a) >= 1.0) {
    qd_real::abort("(qd_real::atanh): Argument out of domain.");
    return qd_real::_nan;
  }

  return mul_pwr2(log((1.0 + a) / (1.0 - a)), 0.5);
}


QD_API qd_real qdrand() {
  static const double m_const = 4.6566128730773926e-10;  /* = 2^{-31} */
  double m = m_const;
  qd_real r = 0.0;
  double d;

  /* Strategy:  Generate 31 bits at a time, using lrand48 
     random number generator.  Shift the bits, and repeat
     7 times. */

  for (int i = 0; i < 7; i++, m *= m_const) {
    d = std::rand() * m;
    r += d;
  }

  return r;
}


/* polyeval(c, n, x)
   Evaluates the given n-th degree polynomial at x.
   The polynomial is given by the array of (n+1) coefficients. */
qd_real polyeval(const qd_real *c, int n, const qd_real &x) {
  /* Just use Horner's method of polynomial evaluation. */
  qd_real r = c[n];
  
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
QD_API qd_real polyroot(const qd_real *c, int n, 
    const qd_real &x0, int max_iter, double thresh) {
  qd_real x = x0;
  qd_real f;
  qd_real *d = new qd_real[n];
  bool conv = false;
  int i;
  double max_c = std::abs(to_double(c[0]));
  double v;

  if (thresh == 0.0) thresh = qd_real::_eps;

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
    qd_real::abort("(qd_real::polyroot): Failed to converge.");
    return qd_real::_nan;
  }

  return x;
}

qd_real qd_real::debug_rand() {
  if (std::rand() % 2 == 0)
    return qdrand();

  int expn = 0;
  qd_real a = 0.0;
  double d;
  for (int i = 0; i < 4; i++) {
    d = std::ldexp(std::rand() / static_cast<double>(RAND_MAX), -expn);
    a += d;
    expn = expn + 54 + std::rand() % 200;
  }
  return a;
}

