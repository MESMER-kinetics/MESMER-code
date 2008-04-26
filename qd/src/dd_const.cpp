/*
 * src/qd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines constants used in quad-double package.
 */
#include "config.h"
#include <qd/dd_real.h>

const dd_real dd_real::_2pi = dd_real(6.283185307179586232e+00,
                                      2.449293598294706414e-16);
const dd_real dd_real::_pi = dd_real(3.141592653589793116e+00,
                                     1.224646799147353207e-16);
const dd_real dd_real::_pi2 = dd_real(1.570796326794896558e+00,
                                      6.123233995736766036e-17);
const dd_real dd_real::_pi4 = dd_real(7.853981633974482790e-01,
                                      3.061616997868383018e-17);
const dd_real dd_real::_pi8 = dd_real(3.926990816987241395e-01,
                                      1.530808498934191509e-17);
const dd_real dd_real::_pi16 = dd_real(1.963495408493620697e-01,
                                       7.654042494670957545e-18);
const dd_real dd_real::_3pi4 = dd_real(2.356194490192344837e+00,
                                       9.1848509936051484375e-17);
const dd_real dd_real::_e = dd_real(2.718281828459045091e+00,
                                    1.445646891729250158e-16);
const dd_real dd_real::_log2 = dd_real(6.931471805599452862e-01,
                                       2.319046813846299558e-17);
const dd_real dd_real::_log10 = dd_real(2.302585092994045901e+00,
                                        -2.170756223382249351e-16);
const dd_real dd_real::_nan = dd_real(qd::_d_nan, qd::_d_nan);
const dd_real dd_real::_inf = dd_real(qd::_d_inf, qd::_d_inf);

const double dd_real::_eps = 4.93038065763132e-32;  // 2^-104
const double dd_real::_min_normalized = 2.0041683600089728e-292;  // = 2^(-1022 + 53)
const dd_real dd_real::_max = 
    dd_real(1.79769313486231570815e+308, 9.97920154767359795037e+291);
const dd_real dd_real::_safe_max = 
    dd_real(1.7976931080746007281e+308, 9.97920154767359795037e+291);
const int dd_real::_ndigits = 31;

/* Table of sin(k * pi/16) and cos(k * pi/16). */
const dd_real dd_real::sin_table [] = {
  dd_real(1.950903220161282758e-01, -7.991079068461731263e-18), 
  dd_real(3.826834323650897818e-01, -1.005077269646158761e-17), 
  dd_real(5.555702330196021776e-01,  4.709410940561676821e-17),
  dd_real(7.071067811865475727e-01, -4.833646656726456726e-17)
};

const dd_real dd_real::cos_table [] = {
  dd_real(9.807852804032304306e-01, 1.854693999782500573e-17),
  dd_real(9.238795325112867385e-01, 1.764504708433667706e-17),
  dd_real(8.314696123025452357e-01, 1.407385698472802389e-18),
  dd_real(7.071067811865475727e-01, -4.833646656726456726e-17)
};

