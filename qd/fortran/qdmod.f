!  qdmod.f
!  
!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.
!
!  Copyright (c) 2000-2007
!
!  Fortran-90 module file to use with quad-double numbers.
!
!  Yozo Hida
!  David H Bailey    2007-03-15

module qdmodule
use ddmodule
  implicit none
  
  type qd_real
    sequence
    real*8 :: re(4)
  end type qd_real

  type qd_complex
    sequence
    real*8 :: cmp(8)
  end type qd_complex

  real*8 d_qd_eps
  parameter (d_qd_eps = 1.21543267145725d-63)

  type (qd_real) qd_one, qd_zero, qd_eps, qd_huge, qd_tiny
  parameter (qd_one = qd_real((/1.0d0, 0.0d0, 0.0d0, 0.0d0/)))
  parameter (qd_zero = qd_real((/0.0d0, 0.0d0, 0.0d0, 0.0d0/)))
  parameter (qd_eps = qd_real((/d_qd_eps, 0.0d0, 0.0d0, 0.0d0/)))
  parameter (qd_huge = qd_real((/ &
    1.79769313486231570815d+308, 9.97920154767359795037d+291, &
    5.53956966280111259858d+275, 3.07507889307840487279d+259/)))
  parameter (qd_tiny = qd_real((/3.25194908739046463067d-260, &
    0.0d0, 0.0d0, 0.0d0/)))

  interface assignment (=)
    module procedure assign_qd_str
    module procedure assign_qd
    module procedure assign_qd_d
    module procedure assign_d_qd
    module procedure assign_dd_qd
    module procedure assign_qd_dd
    module procedure assign_qd_i
    module procedure assign_i_qd
    module procedure assign_qdc
    module procedure assign_qdc_qd
    module procedure assign_qd_qdc
    module procedure assign_qdc_d
    module procedure assign_d_qdc
    module procedure assign_qdc_dc
    module procedure assign_dc_qdc
  end interface

  interface operator (+)
    module procedure add_qd
    module procedure add_qd_d
    module procedure add_d_qd
    module procedure add_qd_i
    module procedure add_i_qd
    module procedure add_qdc
    module procedure add_qdc_qd
    module procedure add_qd_qdc
    module procedure add_qdc_d
    module procedure add_d_qdc
  end interface

  interface operator (-)
    module procedure sub_qd
    module procedure sub_qd_d
    module procedure sub_d_qd
    module procedure neg_qd
    module procedure sub_qdc
    module procedure sub_qdc_qd
    module procedure sub_qd_qdc
    module procedure sub_qdc_d
    module procedure sub_d_qdc
    module procedure neg_qdc
  end interface

  interface operator (*)
    module procedure mul_qd
    module procedure mul_qd_d
    module procedure mul_d_qd
    module procedure mul_qd_i
    module procedure mul_i_qd
    module procedure mul_qdc
    module procedure mul_qdc_qd
    module procedure mul_qd_qdc
    module procedure mul_qdc_d
    module procedure mul_d_qdc
  end interface

  interface operator (/)
    module procedure div_qd
    module procedure div_qd_d
    module procedure div_d_qd
    module procedure div_qd_i
    module procedure div_i_qd
    module procedure div_qdc
    module procedure div_qdc_qd
    module procedure div_qd_qdc
    module procedure div_qdc_d
  end interface

  interface operator (**)
    module procedure pwr_qd
    module procedure pwr_qd_i
    module procedure pwr_d_qd
    module procedure pwr_qdc_i
  end interface

  interface qdreal
    module procedure to_qd_i
    module procedure to_qd_d
    module procedure to_qd_dd
    module procedure to_qd_qd
    module procedure to_qd_str
    module procedure to_qd_qdc
  end interface

  interface ddreal
    module procedure to_dd_qd
  end interface

  interface real
    module procedure to_d_qd
    module procedure to_d_qdc
  end interface

  interface qdcomplex
    module procedure to_qdc_qd
    module procedure to_qdc_d
    module procedure to_qdc_dc
  end interface

  interface int
    module procedure to_int_qd
  end interface

  interface sin
    module procedure qdsin
  end interface
  interface cos
    module procedure qdcos
  end interface
  interface tan
    module procedure qdtan
  end interface
  interface sincos
    module procedure qdsincos
  end interface
  interface qdcssnf
    module procedure qdcossin
  end interface

  interface asin
    module procedure qdasin
  end interface
  interface acos
    module procedure qdacos
  end interface
  interface atan
    module procedure qdatan
  end interface
  interface atan2
    module procedure qdatan2
  end interface

  interface exp
    module procedure qdexp
    module procedure qdcexp
  end interface
  interface log
    module procedure qdlog
    module procedure qdclog
  end interface
  interface log10
    module procedure qdlog10
  end interface

  interface sqrt
    module procedure qdsqrt
  end interface
  interface sqr
    module procedure qdsqr
  end interface
  interface nroot
    module procedure qdnroot
  end interface

  interface sinh
    module procedure qdsinh
  end interface
  interface cosh
    module procedure qdcosh
  end interface
  interface tanh
    module procedure qdtanh
  end interface
  interface sincosh
    module procedure qdsincosh
  end interface
  interface qdcsshf
    module procedure qdcossinh
  end interface

  interface asinh
    module procedure qdasinh
  end interface
  interface acosh
    module procedure qdacosh
  end interface
  interface atanh
    module procedure qdatanh
  end interface

  interface aint
    module procedure qdaint
  end interface

  interface nint
    module procedure qdnint
  end interface

  interface anint
    module procedure qdanint
  end interface

  interface abs
    module procedure qdabs
    module procedure qdcabs
  end interface

  interface sign
    module procedure qdsign
    module procedure qdsign_dd_d
  end interface

  interface random_number
    module procedure qdrand
  end interface

  interface operator (==)
    module procedure eq_qd
    module procedure eq_qd_d
    module procedure eq_d_qd
    module procedure eq_qd_i
    module procedure eq_i_qd
    module procedure eq_qdc
    module procedure eq_qdc_qd
    module procedure eq_qd_qdc
  end interface

  interface operator (/=)
    module procedure ne_qd
    module procedure ne_qd_d
    module procedure ne_d_qd
    module procedure ne_qd_i
    module procedure ne_i_qd
    module procedure ne_qdc
    module procedure ne_qdc_qd
    module procedure ne_qd_qdc
  end interface

  interface operator (>)
    module procedure gt_qd
    module procedure gt_qd_d
    module procedure gt_d_qd
    module procedure gt_qd_i
    module procedure gt_i_qd
  end interface

  interface operator (<)
    module procedure lt_qd
    module procedure lt_qd_d
    module procedure lt_d_qd
    module procedure lt_qd_i
    module procedure lt_i_qd
  end interface

  interface operator (>=)
    module procedure ge_qd
    module procedure ge_qd_d
    module procedure ge_d_qd
    module procedure ge_qd_i
    module procedure ge_i_qd
  end interface

  interface operator (<=)
    module procedure le_qd
    module procedure le_qd_d
    module procedure le_d_qd
    module procedure le_qd_i
    module procedure le_i_qd
  end interface

  interface qdread
    module procedure qdinpq
  end interface

  interface qdwrite
    module procedure qdoutq
  end interface

  interface qdcread
    module procedure qdcinpq
  end interface

  interface qdcwrite
    module procedure qdcoutq
  end interface

  interface dble
    module procedure to_d_qd
    module procedure to_d_qdc
  end interface

  interface cmplx
    module procedure to_dc_qdc
  end interface

  interface conjg
    module procedure qdcconjg
  end interface

  interface min
    module procedure qdmin
    module procedure qdmin2
  end interface
  interface max
    module procedure qdmax
    module procedure qdmax2
  end interface

  interface qdpi
    module procedure qd_pi
  end interface

  interface huge
    module procedure qdhuge
  end interface

  interface safe_huge
    module procedure qd_safe_huge
  end interface

  interface tiny
    module procedure qdtiny
  end interface

  interface epsilon
    module procedure qdepsilon
  end interface

  interface radix
    module procedure qd_radix
  end interface

  interface digits
    module procedure qd_digits
  end interface

  interface maxexponent
    module procedure qd_max_expn
  end interface

  interface minexponent
    module procedure qd_min_expn
  end interface

  interface precision
    module procedure qd_precision
  end interface

  interface range
    module procedure qd_range
  end interface

  interface nan
    module procedure qd_nan
  end interface

contains

! Assignments
  subroutine assign_qd_str(a, s)
    type (qd_real), intent(inout) :: a
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call qdinpc (t, a%re(1))
  end subroutine assign_qd_str

  subroutine assign_qd (a, b)
    type (qd_real), intent(inout) :: a
    type (qd_real), intent(in) :: b
    a%re(1) = b%re(1)
    a%re(2) = b%re(2)
    a%re(3) = b%re(3)
    a%re(4) = b%re(4)
  end subroutine assign_qd

  subroutine assign_qd_d(a, d)
    type (qd_real), intent(inout) :: a
    real*8, intent(in) :: d
    a%re(1) = d
    a%re(2) = 0.0d0
    a%re(3) = 0.0d0
    a%re(4) = 0.0d0
  end subroutine assign_qd_d

  subroutine assign_d_qd(d, a)
    real*8, intent(inout) :: d
    type (qd_real), intent(in) :: a
    d = a%re(1)
  end subroutine assign_d_qd

  subroutine assign_qd_i(a, i)
    type (qd_real), intent(inout) :: a
    integer, intent(in) :: i
    a%re(1) = i
    a%re(2) = 0.0d0
    a%re(3) = 0.0d0
    a%re(4) = 0.0d0
  end subroutine assign_qd_i

  subroutine assign_i_qd(i, a)
    integer, intent(inout) :: i
    type (qd_real), intent(in) :: a
    i = a%re(1)
  end subroutine assign_i_qd

  subroutine assign_dd_qd(dd, qd)
    type (dd_real), intent(inout) :: dd
    type (qd_real), intent(in) :: qd
    dd%re(1) = qd%re(1)
    dd%re(2) = qd%re(2)
  end subroutine assign_dd_qd

  subroutine assign_qd_dd(qd, dd)
    type (qd_real), intent(inout) :: qd
    type (dd_real), intent(in) :: dd
    qd%re(1) = dd%re(1)
    qd%re(2) = dd%re(2)
    qd%re(3) = 0.d0
    qd%re(4) = 0.d0
  end subroutine assign_qd_dd

  subroutine assign_qdc (a, b)
    type (qd_complex), intent(inout) :: a
    type (qd_complex), intent(in) :: b
    a%cmp(1) = b%cmp(1)
    a%cmp(2) = b%cmp(2)
    a%cmp(3) = b%cmp(3)
    a%cmp(4) = b%cmp(4)
    a%cmp(5) = b%cmp(5)
    a%cmp(6) = b%cmp(6)
    a%cmp(7) = b%cmp(7)
    a%cmp(8) = b%cmp(8)
  end subroutine assign_qdc

  subroutine assign_qdc_qd (qdc, qd)
    type (qd_complex), intent (inout) :: qdc
    type (qd_real), intent(in) :: qd
    qdc%cmp(1) = qd%re(1)
    qdc%cmp(2) = qd%re(2)
    qdc%cmp(3) = qd%re(3)
    qdc%cmp(4) = qd%re(4)
    qdc%cmp(5) = 0.d0
    qdc%cmp(6) = 0.d0
    qdc%cmp(7) = 0.d0
    qdc%cmp(8) = 0.d0
  end subroutine assign_qdc_qd

  subroutine assign_qd_qdc (qd, qdc)
    type (qd_real), intent (inout) :: qd
    type (qd_complex), intent(in) :: qdc
    qd%re(1) = qdc%cmp(1)
    qd%re(2) = qdc%cmp(2)
    qd%re(3) = qdc%cmp(3)
    qd%re(4) = qdc%cmp(4)
  end subroutine assign_qd_qdc

  subroutine assign_qdc_d (qdc, d)
    type (qd_complex), intent (inout) :: qdc
    real*8, intent(in) :: d
    qdc%cmp(1) = d
    qdc%cmp(2) = 0.d0
    qdc%cmp(3) = 0.d0
    qdc%cmp(4) = 0.d0
    qdc%cmp(5) = 0.d0
    qdc%cmp(6) = 0.d0
    qdc%cmp(7) = 0.d0
    qdc%cmp(8) = 0.d0
  end subroutine assign_qdc_d

  subroutine assign_d_qdc (d, qdc)
    real*8, intent(inout) :: d
    type (qd_complex), intent (in) :: qdc
    d = qdc%cmp(1)
  end subroutine assign_d_qdc

  subroutine assign_qdc_dc (qdc, dc)
    type (qd_complex), intent (inout) :: qdc
    complex (kind (0.d0)), intent (in) :: dc
    qdc%cmp(1) = dble (dc)
    qdc%cmp(2) = 0.d0
    qdc%cmp(3) = 0.d0
    qdc%cmp(4) = 0.d0
    qdc%cmp(5) = aimag (dc)
    qdc%cmp(6) = 0.d0
    qdc%cmp(7) = 0.d0
    qdc%cmp(8) = 0.d0
  end subroutine assign_qdc_dc

  subroutine assign_dc_qdc (dc, qdc)
    complex (kind (0.D0)), intent (inout) :: dc
    type (qd_complex), intent (in) :: qdc
    dc = cmplx (qdc%cmp(1), qdc%cmp(5), kind (0.d0))
  end subroutine assign_dc_qdc

! Conversions

  type (qd_real) function to_qd_i(ia)
    integer, intent(in) :: ia
    to_qd_i%re(1) = ia
    to_qd_i%re(2) = 0.d0
    to_qd_i%re(3) = 0.d0
    to_qd_i%re(4) = 0.d0
  end function to_qd_i

  type (qd_real) function to_qd_d(d)
    real*8, intent(in) :: d
    to_qd_d%re(1) = d
    to_qd_d%re(2) = 0.0d0
    to_qd_d%re(3) = 0.0d0
    to_qd_d%re(4) = 0.0d0
  end function to_qd_d

  real*8 function to_d_qd(qd) 
    type (qd_real), intent(in) :: qd
    to_d_qd = qd%re(1)
  end function to_d_qd

  integer function to_int_qd(a) 
    type (qd_real), intent(in) :: a
    to_int_qd = a%re(1)
  end function to_int_qd

  type (qd_real) function to_qd_dd (dd)
    type (dd_real), intent(in) :: dd
    to_qd_dd%re(1) = dd%re(1)
    to_qd_dd%re(2) = dd%re(2)
    to_qd_dd%re(3) = 0.d0
    to_qd_dd%re(4) = 0.d0
  end function to_qd_dd

  type (qd_real) function to_qd_qd (qd)
    type (qd_real), intent(in) :: qd
    to_qd_qd%re(1) = qd%re(1)
    to_qd_qd%re(2) = qd%re(2)
    to_qd_qd%re(3) = qd%re(3)
    to_qd_qd%re(4) = qd%re(4)
  end function to_qd_qd

  type (dd_real) function to_dd_qd (qd)
    type (qd_real), intent(in) :: qd
    to_dd_qd%re(1) = qd%re(1)
    to_dd_qd%re(2) = qd%re(2)
  end function to_dd_qd

  type (qd_real) function to_qd_str(s)
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call qdinpc (t, to_qd_str%re(1))
  end function to_qd_str

  type (qd_real) function to_qd_qdc(qdc)
    type (qd_complex), intent(in) :: qdc
    to_qd_qdc%re(1) = qdc%cmp(1)
    to_qd_qdc%re(2) = qdc%cmp(2)
    to_qd_qdc%re(3) = qdc%cmp(3)
    to_qd_qdc%re(4) = qdc%cmp(4)
  end function to_qd_qdc

  type (qd_complex) function to_qdc_qd(qd)
    type (qd_real), intent(in) :: qd
    to_qdc_qd%cmp(1) = qd%re(1)
    to_qdc_qd%cmp(2) = qd%re(2)
    to_qdc_qd%cmp(3) = qd%re(3)
    to_qdc_qd%cmp(4) = qd%re(4)
    to_qdc_qd%cmp(5) = 0.d0
    to_qdc_qd%cmp(6) = 0.d0
    to_qdc_qd%cmp(7) = 0.d0
    to_qdc_qd%cmp(8) = 0.d0
  end function to_qdc_qd

  type (qd_complex) function to_qdc_d(d)
    real*8, intent(in) :: d
    to_qdc_d%cmp(1) = d
    to_qdc_d%cmp(2) = 0.d0
    to_qdc_d%cmp(3) = 0.d0
    to_qdc_d%cmp(4) = 0.d0
    to_qdc_d%cmp(5) = 0.d0
    to_qdc_d%cmp(6) = 0.d0
    to_qdc_d%cmp(7) = 0.d0
    to_qdc_d%cmp(8) = 0.d0
  end function to_qdc_d

  complex (kind (0.D0)) function to_dc_qdc (qdc)
    type (qd_complex), intent (in) :: qdc
    to_dc_qdc = cmplx (qdc%cmp(1), qdc%cmp(5), kind (0.d0))
  end function to_dc_qdc

  type (qd_complex) function to_qdc_dc (dc)
    complex (kind (0.d0)), intent(in) :: dc
    to_qdc_dc%cmp(1) = dble (dc)
    to_qdc_dc%cmp(2) = 0.d0
    to_qdc_dc%cmp(3) = 0.d0
    to_qdc_dc%cmp(4) = 0.d0
    to_qdc_dc%cmp(5) = aimag (dc)
    to_qdc_dc%cmp(6) = 0.d0
    to_qdc_dc%cmp(7) = 0.d0
    to_qdc_dc%cmp(8) = 0.d0
  end function to_qdc_dc

  real*8 function to_d_qdc(qdc)
    type (qd_complex), intent(in) :: qdc
    to_d_qdc = qdc%cmp(1)
  end function to_d_qdc

!  Complex conjugation

  type (qd_complex) function qdcconjg (qdc)
    type (qd_complex), intent(in) :: qdc
    qdcconjg%cmp(1) = qdc%cmp(1)
    qdcconjg%cmp(2) = qdc%cmp(2)
    qdcconjg%cmp(3) = qdc%cmp(3)
    qdcconjg%cmp(4) = qdc%cmp(4)
    qdcconjg%cmp(5) = - qdc%cmp(5)
    qdcconjg%cmp(6) = - qdc%cmp(6)
    qdcconjg%cmp(7) = - qdc%cmp(7)
    qdcconjg%cmp(8) = - qdc%cmp(8)
  end function qdcconjg

! Additions
  type (qd_real) function add_qd(a, b)
    type (qd_real), intent(in) :: a, b
    call f_qd_add(a, b, add_qd)
  end function add_qd

  type (qd_real) function add_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_qd_add_qd_d(a, b, add_qd_d)
  end function add_qd_d

  type (qd_real) function add_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_add_d_qd(a, b, add_d_qd)
  end function add_d_qd

  type (qd_real) function add_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_add_d_qd(dble(a), b, add_i_qd)
  end function add_i_qd

  type (qd_real) function add_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    call f_qd_add_qd_d(a, dble(b), add_qd_i)
  end function add_qd_i

  type (qd_complex) function add_qdc(a, b)
    type (qd_complex), intent(in) :: a, b
    call f_qd_add (a%cmp(1), b%cmp(1), add_qdc%cmp(1))
    call f_qd_add (a%cmp(5), b%cmp(5), add_qdc%cmp(5))
  end function add_qdc

  type (qd_complex) function add_qdc_qd(a, b)
    type (qd_complex), intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_add (a%cmp(1), b%re(1), add_qdc_qd%cmp(1))
    add_qdc_qd%cmp(5) = a%cmp(5)
    add_qdc_qd%cmp(6) = a%cmp(6)
    add_qdc_qd%cmp(7) = a%cmp(7)
    add_qdc_qd%cmp(8) = a%cmp(8)
  end function add_qdc_qd

  type (qd_complex) function add_qd_qdc(a, b)
    type (qd_real), intent(in) :: a
    type (qd_complex), intent(in) :: b
    call f_qd_add (a%re(1), b%cmp(1), add_qd_qdc%cmp(1))
    add_qd_qdc%cmp(5) = b%cmp(5)
    add_qd_qdc%cmp(6) = b%cmp(6)
    add_qd_qdc%cmp(7) = b%cmp(7)
    add_qd_qdc%cmp(8) = b%cmp(8)
  end function add_qd_qdc

  type (qd_complex) function add_qdc_d(a, b)
    type (qd_complex), intent(in) :: a
    real*8, intent(in) :: b
    type (qd_real) :: qdb
    qdb%re(1) = b
    qdb%re(2) = 0.d0
    qdb%re(3) = 0.d0
    qdb%re(4) = 0.d0
    call f_qd_add (a%cmp(1), qdb%re(1), add_qdc_d%cmp(1))
    add_qdc_d%cmp(5) = a%cmp(5)
    add_qdc_d%cmp(6) = a%cmp(6)
    add_qdc_d%cmp(7) = a%cmp(7)
    add_qdc_d%cmp(8) = a%cmp(8)
  end function add_qdc_d

  type (qd_complex) function add_d_qdc(a, b)
    real*8, intent(in) :: a
    type (qd_complex), intent(in) :: b
    type (qd_real) qda
    qda%re(1) = a
    qda%re(2) = 0.d0
    qda%re(3) = 0.d0
    qda%re(4) = 0.d0
    call f_qd_add (qda%re(1), b%cmp(1), add_d_qdc%cmp(1))
    add_d_qdc%cmp(5) = b%cmp(5)
    add_d_qdc%cmp(6) = b%cmp(6)
    add_d_qdc%cmp(7) = b%cmp(7)
    add_d_qdc%cmp(8) = b%cmp(8)
  end function add_d_qdc

! Subtractions
  type (qd_real) function sub_qd(a, b)
    type (qd_real), intent(in) :: a, b
    call f_qd_sub(a, b, sub_qd)
  end function sub_qd

  type (qd_real) function sub_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_qd_sub_qd_d(a, b, sub_qd_d)
  end function sub_qd_d

  type (qd_real) function sub_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_sub_d_qd(a, b, sub_d_qd)
  end function sub_d_qd

  type (qd_complex) function sub_qdc(a, b)
    type (qd_complex), intent(in) :: a, b
    call f_qd_sub (a%cmp(1), b%cmp(1), sub_qdc%cmp(1))
    call f_qd_sub (a%cmp(5), b%cmp(5), sub_qdc%cmp(5))
  end function sub_qdc

  type (qd_complex) function sub_qdc_qd(a, b)
    type (qd_complex), intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_sub (a%cmp(1), b%re(1), sub_qdc_qd%cmp(1))
    sub_qdc_qd%cmp(5) = a%cmp(5)
    sub_qdc_qd%cmp(6) = a%cmp(6)
    sub_qdc_qd%cmp(7) = a%cmp(7)
    sub_qdc_qd%cmp(8) = a%cmp(8)
  end function sub_qdc_qd

  type (qd_complex) function sub_qd_qdc(a, b)
    type (qd_real), intent(in) :: a
    type (qd_complex), intent(in) :: b
    call f_qd_sub (a%re(1), b%cmp(1), sub_qd_qdc%cmp(1))
    sub_qd_qdc%cmp(5) = - b%cmp(5)
    sub_qd_qdc%cmp(6) = - b%cmp(6)
    sub_qd_qdc%cmp(7) = - b%cmp(7)
    sub_qd_qdc%cmp(8) = - b%cmp(8)
  end function sub_qd_qdc

  type (qd_complex) function sub_qdc_d(a, b)
    type (qd_complex), intent(in) :: a
    real*8, intent(in) :: b
    type (qd_real) qdb
    qdb%re(1) = b
    qdb%re(2) = 0.d0
    qdb%re(3) = 0.d0
    qdb%re(4) = 0.d0
    call f_qd_sub (a%cmp(1), qdb%re(1), sub_qdc_d%cmp(1))
    sub_qdc_d%cmp(5) = a%cmp(5)
    sub_qdc_d%cmp(6) = a%cmp(6)
    sub_qdc_d%cmp(7) = a%cmp(7)
    sub_qdc_d%cmp(8) = a%cmp(8)
  end function sub_qdc_d

  type (qd_complex) function sub_d_qdc(a, b)
    real*8, intent(in) :: a
    type (qd_complex), intent(in) :: b
    type (qd_real) qda
    qda%re(1) = a
    qda%re(2) = 0.d0
    qda%re(3) = 0.d0
    qda%re(4) = 0.d0
    call f_qd_sub (qda%re(1), b%cmp(1), sub_d_qdc%cmp(1))
    sub_d_qdc%cmp(5) = - b%cmp(5)
    sub_d_qdc%cmp(6) = - b%cmp(6)
    sub_d_qdc%cmp(7) = - b%cmp(7)
    sub_d_qdc%cmp(8) = - b%cmp(8)
  end function sub_d_qdc

! Unary Minus
  type (qd_real) function neg_qd(a)
    type (qd_real), intent(in) :: a
    call f_qd_neg(a, neg_qd)
  end function neg_qd

  type (qd_complex) function neg_qdc(a)
    type (qd_complex), intent(in) :: a
    integer j
    do j = 1, 8
      neg_qdc%cmp(j) = - a%cmp(j)
    enddo
  end function neg_qdc

! Multiplications
  type (qd_real) function mul_qd(a, b)
    type (qd_real), intent(in) :: a, b
    call f_qd_mul(a, b, mul_qd)
  end function mul_qd

  type (qd_real) function mul_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_qd_mul_qd_d(a, b, mul_qd_d)
  end function mul_qd_d

  type (qd_real) function mul_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_mul_d_qd(a, b, mul_d_qd)
  end function mul_d_qd

  type (qd_real) function mul_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    call f_qd_mul_qd_d(a, dble(b), mul_qd_i)
  end function mul_qd_i

  type (qd_real) function mul_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_mul_d_qd(dble(a), b, mul_i_qd)
  end function mul_i_qd

  type (qd_complex) function mul_qdc(a, b)
    type (qd_complex), intent(in) :: a, b
    type (qd_real) t1, t2
    call f_qd_mul (a%cmp(1), b%cmp(1), t1%re(1))
    call f_qd_mul (a%cmp(5), b%cmp(5), t2%re(1))
    call f_qd_sub (t1%re(1), t2%re(1), mul_qdc%cmp(1))
    call f_qd_mul (a%cmp(1), b%cmp(5), t1%re(1))
    call f_qd_mul (a%cmp(5), b%cmp(1), t2%re(1))
    call f_qd_add (t1%re(1), t2%re(1), mul_qdc%cmp(5))
  end function mul_qdc

  type (qd_complex) function mul_qdc_qd(a, b)
    type (qd_complex), intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_mul (a%cmp(1), b%re(1), mul_qdc_qd%cmp(1))
    call f_qd_mul (a%cmp(5), b%re(1), mul_qdc_qd%cmp(5))
  end function mul_qdc_qd

  type (qd_complex) function mul_qd_qdc(a, b)
    type (qd_real), intent(in) :: a
    type (qd_complex), intent(in) :: b
    call f_qd_mul (a%re(1), b%cmp(1), mul_qd_qdc%cmp(1))
    call f_qd_mul (a%re(1), b%cmp(5), mul_qd_qdc%cmp(5))
  end function mul_qd_qdc

  type (qd_complex) function mul_qdc_d(a, b)
    type (qd_complex), intent(in) :: a
    real*8, intent(in) :: b
    call f_qd_mul_qd_d (a%cmp(1), b, mul_qdc_d%cmp(1))
    call f_qd_mul_qd_d (a%cmp(5), b, mul_qdc_d%cmp(5))
  end function mul_qdc_d

  type (qd_complex) function mul_d_qdc(a, b)
    real*8, intent(in) :: a
    type (qd_complex), intent(in) :: b
    call f_qd_mul_d_qd (a, b%cmp(1), mul_d_qdc%cmp(1))
    call f_qd_mul_d_qd (a, b%cmp(5), mul_d_qdc%cmp(5))
  end function mul_d_qdc

! Divisions
  type (qd_real) function div_qd(a, b)
    type (qd_real), intent(in) :: a, b
    call f_qd_div(a, b, div_qd)
  end function div_qd

  type (qd_real) function div_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_qd_div_qd_d(a, b, div_qd_d)
  end function div_qd_d

  type (qd_real) function div_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_div_d_qd(a, b, div_d_qd)
  end function div_d_qd

  type (qd_real) function div_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    call f_qd_div_qd_d(a, dble(b), div_qd_i)
  end function div_qd_i

  type (qd_real) function div_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_div_d_qd(dble(a), b, div_i_qd)
  end function div_i_qd

  type (qd_complex) function div_qdc(a, b)
    type (qd_complex), intent(in) :: a, b
    type (qd_real) t1, t2, t3, t4, t5
    call f_qd_mul (a%cmp(1), b%cmp(1), t1%re(1))
    call f_qd_mul (a%cmp(5), b%cmp(5), t2%re(1))
    call f_qd_add (t1%re(1), t2%re(1), t3%re(1))
    call f_qd_mul (a%cmp(1), b%cmp(5), t1%re(1))
    call f_qd_mul (a%cmp(5), b%cmp(1), t2%re(1))
    call f_qd_sub (t2%re(1), t1%re(1), t4%re(1))
    call f_qd_mul (b%cmp(1), b%cmp(1), t1%re(1))
    call f_qd_mul (b%cmp(5), b%cmp(5), t2%re(1))
    call f_qd_add (t1%re(1), t2%re(1), t5%re(1))
    call f_qd_div (t3%re(1), t5%re(1), div_qdc%cmp(1))
    call f_qd_div (t4%re(1), t5%re(1), div_qdc%cmp(5))
  end function div_qdc

  type (qd_complex) function div_qdc_qd(a, b)
    type (qd_complex), intent(in) :: a
    type (qd_real), intent(in) :: b
    call f_qd_div (a%cmp(1), b%re(1), div_qdc_qd%cmp(1))
    call f_qd_div (a%cmp(5), b%re(1), div_qdc_qd%cmp(5))
  end function div_qdc_qd

  type (qd_complex) function div_qd_qdc(a, b)
    type (qd_real), intent(in) :: a
    type (qd_complex), intent(in) :: b
    type (qd_real) t1, t2, t3, t4, t5
    call f_qd_mul (a%re(1), b%cmp(1), t1%re(1))
    call f_qd_mul (a%re(1), b%cmp(5), t2%re(1))
    t2%re(1) = - t2%re(1)
    t2%re(2) = - t2%re(2)
    t2%re(3) = - t2%re(3)
    t2%re(4) = - t2%re(4)
    call f_qd_mul (b%cmp(1), b%cmp(1), t3%re(1))
    call f_qd_mul (b%cmp(5), b%cmp(5), t4%re(1))
    call f_qd_add (t3%re(1), t4%re(1), t5%re(1))
    call f_qd_div (t1%re(1), t5%re(1), div_qd_qdc%cmp(1))
    call f_qd_div (t2%re(1), t5%re(1), div_qd_qdc%cmp(5))
  end function div_qd_qdc

  type (qd_complex) function div_qdc_d(a,b)
    type (qd_complex), intent(in) :: a
    real*8, intent(in) :: b
    call f_qd_div_qd_d(a%cmp(1), b, div_qdc_d%cmp(1))
    call f_qd_div_qd_d(a%cmp(5), b, div_qdc_d%cmp(5))
  end function div_qdc_d

! Power
  type (qd_real) function pwr_qd (a, b)
    type (qd_real), intent(in) :: a, b
    type (qd_real) q1, q2
    call f_qd_log(a, q1)
    call f_qd_mul(q1, b, q2)
    call f_qd_exp(q2, pwr_qd)
  end function pwr_qd

  type (qd_real) function pwr_qd_i(a, n)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: n
    call f_qd_npwr(a, n, pwr_qd_i)
  end function pwr_qd_i

  type (qd_real) function pwr_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    type (qd_real) q1, q2, q3
    q1%re(1) = a
    q1%re(2) = 0.d0
    q1%re(3) = 0.d0
    q1%re(4) = 0.d0
    call f_qd_log(q1, q2)
    call f_qd_mul(q2, b, q3)
    call f_qd_exp(q3, pwr_d_qd)
  end function pwr_d_qd

  type (qd_complex) function pwr_qdc_i(a, n)
    type (qd_complex), intent(in) :: a
    integer, intent(in) :: n
    integer i, i2, j, n1
    type (qd_real) t1, t2, t3
    type (qd_complex) c1, c2

    if (n == 0) then
      if (a%cmp(1) == 0.d0 .and. a%cmp(2) == 0.d0 .and. &
        a%cmp(3) == 0.d0 .and. a%cmp(4) == 0.d0 .and. a%cmp(5) == 0.d0 .and. &
        a%cmp(6) == 0.d0 .and. a%cmp(7) == 0.d0 .and. a%cmp(8) == 0.d0) then
        write (6, *) 'pwr_qdc_i: a = 0 and n = 0'
        call f_qd_nan(pwr_qdc_i%cmp(1))
        call f_qd_nan(pwr_qdc_i%cmp(5))
        return
      endif
      pwr_qdc_i%cmp(1) = 1.d0
      pwr_qdc_i%cmp(2) = 0.d0
      pwr_qdc_i%cmp(3) = 0.d0
      pwr_qdc_i%cmp(4) = 0.d0
      pwr_qdc_i%cmp(5) = 0.d0
      pwr_qdc_i%cmp(6) = 0.d0
      pwr_qdc_i%cmp(7) = 0.d0
      pwr_qdc_i%cmp(8) = 0.d0
      return
    endif
    n1 = abs (n)
    i2 = 1
    do i = 1, 31
      i2 = 2 * i2
      if (i2 > n1) goto 100
    enddo

100 continue

    i2 = i2 / 2
    c1%cmp(1) = 1.d0
    c1%cmp(2) = 0.d0
    c1%cmp(3) = 0.d0
    c1%cmp(4) = 0.d0
    c1%cmp(5) = 0.d0
    c1%cmp(6) = 0.d0
    c1%cmp(7) = 0.d0
    c1%cmp(8) = 0.d0

110 continue

    if (n1 >= i2) then
      call f_qd_mul (a%cmp(1), c1%cmp(1), t1%re(1))
      call f_qd_mul (a%cmp(5), c1%cmp(5), t2%re(1))
      call f_qd_sub (t1%re(1), t2%re(1), c2%cmp(1))
      call f_qd_mul (a%cmp(1), c1%cmp(5), t1%re(1))
      call f_qd_mul (a%cmp(5), c1%cmp(1), t2%re(1))
      call f_qd_add (t1%re(1), t2%re(1), c2%cmp(5))
      do j = 1, 8
        c1%cmp(j) = c2%cmp(j)
      enddo
      n1 = n1 - i2
    endif
    i2 = i2 / 2
    if (i2 >= 1) then
      call f_qd_mul (c1%cmp(1), c1%cmp(1), t1%re(1))
      call f_qd_mul (c1%cmp(5), c1%cmp(5), t2%re(1))
      call f_qd_sub (t1%re(1), t2%re(1), c2%cmp(1))
      call f_qd_mul (c1%cmp(1), c1%cmp(5), t1%re(1))
      c2%cmp(5) = 2.d0 * t1%re(1)
      c2%cmp(6) = 2.d0 * t1%re(2)
      c2%cmp(7) = 2.d0 * t1%re(3)
      c2%cmp(8) = 2.d0 * t1%re(4)
      do j = 1, 8
        c1%cmp(j) = c2%cmp(j)
      enddo
      goto 110
    endif

    if (n > 0) then
      do j = 1, 8
        pwr_qdc_i%cmp(j) = c1%cmp(j)
      enddo
    else
      c1%cmp(5) = - c1%cmp(5)
      c1%cmp(6) = - c1%cmp(6)
      c1%cmp(7) = - c1%cmp(7)
      c1%cmp(8) = - c1%cmp(8)
      call f_qd_mul (c1%cmp(1), c1%cmp(1), t1%re(1))
      call f_qd_mul (c1%cmp(5), c1%cmp(5), t2%re(1))
      call f_qd_add (t1%re(1), t2%re(1), t3%re(1))
      call f_qd_div (c1%cmp(1), t3%re(1), pwr_qdc_i%cmp(1))
      call f_qd_div (c1%cmp(5), t3%re(1), pwr_qdc_i%cmp(5))
    endif

    return
  end function pwr_qdc_i


! Trigonometric Functions
  type (qd_real) function qdsin(a)
    type (qd_real), intent(in) :: a
    call f_qd_sin(a, qdsin)
  end function qdsin

  type (qd_real) function qdcos(a)
    type (qd_real), intent(in) :: a
    call f_qd_cos(a, qdcos)
  end function qdcos

  type (qd_real) function qdtan(a)
    type (qd_real), intent(in) :: a
    call f_qd_tan(a, qdtan)
  end function qdtan

  subroutine qdsincos(a, s, c)
    type (qd_real), intent(in) :: a
    type (qd_real), intent(out) :: s, c
    call f_qd_sincos(a, s, c)
  end subroutine qdsincos

  subroutine qdcossin(a, c, s)
    type (qd_real), intent(in) :: a
    type (qd_real), intent(out) :: s, c
    call f_qd_sincos(a, s, c)
  end subroutine qdcossin


! Inverse Trigonometric Functions
  type (qd_real) function qdasin(a)
    type (qd_real), intent(in) :: a
    call f_qd_asin(a, qdasin)
  end function qdasin

  type (qd_real) function qdacos(a)
    type (qd_real), intent(in) :: a
    call f_qd_acos(a, qdacos)
  end function qdacos

  type (qd_real) function qdatan(a)
    type (qd_real), intent(in) :: a
    call f_qd_atan(a, qdatan)
  end function qdatan

  type (qd_real) function qdatan2(a, b)
    type (qd_real), intent(in) :: a, b
    call f_qd_atan2(a, b, qdatan2)
  end function qdatan2

! Exponential and Logarithms
  type (qd_real) function qdexp(a)
    type (qd_real), intent(in) :: a
    call f_qd_exp(a, qdexp)
  end function qdexp

  type (qd_complex) function qdcexp (a)
    type (qd_complex), intent(in) :: a
    type (qd_real) t1, t2, t3
    call f_qd_exp (a%cmp(1), t1)
    call f_qd_sincos (a%cmp(5), t3, t2)
    call f_qd_mul (t1, t2, qdcexp%cmp(1))
    call f_qd_mul (t1, t3, qdcexp%cmp(5))
  end function qdcexp

  type (qd_real) function qdlog(a)
    type (qd_real), intent(in) :: a
    call f_qd_log(a, qdlog)
  end function qdlog

  type (qd_complex) function qdclog (a)
    type (qd_complex), intent(in) :: a
    type (qd_real) t1, t2, t3
    call f_qd_mul (a%cmp(1), a%cmp(1), t1)
    call f_qd_mul (a%cmp(5), a%cmp(5), t2)
    call f_qd_add (t1, t2, t3)
    call f_qd_log (t3, t1)
    qdclog%cmp(1) = 0.5d0 * t1%re(1)
    qdclog%cmp(2) = 0.5d0 * t1%re(2)
    qdclog%cmp(3) = 0.5d0 * t1%re(3)
    qdclog%cmp(4) = 0.5d0 * t1%re(4)
    call f_qd_atan2 (a%cmp(5), a%cmp(1), qdclog%cmp(5))
  end function qdclog

  type (qd_real) function qdlog10(a)
    type (qd_real), intent(in) :: a
    call f_qd_log10(a, qdlog10)
  end function qdlog10


! SQRT, etc.
  type (qd_real) function qdsqrt(a)
    type (qd_real), intent(in) :: a
    call f_qd_sqrt(a, qdsqrt)
  end function qdsqrt

  type (qd_real) function qdsqr(a)
    type (qd_real), intent(in) :: a
    call f_qd_sqr(a, qdsqr)
  end function qdsqr

  type (qd_real) function qdnroot(a, n)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: n
    call f_qd_nroot(a, n, qdnroot)
  end function qdnroot


! Hyperbolic Functions
  type (qd_real) function qdsinh(a)
    type (qd_real), intent(in) :: a
    call f_qd_sinh(a, qdsinh)
  end function qdsinh

  type (qd_real) function qdcosh(a)
    type (qd_real), intent(in) :: a
    call f_qd_cosh(a, qdcosh)
  end function qdcosh

  type (qd_real) function qdtanh(a)
    type (qd_real), intent(in) :: a
    call f_qd_tanh(a, qdtanh)
  end function qdtanh

  subroutine qdsincosh(a, s, c)
    type (qd_real), intent(in) :: a
    type (qd_real), intent(out) :: s, c
    call f_qd_sincosh(a, s, c)
  end subroutine qdsincosh

  subroutine qdcossinh(a, c, s)
    type (qd_real), intent(in) :: a
    type (qd_real), intent(out) :: s, c
    call f_qd_sincosh(a, s, c)
  end subroutine qdcossinh

! Inverse Hyperbolic Functions
  type (qd_real) function qdasinh(a)
    type (qd_real), intent(in) :: a
    call f_qd_asinh(a, qdasinh)
  end function qdasinh

  type (qd_real) function qdacosh(a)
    type (qd_real), intent(in) :: a
    call f_qd_acosh(a, qdacosh)
  end function qdacosh

  type (qd_real) function qdatanh(a)
    type (qd_real), intent(in) :: a
    call f_qd_atanh(a, qdatanh)
  end function qdatanh


! Rounding
  type (qd_real) function qdaint(a)
    type (qd_real), intent(in) :: a
    call f_qd_aint(a, qdaint)
  end function qdaint

  type (qd_real) function qdanint(a)
    type (qd_real), intent(in) :: a
    call f_qd_nint(a, qdanint)
  end function qdanint

  integer function qdnint(a)
    type (qd_real), intent(in) :: a
    qdnint = to_int_qd(qdaint(a));
  end function qdnint


! Random Number Generator
  subroutine qdrand(harvest)
    type (qd_real), intent(out) :: harvest
    call f_qd_rand(harvest)
  end subroutine qdrand


! Equality
  logical function eq_qd(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r == 0) then
      eq_qd = .true.
    else
      eq_qd = .false.
    end if
  end function eq_qd

  logical function eq_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_qd_comp_qd_d(a, b, r)
    if (r == 0) then
      eq_qd_d = .true.
    else
      eq_qd_d = .false.
    end if
  end function eq_qd_d

  logical function eq_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: r
    call f_qd_comp_d_qd(a, b, r)
    if (r == 0) then
      eq_d_qd = .true.
    else
      eq_d_qd = .false.
    end if
  end function eq_d_qd

  logical function eq_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    eq_qd_i = eq_qd_d(a, dble(b))
  end function eq_qd_i

  logical function eq_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    eq_i_qd = eq_d_qd(dble(a), b)
  end function eq_i_qd

  logical function eq_qdc (a, b)
    type (qd_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_qd_comp (a%cmp(1), b%cmp(1), i1)
    call f_qd_comp (a%cmp(5), b%cmp(5), i2)
    if (i1 == 0 .and. i2 == 0) then
      eq_qdc = .true.
    else
      eq_qdc = .false.
    endif
  end function eq_qdc

  logical function eq_qdc_qd (a, b)
    type (qd_complex), intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: i1
    call f_qd_comp (a%cmp(1), b%re(1), i1)
    if (i1 == 0 .and. a%cmp(5) == 0.d0 .and. a%cmp(6) == 0.d0 .and. &
      a%cmp(7) == 0.d0 .and. a%cmp(8) == 0.d0) then
      eq_qdc_qd = .true.
    else
      eq_qdc_qd = .false.
    endif
  end function eq_qdc_qd

  logical function eq_qd_qdc (a, b)
    type (qd_real), intent(in) :: a
    type (qd_complex), intent(in) :: b
    integer :: i1
    call f_qd_comp (a%re(1), b%cmp(1), i1)
    if (i1 == 0 .and. b%cmp(5) == 0.d0 .and. b%cmp(6) == 0.d0 .and. &
      b%cmp(7) == 0.d0 .and. b%cmp(8) == 0.d0) then
      eq_qd_qdc = .true.
    else
      eq_qd_qdc = .false.
    endif
  end function eq_qd_qdc


! Non-Equality
  logical function ne_qd(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r == 0) then
      ne_qd = .false.
    else
      ne_qd = .true.
    end if
  end function ne_qd

  logical function ne_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_qd_comp_qd_d(a, b, r)
    if (r == 0) then
      ne_qd_d = .false.
    else
      ne_qd_d = .true.
    end if
  end function ne_qd_d

  logical function ne_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: r
    call f_qd_comp_d_qd(a, b, r)
    if (r == 0) then
      ne_d_qd = .false.
    else
      ne_d_qd = .true.
    end if
  end function ne_d_qd

  logical function ne_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    ne_qd_i = ne_qd_d(a, dble(b))
  end function ne_qd_i

  logical function ne_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    ne_i_qd = ne_d_qd(dble(a), b)
  end function ne_i_qd

  logical function ne_qdc (a, b)
    type (qd_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_qd_comp (a%cmp(1), b%cmp(1), i1)
    call f_qd_comp (a%cmp(5), b%cmp(5), i2)
    if (i1 /= 0 .or. i2 /= 0) then
      ne_qdc = .true.
    else
      ne_qdc = .false.
    endif
  end function ne_qdc

  logical function ne_qdc_qd (a, b)
    type (qd_complex), intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: i1
    call f_qd_comp (a%cmp(1), b%re(1), i1)
    if (i1 /= 0 .or. a%cmp(5) /= 0.d0 .or. a%cmp(6) /= 0.d0 .or. &
      a%cmp(7) /= 0.d0 .or. a%cmp(8) /= 0.d0) then
      ne_qdc_qd = .true.
    else
      ne_qdc_qd = .false.
    endif
  end function ne_qdc_qd

  logical function ne_qd_qdc (a, b)
    type (qd_real), intent(in) :: a
    type (qd_complex), intent(in) :: b
    integer :: i1
    call f_qd_comp (a%re(1), b%cmp(1), i1)
    if (i1 /= 0 .or. b%cmp(5) /= 0.d0 .or. b%cmp(6) /= 0.d0 .or. &
      b%cmp(7) /= 0.d0 .or. b%cmp(8) /= 0.d0) then
      ne_qd_qdc = .true.
    else
      ne_qd_qdc = .false.
    endif
  end function ne_qd_qdc


! Greater-Than
  logical function gt_qd(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r == 1) then
      gt_qd = .true.
    else
      gt_qd = .false.
    end if
  end function gt_qd

  logical function gt_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_qd_comp_qd_d(a, b, r)
    if (r == 1) then
      gt_qd_d = .true.
    else
      gt_qd_d = .false.
    end if
  end function gt_qd_d

  logical function gt_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: r
    call f_qd_comp_d_qd(a, b, r)
    if (r == 1) then
      gt_d_qd = .true.
    else
      gt_d_qd = .false.
    end if
  end function gt_d_qd

  logical function gt_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    gt_qd_i = gt_qd_d(a, dble(b))
  end function gt_qd_i

  logical function gt_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    gt_i_qd = gt_d_qd(dble(a), b)
  end function gt_i_qd

! Less-Than
  logical function lt_qd(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r == -1) then
      lt_qd = .true.
    else
      lt_qd = .false.
    end if
  end function lt_qd

  logical function lt_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_qd_comp_qd_d(a, b, r)
    if (r == -1) then
      lt_qd_d = .true.
    else
      lt_qd_d = .false.
    end if
  end function lt_qd_d

  logical function lt_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: r
    call f_qd_comp_d_qd(a, b, r)
    if (r == -1) then
      lt_d_qd = .true.
    else
      lt_d_qd = .false.
    end if
  end function lt_d_qd

  logical function lt_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    lt_qd_i = lt_qd_d(a, dble(b))
  end function lt_qd_i

  logical function lt_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    lt_i_qd = lt_d_qd(dble(a), b)
  end function lt_i_qd

! Greater-Than-Or-Equal-To
  logical function ge_qd(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r >= 0) then
      ge_qd = .true.
    else
      ge_qd = .false.
    end if
  end function ge_qd

  logical function ge_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_qd_comp_qd_d(a, b, r)
    if (r >= 0) then
      ge_qd_d = .true.
    else
      ge_qd_d = .false.
    end if
  end function ge_qd_d

  logical function ge_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: r
    call f_qd_comp_d_qd(a, b, r)
    if (r >= 0) then
      ge_d_qd = .true.
    else
      ge_d_qd = .false.
    end if
  end function ge_d_qd

  logical function ge_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    ge_qd_i = ge_qd_d(a, dble(b))
  end function ge_qd_i

  logical function ge_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    ge_i_qd = ge_d_qd(dble(a), b)
  end function ge_i_qd

! Less-Than-Or-Equal-To
  logical function le_qd(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r <= 0) then
      le_qd = .true.
    else
      le_qd = .false.
    end if
  end function le_qd

  logical function le_qd_d(a, b)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_qd_comp_qd_d(a, b, r)
    if (r <= 0) then
      le_qd_d = .true.
    else
      le_qd_d = .false.
    end if
  end function le_qd_d

  logical function le_d_qd(a, b)
    real*8, intent(in) :: a
    type (qd_real), intent(in) :: b
    integer :: r
    call f_qd_comp_d_qd(a, b, r)
    if (r <= 0) then
      le_d_qd = .true.
    else
      le_d_qd = .false.
    end if
  end function le_d_qd

  logical function le_qd_i(a, b)
    type (qd_real), intent(in) :: a
    integer, intent(in) :: b
    le_qd_i = le_qd_d(a, dble(b))
  end function le_qd_i

  logical function le_i_qd(a, b)
    integer, intent(in) :: a
    type (qd_real), intent(in) :: b
    le_i_qd = le_d_qd(dble(a), b)
  end function le_i_qd


! Absolute Value
  type (qd_real) function qdabs(a)
    type (qd_real), intent(in) :: a
    call f_qd_abs(a, qdabs)
  end function qdabs

  type (qd_real) function qdcabs (qdc)
    type (qd_complex), intent(in) :: qdc
    type (qd_real) t1, t2, t3
    call f_qd_mul (qdc%cmp(1), qdc%cmp(1), t1%re(1))
    call f_qd_mul (qdc%cmp(5), qdc%cmp(5), t2%re(1))
    call f_qd_add (t1%re(1), t2%re(1), t3%re(1))
    call f_qd_sqrt (t3%re(1), qdcabs%re(1))
  end function qdcabs

! Sign transfer
  type (qd_real) function qdsign(a, b) result (c)
    type (qd_real), intent(in) :: a, b
    if (b%re(1) .gt. 0.0d0) then
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
        c%re(3) = a%re(3)
        c%re(4) = a%re(4)
      else
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
        c%re(3) = -a%re(3)
        c%re(4) = -a%re(4)
      end if
    else
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
        c%re(3) = -a%re(3)
        c%re(4) = -a%re(4)
      else
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
        c%re(3) = a%re(3)
        c%re(4) = a%re(4)
      end if
    endif
  end function qdsign

  type (qd_real) function qdsign_dd_d(a, b) result (c)
    type (qd_real), intent(in) :: a
    real*8, intent(in) :: b
    if (b .gt. 0.0d0) then
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
        c%re(3) = a%re(3)
        c%re(4) = a%re(4)
      else
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
        c%re(3) = -a%re(3)
        c%re(4) = -a%re(4)
      end if
    else
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
        c%re(3) = -a%re(3)
        c%re(4) = -a%re(4)
      else
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
        c%re(3) = a%re(3)
        c%re(4) = a%re(4)
      end if
    endif
  end function qdsign_dd_d

! Input
  subroutine qdinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qd_real), intent(in) :: q1
    type (qd_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qdinp (u, q1%re(1))

    if (present(q2)) then
      call qdinp (u, q2%re(1))
    end if

    if (present(q3)) then
      call qdinp (u, q3%re(1))
    end if

    if (present(q4)) then
      call qdinp (u, q4%re(1))
    end if

    if (present(q5)) then
      call qdinp (u, q5%re(1))
    end if

    if (present(q6)) then
      call qdinp (u, q6%re(1))
    end if

    if (present(q7)) then
      call qdinp (u, q7%re(1))
    end if

    if (present(q8)) then
      call qdinp (u, q8%re(1))
    end if

    if (present(q9)) then
      call qdinp (u, q9%re(1))
    end if

  end subroutine qdinpq

  subroutine qdcinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qd_complex), intent(in) :: q1
    type (qd_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qdinp (u, q1%cmp(1))
    call qdinp (u, q1%cmp(5))

    if (present(q2)) then
      call qdinp (u, q2%cmp(1))
      call qdinp (u, q2%cmp(5))
    end if

    if (present(q3)) then
      call qdinp (u, q3%cmp(1))
      call qdinp (u, q3%cmp(5))
    end if

    if (present(q4)) then
      call qdinp (u, q4%cmp(1))
      call qdinp (u, q4%cmp(5))
    end if

    if (present(q5)) then
      call qdinp (u, q5%cmp(1))
      call qdinp (u, q5%cmp(5))
    end if

    if (present(q6)) then
      call qdinp (u, q6%cmp(1))
      call qdinp (u, q6%cmp(5))
    end if

    if (present(q7)) then
      call qdinp (u, q7%cmp(1))
      call qdinp (u, q7%cmp(5))
    end if

    if (present(q8)) then
      call qdinp (u, q8%cmp(1))
      call qdinp (u, q8%cmp(5))
    end if

    if (present(q9)) then
      call qdinp (u, q9%cmp(1))
      call qdinp (u, q9%cmp(5))
    end if

  end subroutine qdcinpq

! Output
  subroutine qdoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qd_real), intent(in) :: q1
    type (qd_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qdout (u, q1%re(1))

    if (present(q2)) then
      call qdout (u, q2%re(1))
    end if

    if (present(q3)) then
      call qdout (u, q3%re(1))
    end if

    if (present(q4)) then
      call qdout (u, q4%re(1))
    end if

    if (present(q5)) then
      call qdout (u, q5%re(1))
    end if

    if (present(q6)) then
      call qdout (u, q6%re(1))
    end if

    if (present(q7)) then
      call qdout (u, q7%re(1))
    end if

    if (present(q8)) then
      call qdout (u, q8%re(1))
    end if

    if (present(q9)) then
      call qdout (u, q9%re(1))
    end if

  end subroutine qdoutq

  subroutine qdcoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qd_complex), intent(in) :: q1
    type (qd_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qdout (u, q1%cmp(1))
    call qdout (u, q1%cmp(5))

    if (present(q2)) then
      call qdout (u, q2%cmp(1))
      call qdout (u, q2%cmp(5))
    end if

    if (present(q3)) then
      call qdout (u, q3%cmp(1))
      call qdout (u, q3%cmp(5))
    end if

    if (present(q4)) then
      call qdout (u, q4%cmp(1))
      call qdout (u, q4%cmp(5))
    end if

    if (present(q5)) then
      call qdout (u, q5%cmp(1))
      call qdout (u, q5%cmp(5))
    end if

    if (present(q6)) then
      call qdout (u, q6%cmp(1))
      call qdout (u, q6%cmp(5))
    end if

    if (present(q7)) then
      call qdout (u, q7%cmp(1))
      call qdout (u, q7%cmp(5))
    end if

    if (present(q8)) then
      call qdout (u, q8%cmp(1))
      call qdout (u, q8%cmp(5))
    end if

    if (present(q9)) then
      call qdout (u, q9%cmp(1))
      call qdout (u, q9%cmp(5))
    end if

  end subroutine qdcoutq

  real*8 function qd_to_d(a)
    type (qd_real), intent(in) :: a
    qd_to_d = a%re(1)
  end function qd_to_d

  type (qd_real) function qdmin2(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r == 1) then
      qdmin2 = b
    else
      qdmin2 = a
    end if
  end function qdmin2

  type (qd_real) function qdmin(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (qd_real), intent(in) :: a1, a2, a3
    type (qd_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    qdmin = qdmin2(qdmin2(a1, a2), a3)
    if (present(a4)) qdmin = qdmin2(qdmin, a4)
    if (present(a5)) qdmin = qdmin2(qdmin, a5)
    if (present(a6)) qdmin = qdmin2(qdmin, a6)
    if (present(a7)) qdmin = qdmin2(qdmin, a7)
    if (present(a8)) qdmin = qdmin2(qdmin, a8)
    if (present(a9)) qdmin = qdmin2(qdmin, a9)
  end function qdmin

  type (qd_real) function qdmax2(a, b)
    type (qd_real), intent(in) :: a, b
    integer :: r
    call f_qd_comp(a, b, r)
    if (r == -1) then
      qdmax2 = b
    else
      qdmax2 = a
    end if
  end function qdmax2

  type (qd_real) function qdmax(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (qd_real), intent(in) :: a1, a2, a3
    type (qd_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    qdmax = qdmax2(qdmax2(a1, a2), a3)
    if (present(a4)) qdmax = qdmax2(qdmax, a4)
    if (present(a5)) qdmax = qdmax2(qdmax, a5)
    if (present(a6)) qdmax = qdmax2(qdmax, a6)
    if (present(a7)) qdmax = qdmax2(qdmax, a7)
    if (present(a8)) qdmax = qdmax2(qdmax, a8)
    if (present(a9)) qdmax = qdmax2(qdmax, a9)
  end function qdmax

  type (qd_real) function qd_pi()
    call f_qd_pi(qd_pi)
  end function qd_pi

subroutine qdinp (iu, a)

!   This routine reads the DD number A from logical unit IU.  The input
!   value must be placed on a single line of not more than 80 characters.

implicit none
integer iu, ln
parameter (ln = 80)
character*80 cs
real*8 a(4)

read (iu, '(a)', end = 100) cs
call qdinpc (cs, a)
goto 110

100 write (6, 1)
1  format ('*** qdinp: End-of-file encountered.')
! call qdabrt
stop

110 return

end subroutine

subroutine qdinpc (a, b)

!   Converts the CHARACTER*80 array A into the DD number B.

implicit none
integer i, id, ie, inz, ip, is, k, ln, lnn, beg
parameter (ln = 80)
real*8 bi
character*80 a
character*1 ai
character*10 dig
character*16 ca
parameter (dig = '0123456789')
real*8 b(4), f(4), s0(4), s1(4), s2(4)

id = 0
ip = -1
is = 0
inz = 0
s1(1) = 0.d0
s1(2) = 0.d0
s1(3) = 0.d0
s1(4) = 0.d0

beg = 0
do i = 1, 80
  if (a(i:i) /= ' ') then
    beg = i
    goto 80
  end if
end do

goto 210
80 continue

do i = beg, 80
  if (a(i:i) == ' ') then
    lnn = i-1
    goto 90
  end if
enddo

lnn = 80
90 continue

!   Scan for digits, looking for the period also.

do i = beg, lnn
  ai = a(i:i)
  if (ai .eq. '.') then
    if (ip >= 0) goto 210
    ip = id
    inz = 1
  elseif (ai .eq. '+') then
    if (id .ne. 0 .or. ip >= 0 .or. is .ne. 0) goto 210
    is = 1
  elseif (ai .eq. '-') then
    if (id .ne. 0 .or. ip >= 0 .or. is .ne. 0) goto 210
    is = -1
  elseif (ai .eq. 'e' .or. ai .eq. 'E' .or. ai .eq. 'd' .or. ai .eq. 'D') then
    goto 100
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
!    read (ai, '(f1.0)') bi
    bi = index (dig, ai) - 1
    if (inz > 0 .or. bi > 0.d0) then
      inz = 1
      id = id + 1
! call qdmuld (s1, 10.d0, s0)
      call f_qd_mul_qd_d (s1, 10.d0, s0)
      f(1) = bi
      f(2) = 0.d0
      f(3) = 0.d0
      f(4) = 0.d0
!    call qddqc (bi, f)
!    call qdadd (s0, f, s1)
      call f_qd_add (s0, f, s1)
    endif
  endif
enddo

100   continue
if (is .eq. -1) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
  s1(3) = - s1(3)
  s1(4) = - s1(4)
endif
k = i
if (ip == -1) ip = id
ie = 0
is = 0
ca = ' '

do i = k + 1, lnn
  ai = a(i:i)
  if (ai .eq. ' ') then
  elseif (ai .eq. '+') then
    if (ie .ne. 0 .or. is .ne. 0) goto 210
    is = 1
  elseif (ai .eq. '-') then
    if (ie .ne. 0 .or. is .ne. 0) goto 210
    is = -1
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
    ie = ie + 1
    if (ie .gt. 3) goto 210
    ca(ie:ie) = ai
  endif
enddo

! read (ca, '(i4)') ie
ie = dddigin (ca, 4)
if (is .eq. -1) ie = - ie
ie = ie + ip - id
s0(1) = 10.d0
s0(2) = 0.d0
s0(3) = 0.d0
s0(4) = 0.d0
! call qdnpwr (s0, ie, s2)
call f_qd_npwr (s0, ie, s2)
! call qdmul (s1, s2, b)
call f_qd_mul (s1, s2, b)
goto 220

210  write (6, 1) a
1 format ('*** qdinpc: Syntax error in literal string: ', a)
! call qdabrt
stop

220  return

end subroutine

subroutine qdout (iu, a)

!   This routine writes the QD number A on logical unit iu using a standard
!   E format, with lines 72 characters long.

implicit none
integer iu, ln
parameter (ln = 72)
character cs(72)
real*8 a(4)

call qdoutc (a, cs)
write (iu, '  (72a)') cs

return
end subroutine

subroutine qdoutc (a, b)
  implicit none
  real*8 a(4)
  character b(72)

  call f_qd_swrite(a, b)
end subroutine

type (qd_real) function qdhuge(a) 
  type (qd_real), intent(in) :: a
  qdhuge = qd_huge
end function qdhuge

type (qd_real) function qd_safe_huge(a)
  type (qd_real), intent(in) :: a
  qd_safe_huge = qd_real((/ &
    1.7976931080746007281d+308,  9.97920154767359795037d+291, &
    5.53956966280111259858d+275, 3.07507889307840487279d+259/))
end function qd_safe_huge

type (qd_real) function qdtiny(a) 
  type (qd_real), intent(in) :: a
  qdtiny = qd_tiny
end function qdtiny

type (qd_real) function qdepsilon(a) 
  type (qd_real), intent(in) :: a
  qdepsilon = qd_eps
end function qdepsilon

integer function qd_radix(a)
  type (qd_real), intent(in) :: a
  qd_radix = 2
end function qd_radix

integer function qd_digits(a)
  type (qd_real), intent(in) :: a
  qd_digits = 209
end function qd_digits

integer function qd_max_expn(a)
  type (qd_real), intent(in) :: a
  qd_max_expn = 1023
end function qd_max_expn

integer function qd_min_expn(a)
  type (qd_real), intent(in) :: a
  qd_min_expn = -863
end function qd_min_expn

integer function qd_precision(a)
  type (qd_real), intent(in) :: a
  qd_precision = 62
end function qd_precision

integer function qd_range(a)
  type (qd_real), intent(in) :: a
  qd_range = 259
end function qd_range

type (qd_real) function qd_nan(a)
  type (qd_real), intent(in) :: a
  call f_qd_nan(qd_nan)
end function qd_nan

end module qdmodule


