!  ddmod.f
!  
!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.
!
!  Copyright (c) 2000-2007
!
!  Fortran-90 module file to use with double-double numbers.
!
!  Yozo Hida
!  David H Bailey    2007-03-15

module ddmodule
  implicit none
  
  type dd_real
    sequence
    real*8 :: re(2)
  end type dd_real

  type dd_complex
    sequence
    real*8 :: cmp(4)
  end type dd_complex

  real*8 d_dd_eps
  parameter (d_dd_eps = 4.93038065763132d-32)

  type (dd_real) dd_one, dd_zero, dd_eps, dd_huge, dd_tiny
  parameter (dd_one = dd_real((/1.0d0, 0.0d0/)), &
             dd_zero = dd_real((/0.0d0, 0.0d0/)))
  parameter (dd_eps = dd_real((/d_dd_eps, 0.0d0/)))
  parameter (dd_huge = dd_real((/1.79769313486231570815d+308, &
                                 9.97920154767359795037d+291/)))
  parameter (dd_tiny = dd_real((/4.00833672001794555599d-292, 0.0d0/)))


  interface assignment (=)
    module procedure assign_dd_str
    module procedure assign_dd
    module procedure assign_dd_d
    module procedure assign_d_dd
    module procedure assign_dd_i
    module procedure assign_i_dd
    module procedure assign_ddc
    module procedure assign_ddc_dd
    module procedure assign_dd_ddc
    module procedure assign_ddc_d
    module procedure assign_d_ddc
    module procedure assign_ddc_dc
    module procedure assign_dc_ddc
  end interface

  interface operator (+)
    module procedure add_dd
    module procedure add_dd_d
    module procedure add_d_dd
    module procedure add_dd_i
    module procedure add_i_dd
    module procedure add_ddc
    module procedure add_ddc_dd
    module procedure add_dd_ddc
    module procedure add_ddc_d
    module procedure add_d_ddc
  end interface

  interface operator (-)
    module procedure sub_dd
    module procedure sub_dd_d
    module procedure sub_d_dd
    module procedure neg_dd
    module procedure sub_ddc
    module procedure sub_ddc_dd
    module procedure sub_dd_ddc
    module procedure sub_ddc_d
    module procedure sub_d_ddc
    module procedure neg_ddc
  end interface

  interface operator (*)
    module procedure mul_dd
    module procedure mul_dd_d
    module procedure mul_d_dd
    module procedure mul_dd_i
    module procedure mul_i_dd
    module procedure mul_ddc
    module procedure mul_ddc_dd
    module procedure mul_dd_ddc
    module procedure mul_ddc_d
    module procedure mul_d_ddc
  end interface

  interface operator (/)
    module procedure div_dd
    module procedure div_dd_d
    module procedure div_d_dd
    module procedure div_dd_i
    module procedure div_i_dd
    module procedure div_ddc
    module procedure div_ddc_dd
    module procedure div_dd_ddc
    module procedure div_ddc_d
  end interface

  interface operator (**)
    module procedure pwr_dd
    module procedure pwr_dd_i
    module procedure pwr_d_dd
    module procedure pwr_ddc_i
  end interface

  interface ddreal
    module procedure to_dd_i
    module procedure to_dd_d
    module procedure to_dd_dd
    module procedure to_dd_str
    module procedure to_dd_ddc
  end interface

  interface ddcomplex
     module procedure to_ddc_dd
     module procedure to_ddc_d
     module procedure to_ddc_dc
  end interface

  interface real
    module procedure to_d_dd
    module procedure to_d_ddc
  end interface

  interface int
    module procedure to_int_dd
  end interface

  interface sin
    module procedure ddsin
  end interface
  interface cos
    module procedure ddcos
  end interface
  interface tan
    module procedure ddtan
  end interface
  interface sincos
    module procedure ddsincos
  end interface
  interface ddcssnf
    module procedure ddcossin
  end interface

  interface asin
    module procedure ddasin
  end interface
  interface acos
    module procedure ddacos
  end interface
  interface atan
    module procedure ddatan
  end interface
  interface atan2
    module procedure ddatan2
  end interface

  interface exp
    module procedure ddexp
    module procedure ddcexp
  end interface
  interface log
    module procedure ddlog
    module procedure ddclog
  end interface
  interface log10
    module procedure ddlog10
  end interface

  interface sqrt
    module procedure ddsqrt
  end interface
  interface sqr
    module procedure ddsqr
  end interface
  interface nroot
    module procedure ddnroot
  end interface

  interface sinh
    module procedure ddsinh
  end interface
  interface cosh
    module procedure ddcosh
  end interface
  interface tanh
    module procedure ddtanh
  end interface
  interface sincosh
    module procedure ddsincosh
  end interface
  interface ddcsshf
    module procedure ddcossinh
  end interface

  interface asinh
    module procedure ddasinh
  end interface
  interface acosh
    module procedure ddacosh
  end interface
  interface atanh
    module procedure ddatanh
  end interface

  interface aint
    module procedure ddaint
  end interface

  interface anint
    module procedure ddanint
  end interface

  interface nint
    module procedure ddnint
  end interface

  interface abs
    module procedure ddabs
    module procedure ddcabs
  end interface

  interface sign
    module procedure ddsign
    module procedure ddsign_dd_d
  end interface

  interface random_number
    module procedure ddrand
  end interface

  interface operator (==)
    module procedure eq_dd
    module procedure eq_dd_d
    module procedure eq_d_dd
    module procedure eq_dd_i
    module procedure eq_i_dd
    module procedure eq_ddc
    module procedure eq_ddc_dd
    module procedure eq_dd_ddc
  end interface

  interface operator (/=)
    module procedure ne_dd
    module procedure ne_dd_d
    module procedure ne_d_dd
    module procedure ne_dd_i
    module procedure ne_i_dd
    module procedure ne_ddc
    module procedure ne_ddc_dd
    module procedure ne_dd_ddc
  end interface

  interface operator (>)
    module procedure gt_dd
    module procedure gt_dd_d
    module procedure gt_d_dd
    module procedure gt_dd_i
    module procedure gt_i_dd
  end interface

  interface operator (<)
    module procedure lt_dd
    module procedure lt_dd_d
    module procedure lt_d_dd
    module procedure lt_dd_i
    module procedure lt_i_dd
  end interface

  interface operator (>=)
    module procedure ge_dd
    module procedure ge_dd_d
    module procedure ge_d_dd
    module procedure ge_dd_i
    module procedure ge_i_dd
  end interface

  interface operator (<=)
    module procedure le_dd
    module procedure le_dd_d
    module procedure le_d_dd
    module procedure le_dd_i
    module procedure le_i_dd
  end interface

  interface ddread
    module procedure ddinpq
  end interface

  interface ddwrite
    module procedure ddoutq
  end interface

  interface ddcread
    module procedure ddcinpq
  end interface

  interface ddcwrite
    module procedure ddcoutq
  end interface

  interface dble
    module procedure to_d_dd
    module procedure to_d_ddc
  end interface

  interface cmplx
    module procedure to_dc_ddc
  end interface

  interface conjg
    module procedure ddcconjg
  end interface

  interface min
    module procedure ddmin
    module procedure ddmin2
  end interface
  interface max
    module procedure ddmax
    module procedure ddmax2
  end interface

  interface ddpi
    module procedure dd_pi
  end interface

  interface huge
    module procedure ddhuge
  end interface

  interface safe_huge
    module procedure dd_safe_huge
  end interface

  interface tiny
    module procedure ddtiny
  end interface

  interface epsilon
    module procedure ddepsilon
  end interface

  interface radix
    module procedure dd_radix
  end interface

  interface digits
    module procedure dd_digits
  end interface

  interface maxexponent
    module procedure dd_max_expn
  end interface

  interface minexponent
    module procedure dd_min_expn
  end interface

  interface nan
    module procedure dd_nan
  end interface

contains

! Assignments
  subroutine assign_dd_str(a, s)
    type (dd_real), intent(inout) :: a
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call ddinpc (t, a%re(1))
  end subroutine assign_dd_str

  subroutine assign_dd (a, b)
    type (dd_real), intent(inout) :: a
    type (dd_real), intent(in) :: b
    a%re(1) = b%re(1)
    a%re(2) = b%re(2)
  end subroutine assign_dd


  subroutine assign_dd_d(a, d)
    type (dd_real), intent(inout) :: a
    real*8, intent(in) :: d
    a%re(1) = d
    a%re(2) = 0.0d0
  end subroutine assign_dd_d

  subroutine assign_d_dd(d, a)
    real*8, intent(inout) :: d
    type (dd_real), intent(in) :: a
    d = a%re(1)
  end subroutine assign_d_dd

  subroutine assign_dd_i(a, i)
    type (dd_real), intent(inout) :: a
    integer, intent(in) :: i
    a%re(1) = i
    a%re(2) = 0.0d0
  end subroutine assign_dd_i

  subroutine assign_i_dd(i, a)
    integer, intent(inout) :: i
    type (dd_real), intent(in) :: a
    i = a%re(1)
  end subroutine assign_i_dd

  subroutine assign_ddc (a, b)
    type (dd_complex), intent(inout) :: a
    type (dd_complex), intent(in) :: b
    a%cmp(1) = b%cmp(1)
    a%cmp(2) = b%cmp(2)
    a%cmp(3) = b%cmp(3)
    a%cmp(4) = b%cmp(4)
  end subroutine assign_ddc

  subroutine assign_ddc_dd (ddc, dd)
    type (dd_complex), intent (inout) :: ddc
    type (dd_real), intent(in) :: dd
    ddc%cmp(1) = dd%re(1)
    ddc%cmp(2) = dd%re(2)
    ddc%cmp(3) = 0.d0
    ddc%cmp(4) = 0.d0
  end subroutine assign_ddc_dd

  subroutine assign_dd_ddc (dd, ddc)
    type (dd_real), intent (inout) :: dd
    type (dd_complex), intent(in) :: ddc
    dd%re(1) = ddc%cmp(1)
    dd%re(2) = ddc%cmp(2)
  end subroutine assign_dd_ddc

  subroutine assign_ddc_d (ddc, d)
    type (dd_complex), intent (inout) :: ddc
    real*8, intent(in) :: d
    ddc%cmp(1) = d
    ddc%cmp(2) = 0.d0
    ddc%cmp(3) = 0.d0
    ddc%cmp(4) = 0.d0
  end subroutine assign_ddc_d

  subroutine assign_d_ddc (d, ddc)
    real*8, intent(inout) :: d
    type (dd_complex), intent (in) :: ddc
    d = ddc%cmp(1)
  end subroutine assign_d_ddc

  subroutine assign_ddc_dc (ddc, dc)
    type (dd_complex), intent (inout) :: ddc
    complex (kind (0.d0)), intent (in) :: dc
    ddc%cmp(1) = dble (dc)
    ddc%cmp(2) = 0.d0
    ddc%cmp(3) = aimag (dc)
    ddc%cmp(4) = 0.d0
  end subroutine assign_ddc_dc

  subroutine assign_dc_ddc (dc, ddc)
    complex (kind (0.D0)), intent (inout) :: dc
    type (dd_complex), intent (in) :: ddc
    dc = cmplx (ddc%cmp(1), ddc%cmp(3), kind (0.d0))
  end subroutine assign_dc_ddc


! Conversions

  type (dd_real) function to_dd_i(ia)
    integer, intent(in) :: ia
    to_dd_i%re(1) = ia
    to_dd_i%re(2) = 0.d0
  end function to_dd_i

  type (dd_real) function to_dd_d(a)
    real*8, intent(in) :: a
    to_dd_d%re(1) = a
    to_dd_d%re(2) = 0.0d0
  end function to_dd_d

  type (dd_real) function to_dd_dd(a)
    type (dd_real), intent(in) :: a
    to_dd_dd%re(1) = a%re(1)
    to_dd_dd%re(2) = a%re(2)
  end function to_dd_dd

  real*8 function to_d_dd(a) 
    type (dd_real), intent(in) :: a
    to_d_dd = a%re(1)
  end function to_d_dd

  integer function to_int_dd(a) 
    type (dd_real), intent(in) :: a
    to_int_dd = a%re(1)
  end function to_int_dd

  type (dd_real) function to_dd_str(s)
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call ddinpc (t, to_dd_str%re(1))
  end function to_dd_str

  type (dd_real) function to_dd_ddc(ddc)
    type (dd_complex), intent(in) :: ddc
    to_dd_ddc%re(1) = ddc%cmp(1)
    to_dd_ddc%re(2) = ddc%cmp(2)
  end function to_dd_ddc

  type (dd_complex) function to_ddc_dd(dd)
    type (dd_real), intent(in) :: dd
    to_ddc_dd%cmp(1) = dd%re(1)
    to_ddc_dd%cmp(2) = dd%re(2)
    to_ddc_dd%cmp(3) = 0.d0
    to_ddc_dd%cmp(4) = 0.d0
  end function to_ddc_dd

  type (dd_complex) function to_ddc_d(d)
    real*8, intent(in) :: d
    to_ddc_d%cmp(1) = d
    to_ddc_d%cmp(2) = 0.d0
    to_ddc_d%cmp(3) = 0.d0
    to_ddc_d%cmp(4) = 0.d0
  end function to_ddc_d

  complex (kind (0.D0)) function to_dc_ddc (ddc)
    type (dd_complex), intent (in) :: ddc
    to_dc_ddc = cmplx (ddc%cmp(1), ddc%cmp(3), kind (0.d0))
  end function to_dc_ddc

  type (dd_complex) function to_ddc_dc (dc)
    complex (kind (0.d0)), intent(in) :: dc
    to_ddc_dc%cmp(1) = dble (dc)
    to_ddc_dc%cmp(2) = 0.d0
    to_ddc_dc%cmp(3) = aimag (dc)
    to_ddc_dc%cmp(4) = 0.d0
  end function to_ddc_dc

  real*8 function to_d_ddc(ddc)
    type (dd_complex), intent(in) :: ddc
    to_d_ddc = ddc%cmp(1)
  end function to_d_ddc

!  Complex conjugation
  type (dd_complex) function ddcconjg (ddc)
    type (dd_complex), intent(in) :: ddc
    ddcconjg%cmp(1) = ddc%cmp(1)
    ddcconjg%cmp(2) = ddc%cmp(2)
    ddcconjg%cmp(3) = - ddc%cmp(3)
    ddcconjg%cmp(4) = - ddc%cmp(4)
  end function ddcconjg

! Additions
  type (dd_real) function add_dd(a, b)
    type (dd_real), intent(in) :: a, b
    call f_dd_add(a, b, add_dd)
  end function add_dd

  type (dd_real) function add_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_dd_add_dd_d(a, b, add_dd_d)
  end function add_dd_d

  type (dd_real) function add_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_add_d_dd(a, b, add_d_dd)
  end function add_d_dd

  type (dd_real) function add_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_add_d_dd(dble(a), b, add_i_dd)
  end function add_i_dd

  type (dd_real) function add_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    call f_dd_add_dd_d(a, dble(b), add_dd_i)
  end function add_dd_i

  type (dd_complex) function add_ddc(a, b)
    type (dd_complex), intent(in) :: a, b
    call f_dd_add (a%cmp(1), b%cmp(1), add_ddc%cmp(1))
    call f_dd_add (a%cmp(3), b%cmp(3), add_ddc%cmp(3))
  end function add_ddc

  type (dd_complex) function add_ddc_d(a, b)
    type (dd_complex), intent(in) :: a
    real*8, intent(in) :: b
    type (dd_real) :: ddb
    ddb%re(1) = b
    ddb%re(2) = 0.d0
    call f_dd_add (a%cmp(1), ddb%re(1), add_ddc_d%cmp(1))
    add_ddc_d%cmp(3) = a%cmp(3)
    add_ddc_d%cmp(4) = a%cmp(4)
  end function add_ddc_d

  type (dd_complex) function add_d_ddc(a, b)
    real*8, intent(in) :: a
    type (dd_complex), intent(in) :: b
    type (dd_real) dda
    dda%re(1) = a
    dda%re(2) = 0.d0
    call f_dd_add (dda%re(1), b%cmp(1), add_d_ddc%cmp(1))
    add_d_ddc%cmp(3) = b%cmp(3)
    add_d_ddc%cmp(4) = b%cmp(4)
  end function add_d_ddc

  type (dd_complex) function add_ddc_dd(a, b)
    type (dd_complex), intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_add (a%cmp(1), b%re(1), add_ddc_dd%cmp(1))
    add_ddc_dd%cmp(3) = a%cmp(3)
    add_ddc_dd%cmp(4) = a%cmp(4)
  end function add_ddc_dd

  type (dd_complex) function add_dd_ddc(a, b)
    type (dd_real), intent(in) :: a
    type (dd_complex), intent(in) :: b
    call f_dd_add (a%re(1), b%cmp(1), add_dd_ddc%cmp(1))
    add_dd_ddc%cmp(3) = b%cmp(3)
    add_dd_ddc%cmp(4) = b%cmp(4)
  end function add_dd_ddc

! Subtractions
  type (dd_real) function sub_dd(a, b)
    type (dd_real), intent(in) :: a, b
    call f_dd_sub(a, b, sub_dd)
  end function sub_dd

  type (dd_real) function sub_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_dd_sub_dd_d(a, b, sub_dd_d)
  end function sub_dd_d

  type (dd_real) function sub_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_sub_d_dd(a, b, sub_d_dd)
  end function sub_d_dd

  type (dd_complex) function sub_ddc(a, b)
    type (dd_complex), intent(in) :: a, b
    call f_dd_sub (a%cmp(1), b%cmp(1), sub_ddc%cmp(1))
    call f_dd_sub (a%cmp(3), b%cmp(3), sub_ddc%cmp(3))
  end function sub_ddc

  type (dd_complex) function sub_ddc_d(a, b)
    type (dd_complex), intent(in) :: a
    real*8, intent(in) :: b
    type (dd_real) ddb
    ddb%re(1) = b
    ddb%re(2) = 0.d0
    call f_dd_sub (a%cmp(1), ddb%re(1), sub_ddc_d%cmp(1))
    sub_ddc_d%cmp(3) = a%cmp(3)
    sub_ddc_d%cmp(4) = a%cmp(4)
  end function sub_ddc_d

  type (dd_complex) function sub_d_ddc(a, b)
    real*8, intent(in) :: a
    type (dd_complex), intent(in) :: b
    type (dd_real) dda
    dda%re(1) = a
    dda%re(2) = 0.d0
    call f_dd_sub (dda%re(1), b%cmp(1), sub_d_ddc%cmp(1))
    sub_d_ddc%cmp(3) = - b%cmp(3)
    sub_d_ddc%cmp(4) = - b%cmp(4)
  end function sub_d_ddc

  type (dd_complex) function sub_ddc_dd(a, b)
    type (dd_complex), intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_sub (a%cmp(1), b%re(1), sub_ddc_dd%cmp(1))
    sub_ddc_dd%cmp(3) = a%cmp(3)
    sub_ddc_dd%cmp(4) = a%cmp(4)
  end function sub_ddc_dd

  type (dd_complex) function sub_dd_ddc(a, b)
    type (dd_real), intent(in) :: a
    type (dd_complex), intent(in) :: b
    call f_dd_sub (a%re(1), b%cmp(1), sub_dd_ddc%cmp(1))
    sub_dd_ddc%cmp(3) = - b%cmp(3)
    sub_dd_ddc%cmp(4) = - b%cmp(4)
  end function sub_dd_ddc

! Unary Minus
  type (dd_real) function neg_dd(a)
    type (dd_real), intent(in) :: a
    call f_dd_neg(a, neg_dd)
  end function neg_dd

  type (dd_complex) function neg_ddc(a)
    type (dd_complex), intent(in) :: a
    integer j
    do j = 1, 4
      neg_ddc%cmp(j) = - a%cmp(j)
    enddo
  end function neg_ddc

! Multiplications
  type (dd_real) function mul_dd(a, b)
    type (dd_real), intent(in) :: a, b
    call f_dd_mul(a, b, mul_dd)
  end function mul_dd

  type (dd_real) function mul_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_dd_mul_dd_d(a, b, mul_dd_d)
  end function mul_dd_d

  type (dd_real) function mul_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_mul_d_dd(a, b, mul_d_dd)
  end function mul_d_dd

  type (dd_real) function mul_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    call f_dd_mul_dd_d(a, dble(b), mul_dd_i)
  end function mul_dd_i

  type (dd_real) function mul_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_mul_d_dd(dble(a), b, mul_i_dd)
  end function mul_i_dd

  type (dd_complex) function mul_ddc(a, b)
    type (dd_complex), intent(in) :: a, b
    type (dd_real) t1, t2
    call f_dd_mul (a%cmp(1), b%cmp(1), t1%re(1))
    call f_dd_mul (a%cmp(3), b%cmp(3), t2%re(1))
    call f_dd_sub (t1%re(1), t2%re(1), mul_ddc%cmp(1))
    call f_dd_mul (a%cmp(1), b%cmp(3), t1%re(1))
    call f_dd_mul (a%cmp(3), b%cmp(1), t2%re(1))
    call f_dd_add (t1%re(1), t2%re(1), mul_ddc%cmp(3))
  end function mul_ddc

  type (dd_complex) function mul_ddc_d(a, b)
    type (dd_complex), intent(in) :: a
    real*8, intent(in) :: b
    call f_dd_mul_dd_d (a%cmp(1), b, mul_ddc_d%cmp(1))
    call f_dd_mul_dd_d (a%cmp(3), b, mul_ddc_d%cmp(3))
  end function mul_ddc_d

  type (dd_complex) function mul_d_ddc(a, b)
    real*8, intent(in) :: a
    type (dd_complex), intent(in) :: b
    call f_dd_mul_d_dd (a, b%cmp(1), mul_d_ddc%cmp(1))
    call f_dd_mul_d_dd (a, b%cmp(3), mul_d_ddc%cmp(3))
  end function mul_d_ddc

  type (dd_complex) function mul_ddc_dd(a, b)
    type (dd_complex), intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_mul (a%cmp(1), b%re(1), mul_ddc_dd%cmp(1))
    call f_dd_mul (a%cmp(3), b%re(1), mul_ddc_dd%cmp(3))
  end function mul_ddc_dd

  type (dd_complex) function mul_dd_ddc(a, b)
    type (dd_real), intent(in) :: a
    type (dd_complex), intent(in) :: b
    call f_dd_mul (a%re(1), b%cmp(1), mul_dd_ddc%cmp(1))
    call f_dd_mul (a%re(1), b%cmp(3), mul_dd_ddc%cmp(3))
  end function mul_dd_ddc

! Divisions
  type (dd_real) function div_dd(a, b)
    type (dd_real), intent(in) :: a, b
    call f_dd_div(a, b, div_dd)
  end function div_dd

  type (dd_real) function div_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    call f_dd_div_dd_d(a, b, div_dd_d)
  end function div_dd_d

  type (dd_real) function div_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_div_d_dd(a, b, div_d_dd)
  end function div_d_dd

  type (dd_real) function div_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    call f_dd_div_dd_d(a, dble(b), div_dd_i)
  end function div_dd_i

  type (dd_real) function div_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_div_d_dd(dble(a), b, div_i_dd)
  end function div_i_dd

  type (dd_complex) function div_ddc(a, b)
    type (dd_complex), intent(in) :: a, b
    type (dd_real) t1, t2, t3, t4, t5
    call f_dd_mul (a%cmp(1), b%cmp(1), t1%re(1))
    call f_dd_mul (a%cmp(3), b%cmp(3), t2%re(1))
    call f_dd_add (t1%re(1), t2%re(1), t3%re(1))
    call f_dd_mul (a%cmp(1), b%cmp(3), t1%re(1))
    call f_dd_mul (a%cmp(3), b%cmp(1), t2%re(1))
    call f_dd_sub (t2%re(1), t1%re(1), t4%re(1))
    call f_dd_mul (b%cmp(1), b%cmp(1), t1%re(1))
    call f_dd_mul (b%cmp(3), b%cmp(3), t2%re(1))
    call f_dd_add (t1%re(1), t2%re(1), t5%re(1))
    call f_dd_div (t3%re(1), t5%re(1), div_ddc%cmp(1))
    call f_dd_div (t4%re(1), t5%re(1), div_ddc%cmp(3))
  end function div_ddc

  type (dd_complex) function div_ddc_d(a,b)
    type (dd_complex), intent(in) :: a
    real*8, intent(in) :: b
    call f_dd_div_dd_d(a%cmp(1), b, div_ddc_d%cmp(1))
    call f_dd_div_dd_d(a%cmp(3), b, div_ddc_d%cmp(3))
  end function div_ddc_d

  type (dd_complex) function div_ddc_dd(a, b)
    type (dd_complex), intent(in) :: a
    type (dd_real), intent(in) :: b
    call f_dd_div (a%cmp(1), b%re(1), div_ddc_dd%cmp(1))
    call f_dd_div (a%cmp(3), b%re(1), div_ddc_dd%cmp(3))
  end function div_ddc_dd

  type (dd_complex) function div_dd_ddc(a, b)
    type (dd_real), intent(in) :: a
    type (dd_complex), intent(in) :: b
    type (dd_real) t1, t2, t3, t4, t5
    call f_dd_mul (a%re(1), b%cmp(1), t1%re(1))
    call f_dd_mul (a%re(1), b%cmp(3), t2%re(1))
    t2%re(1) = - t2%re(1)
    t2%re(2) = - t2%re(2)
    call f_dd_mul (b%cmp(1), b%cmp(1), t3%re(1))
    call f_dd_mul (b%cmp(3), b%cmp(3), t4%re(1))
    call f_dd_add (t3%re(1), t4%re(1), t5%re(1))
    call f_dd_div (t1%re(1), t5%re(1), div_dd_ddc%cmp(1))
    call f_dd_div (t2%re(1), t5%re(1), div_dd_ddc%cmp(3))
  end function div_dd_ddc

! Power
  type (dd_real) function pwr_dd (a, b)
    type (dd_real), intent(in) :: a, b
    type (dd_real) q1, q2
    call f_dd_log(a, q1)
    call f_dd_mul(q1, b, q2)
    call f_dd_exp(q2, pwr_dd)
  end function pwr_dd

  type (dd_real) function pwr_dd_i(a, n)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: n
    call f_dd_npwr(a, n, pwr_dd_i)
  end function pwr_dd_i

  type (dd_real) function pwr_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    type (dd_real) q1, q2, q3
    q1%re(1) = a
    q1%re(2) = 0.d0
    call f_dd_log(q1, q2)
    call f_dd_mul(q2, b, q3)
    call f_dd_exp(q3, pwr_d_dd)
  end function pwr_d_dd

  type (dd_complex) function pwr_ddc_i(a, n)
    type (dd_complex), intent(in) :: a
    integer, intent(in) :: n
    integer i, i2, j, n1
    type (dd_real) t1, t2, t3
    type (dd_complex) c1, c2

    if (n == 0) then
      if (a%cmp(1) == 0.d0 .and. a%cmp(2) == 0.d0 .and. &
        a%cmp(3) == 0.d0 .and. a%cmp(4) == 0.d0) then
        write (6, *) 'pwr_ddc_i: a = 0 and n = 0'
        call f_dd_nan(pwr_ddc_i%cmp(1))
        call f_dd_nan(pwr_ddc_i%cmp(3))
        return
      endif
      pwr_ddc_i%cmp(1) = 1.d0
      pwr_ddc_i%cmp(2) = 0.d0
      pwr_ddc_i%cmp(3) = 0.d0
      pwr_ddc_i%cmp(4) = 0.d0
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

110 continue

    if (n1 >= i2) then
      call f_dd_mul (a%cmp(1), c1%cmp(1), t1%re(1))
      call f_dd_mul (a%cmp(3), c1%cmp(3), t2%re(1))
      call f_dd_sub (t1%re(1), t2%re(1), c2%cmp(1))
      call f_dd_mul (a%cmp(1), c1%cmp(3), t1%re(1))
      call f_dd_mul (a%cmp(3), c1%cmp(1), t2%re(1))
      call f_dd_add (t1%re(1), t2%re(1), c2%cmp(3))
      do j = 1, 4
        c1%cmp(j) = c2%cmp(j)
      enddo
      n1 = n1 - i2
    endif
    i2 = i2 / 2
    if (i2 >= 1) then
      call f_dd_mul (c1%cmp(1), c1%cmp(1), t1%re(1))
      call f_dd_mul (c1%cmp(3), c1%cmp(3), t2%re(1))
      call f_dd_sub (t1%re(1), t2%re(1), c2%cmp(1))
      call f_dd_mul (c1%cmp(1), c1%cmp(3), t1%re(1))
      c2%cmp(3) = 2.d0 * t1%re(1)
      c2%cmp(4) = 2.d0 * t1%re(2)
      do j = 1, 4
        c1%cmp(j) = c2%cmp(j)
      enddo
      goto 110
    endif

    if (n > 0) then
      do j = 1, 4
        pwr_ddc_i%cmp(j) = c1%cmp(j)
      enddo
    else
      c1%cmp(3) = - c1%cmp(3)
      c1%cmp(4) = - c1%cmp(4)
      call f_dd_mul (c1%cmp(1), c1%cmp(1), t1%re(1))
      call f_dd_mul (c1%cmp(3), c1%cmp(3), t2%re(1))
      call f_dd_add (t1%re(1), t2%re(1), t3%re(1))
      call f_dd_div (c1%cmp(1), t3%re(1), pwr_ddc_i%cmp(1))
      call f_dd_div (c1%cmp(3), t3%re(1), pwr_ddc_i%cmp(3))
    endif

    return
  end function pwr_ddc_i

! Trigonometric Functions
  type (dd_real) function ddsin(a)
    type (dd_real), intent(in) :: a
    call f_dd_sin(a, ddsin)
  end function ddsin

  type (dd_real) function ddcos(a)
    type (dd_real), intent(in) :: a
    call f_dd_cos(a, ddcos)
  end function ddcos

  type (dd_real) function ddtan(a)
    type (dd_real), intent(in) :: a
    call f_dd_tan(a, ddtan)
  end function ddtan

  subroutine ddsincos(a, s, c)
    type (dd_real), intent(in) :: a
    type (dd_real), intent(out) :: s, c
    call f_dd_sincos(a, s, c)
  end subroutine ddsincos

  subroutine ddcossin(a, c, s)
    type (dd_real), intent(in) :: a
    type (dd_real), intent(out) :: s, c
    call f_dd_sincos(a, s, c)
  end subroutine ddcossin


! Inverse Trigonometric Functions
  type (dd_real) function ddasin(a)
    type (dd_real), intent(in) :: a
    call f_dd_asin(a, ddasin)
  end function ddasin

  type (dd_real) function ddacos(a)
    type (dd_real), intent(in) :: a
    call f_dd_acos(a, ddacos)
  end function ddacos

  type (dd_real) function ddatan(a)
    type (dd_real), intent(in) :: a
    call f_dd_atan(a, ddatan)
  end function ddatan

  type (dd_real) function ddatan2(a, b)
    type (dd_real), intent(in) :: a, b
    call f_dd_atan2(a, b, ddatan2)
  end function ddatan2

! Exponential and Logarithms
  type (dd_real) function ddexp(a)
    type (dd_real), intent(in) :: a
    call f_dd_exp(a, ddexp)
  end function ddexp

  type (dd_complex) function ddcexp (a)
    type (dd_complex), intent(in) :: a
    type (dd_real) t1, t2, t3
    call f_dd_exp (a%cmp(1), t1)
    call f_dd_sincos (a%cmp(3), t3, t2)
    call f_dd_mul (t1, t2, ddcexp%cmp(1))
    call f_dd_mul (t1, t3, ddcexp%cmp(3))
  end function ddcexp

  type (dd_real) function ddlog(a)
    type (dd_real), intent(in) :: a
    call f_dd_log(a, ddlog)
  end function ddlog

  type (dd_complex) function ddclog (a)
    type (dd_complex), intent(in) :: a
    type (dd_real) t1, t2, t3
    call f_dd_mul (a%cmp(1), a%cmp(1), t1)
    call f_dd_mul (a%cmp(3), a%cmp(3), t2)
    call f_dd_add (t1, t2, t3)
    call f_dd_log (t3, t1)
    ddclog%cmp(1) = 0.5d0 * t1%re(1)
    ddclog%cmp(2) = 0.5d0 * t1%re(2)
    call f_dd_atan2 (a%cmp(3), a%cmp(1), ddclog%cmp(3))
  end function ddclog


  type (dd_real) function ddlog10(a)
    type (dd_real), intent(in) :: a
    call f_dd_log10(a, ddlog10)
  end function ddlog10

! SQRT, etc.
  type (dd_real) function ddsqrt(a)
    type (dd_real), intent(in) :: a
    call f_dd_sqrt(a, ddsqrt)
  end function ddsqrt

  type (dd_real) function ddsqr(a)
    type (dd_real), intent(in) :: a
    call f_dd_sqr(a, ddsqr)
  end function ddsqr

  type (dd_real) function ddnroot(a, n)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: n
    call f_dd_nroot(a, n, ddnroot)
  end function ddnroot


! Hyperbolic Functions
  type (dd_real) function ddsinh(a)
    type (dd_real), intent(in) :: a
    call f_dd_sinh(a, ddsinh)
  end function ddsinh

  type (dd_real) function ddcosh(a)
    type (dd_real), intent(in) :: a
    call f_dd_cosh(a, ddcosh)
  end function ddcosh

  type (dd_real) function ddtanh(a)
    type (dd_real), intent(in) :: a
    call f_dd_tanh(a, ddtanh)
  end function ddtanh

  subroutine ddsincosh(a, s, c)
    type (dd_real), intent(in) :: a
    type (dd_real), intent(out) :: s, c
    call f_dd_sincosh(a, s, c)
  end subroutine ddsincosh

  subroutine ddcossinh(a, c, s)
    type (dd_real), intent(in) :: a
    type (dd_real), intent(out) :: s, c
    call f_dd_sincosh(a, s, c)
  end subroutine ddcossinh

! Inverse Hyperbolic Functions
  type (dd_real) function ddasinh(a)
    type (dd_real), intent(in) :: a
    call f_dd_asinh(a, ddasinh)
  end function ddasinh

  type (dd_real) function ddacosh(a)
    type (dd_real), intent(in) :: a
    call f_dd_acosh(a, ddacosh)
  end function ddacosh

  type (dd_real) function ddatanh(a)
    type (dd_real), intent(in) :: a
    call f_dd_atanh(a, ddatanh)
  end function ddatanh


! Rounding
  type (dd_real) function ddaint(a)
    type (dd_real), intent(in) :: a
    call f_dd_aint(a, ddaint)
  end function ddaint

  type (dd_real) function ddanint(a)
    type (dd_real), intent(in) :: a
    call f_dd_nint(a, ddanint)
  end function ddanint

  integer function ddnint(a)
    type (dd_real), intent(in) :: a
    ddnint = to_int_dd(ddaint(a));
  end function ddnint


! Random Number Generator
  subroutine ddrand(harvest)
    type (dd_real), intent(out) :: harvest
    call f_dd_rand(harvest)
  end subroutine ddrand


! Equality
  logical function eq_dd(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r == 0) then
      eq_dd = .true.
    else
      eq_dd = .false.
    end if
  end function eq_dd

  logical function eq_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_dd_comp_dd_d(a, b, r)
    if (r == 0) then
      eq_dd_d = .true.
    else
      eq_dd_d = .false.
    end if
  end function eq_dd_d

  logical function eq_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: r
    call f_dd_comp_d_dd(a, b, r)
    if (r == 0) then
      eq_d_dd = .true.
    else
      eq_d_dd = .false.
    end if
  end function eq_d_dd

  logical function eq_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    eq_dd_i = eq_dd_d(a, dble(b))
  end function eq_dd_i

  logical function eq_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    eq_i_dd = eq_d_dd(dble(a), b)
  end function eq_i_dd

  logical function eq_ddc (a, b)
    type (dd_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_dd_comp (a%cmp(1), b%cmp(1), i1)
    call f_dd_comp (a%cmp(3), b%cmp(3), i2)
    if (i1 == 0 .and. i2 == 0) then
      eq_ddc = .true.
    else
      eq_ddc = .false.
    endif
  end function eq_ddc

  logical function eq_ddc_dd (a, b)
    type (dd_complex), intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: i1
    call f_dd_comp (a%cmp(1), b%re(1), i1)
    if (i1 == 0 .and. a%cmp(3) == 0.d0 .and. a%cmp(4) == 0.d0) then
      eq_ddc_dd = .true.
    else
      eq_ddc_dd = .false.
    endif
  end function eq_ddc_dd

  logical function eq_dd_ddc (a, b)
    type (dd_real), intent(in) :: a
    type (dd_complex), intent(in) :: b
    integer :: i1
    call f_dd_comp (a%re(1), b%cmp(1), i1)
    if (i1 == 0 .and. b%cmp(3) == 0.d0 .and. b%cmp(4) == 0.d0) then
      eq_dd_ddc = .true.
    else
      eq_dd_ddc = .false.
    endif
  end function eq_dd_ddc


! Non-Equality
  logical function ne_dd(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r == 0) then
      ne_dd = .false.
    else
      ne_dd = .true.
    end if
  end function ne_dd

  logical function ne_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_dd_comp_dd_d(a, b, r)
    if (r == 0) then
      ne_dd_d = .false.
    else
      ne_dd_d = .true.
    end if
  end function ne_dd_d

  logical function ne_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: r
    call f_dd_comp_d_dd(a, b, r)
    if (r == 0) then
      ne_d_dd = .false.
    else
      ne_d_dd = .true.
    end if
  end function ne_d_dd

  logical function ne_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    ne_dd_i = ne_dd_d(a, dble(b))
  end function ne_dd_i

  logical function ne_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    ne_i_dd = ne_d_dd(dble(a), b)
  end function ne_i_dd

  logical function ne_ddc (a, b)
    type (dd_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_dd_comp (a%cmp(1), b%cmp(1), i1)
    call f_dd_comp (a%cmp(3), b%cmp(3), i2)
    if (i1 /= 0 .or. i2 /= 0) then
      ne_ddc = .true.
    else
      ne_ddc = .false.
    endif
  end function ne_ddc

  logical function ne_ddc_dd (a, b)
    type (dd_complex), intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: i1
    call f_dd_comp (a%cmp(1), b%re(1), i1)
    if (i1 /= 0 .or. a%cmp(3) /= 0.d0 .or. a%cmp(4) /= 0.d0) then
      ne_ddc_dd = .true.
    else
      ne_ddc_dd = .false.
    endif
  end function ne_ddc_dd

  logical function ne_dd_ddc (a, b)
    type (dd_real), intent(in) :: a
    type (dd_complex), intent(in) :: b
    integer :: i1
    call f_dd_comp (a%re(1), b%cmp(1), i1)
    if (i1 /= 0 .or. b%cmp(3) /= 0.d0 .or. b%cmp(4) /= 0.d0) then
      ne_dd_ddc = .true.
    else
      ne_dd_ddc = .false.
    endif
  end function ne_dd_ddc


! Greater-Than
  logical function gt_dd(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r == 1) then
      gt_dd = .true.
    else
      gt_dd = .false.
    end if
  end function gt_dd

  logical function gt_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_dd_comp_dd_d(a, b, r)
    if (r == 1) then
      gt_dd_d = .true.
    else
      gt_dd_d = .false.
    end if
  end function gt_dd_d

  logical function gt_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: r
    call f_dd_comp_d_dd(a, b, r)
    if (r == 1) then
      gt_d_dd = .true.
    else
      gt_d_dd = .false.
    end if
  end function gt_d_dd

  logical function gt_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    gt_dd_i = gt_dd_d(a, dble(b))
  end function gt_dd_i

  logical function gt_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    gt_i_dd = gt_d_dd(dble(a), b)
  end function gt_i_dd


! Less-Than
  logical function lt_dd(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r == -1) then
      lt_dd = .true.
    else
      lt_dd = .false.
    end if
  end function lt_dd

  logical function lt_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_dd_comp_dd_d(a, b, r)
    if (r == -1) then
      lt_dd_d = .true.
    else
      lt_dd_d = .false.
    end if
  end function lt_dd_d

  logical function lt_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: r
    call f_dd_comp_d_dd(a, b, r)
    if (r == -1) then
      lt_d_dd = .true.
    else
      lt_d_dd = .false.
    end if
  end function lt_d_dd

  logical function lt_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    lt_dd_i = lt_dd_d(a, dble(b))
  end function lt_dd_i

  logical function lt_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    lt_i_dd = lt_d_dd(dble(a), b)
  end function lt_i_dd

! Greater-Than-Or-Equal-To
  logical function ge_dd(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r >= 0) then
      ge_dd = .true.
    else
      ge_dd = .false.
    end if
  end function ge_dd

  logical function ge_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_dd_comp_dd_d(a, b, r)
    if (r >= 0) then
      ge_dd_d = .true.
    else
      ge_dd_d = .false.
    end if
  end function ge_dd_d

  logical function ge_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: r
    call f_dd_comp_d_dd(a, b, r)
    if (r >= 0) then
      ge_d_dd = .true.
    else
      ge_d_dd = .false.
    end if
  end function ge_d_dd

  logical function ge_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    ge_dd_i = ge_dd_d(a, dble(b))
  end function ge_dd_i

  logical function ge_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    ge_i_dd = ge_d_dd(dble(a), b)
  end function ge_i_dd

! Less-Than-Or-Equal-To
  logical function le_dd(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r <= 0) then
      le_dd = .true.
    else
      le_dd = .false.
    end if
  end function le_dd

  logical function le_dd_d(a, b)
    type (dd_real), intent(in) :: a
    real*8, intent(in) :: b
    integer :: r
    call f_dd_comp_dd_d(a, b, r)
    if (r <= 0) then
      le_dd_d = .true.
    else
      le_dd_d = .false.
    end if
  end function le_dd_d

  logical function le_d_dd(a, b)
    real*8, intent(in) :: a
    type (dd_real), intent(in) :: b
    integer :: r
    call f_dd_comp_d_dd(a, b, r)
    if (r <= 0) then
      le_d_dd = .true.
    else
      le_d_dd = .false.
    end if
  end function le_d_dd

  logical function le_dd_i(a, b)
    type (dd_real), intent(in) :: a
    integer, intent(in) :: b
    le_dd_i = le_dd_d(a, dble(b))
  end function le_dd_i

  logical function le_i_dd(a, b)
    integer, intent(in) :: a
    type (dd_real), intent(in) :: b
    le_i_dd = le_d_dd(dble(a), b)
  end function le_i_dd

! Absolute Value
  type (dd_real) function ddabs(a)
    type (dd_real), intent(in) :: a
    call f_dd_abs(a, ddabs)
  end function ddabs

  type (dd_real) function ddcabs (ddc)
    type (dd_complex), intent(in) :: ddc
    type (dd_real) t1, t2, t3
    call f_dd_mul (ddc%cmp(1), ddc%cmp(1), t1%re(1))
    call f_dd_mul (ddc%cmp(3), ddc%cmp(3), t2%re(1))
    call f_dd_add (t1%re(1), t2%re(1), t3%re(1))
    call f_dd_sqrt (t3%re(1), ddcabs%re(1))
  end function ddcabs

! Sign transfer
  type (dd_real) function ddsign(a, b) result (c)
    type (dd_real), intent(in) :: a, b
    if (b%re(1) .gt. 0.0d0) then
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
      else
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
      end if
    else
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
      else
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
      end if
    endif
  end function ddsign

  type (dd_real) function ddsign_dd_d(a, b) result (c)
    type (dd_real), intent(in) :: a
    real*8, intent(in) ::  b
    if (b .gt. 0.0d0) then
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
      else
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
      end if
    else
      if (a%re(1) .gt. 0.0d0) then
        c%re(1) = -a%re(1)
        c%re(2) = -a%re(2)
      else
        c%re(1) = a%re(1)
        c%re(2) = a%re(2)
      end if
    endif
  end function ddsign_dd_d

! Input
  subroutine ddinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (dd_real), intent(in) :: q1
    type (dd_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9
!    CHARACTER (LEN=72) :: str

    call ddinp (u, q1%re(1))

    if (present(q2)) then
      call ddinp (u, q2%re(1))
    end if

    if (present(q3)) then
      call ddinp (u, q3%re(1))
    end if

    if (present(q4)) then
      call ddinp (u, q4%re(1))
    end if

    if (present(q5)) then
      call ddinp (u, q5%re(1))
    end if

    if (present(q6)) then
      call ddinp (u, q6%re(1))
    end if

    if (present(q7)) then
      call ddinp (u, q7%re(1))
    end if

    if (present(q8)) then
      call ddinp (u, q8%re(1))
    end if

    if (present(q9)) then
      call ddinp (u, q9%re(1))
   end if

  end subroutine ddinpq

  subroutine ddcinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (dd_complex), intent(in) :: q1
    type (dd_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call ddinp (u, q1%cmp(1))
    call ddinp (u, q1%cmp(3))

    if (present(q2)) then
      call ddinp (u, q2%cmp(1))
      call ddinp (u, q2%cmp(3))
    end if

    if (present(q3)) then
      call ddinp (u, q3%cmp(1))
      call ddinp (u, q3%cmp(3))
    end if

    if (present(q4)) then
      call ddinp (u, q4%cmp(1))
      call ddinp (u, q4%cmp(3))
    end if

    if (present(q5)) then
      call ddinp (u, q5%cmp(1))
      call ddinp (u, q5%cmp(3))
    end if

    if (present(q6)) then
      call ddinp (u, q6%cmp(1))
      call ddinp (u, q6%cmp(3))
    end if

    if (present(q7)) then
      call ddinp (u, q7%cmp(1))
      call ddinp (u, q7%cmp(3))
    end if

    if (present(q8)) then
      call ddinp (u, q8%cmp(1))
      call ddinp (u, q8%cmp(3))
    end if

    if (present(q9)) then
      call ddinp (u, q9%cmp(1))
      call ddinp (u, q9%cmp(3))
    end if

  end subroutine ddcinpq


! Output
  subroutine ddoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (dd_real), intent(in) :: q1
    type (dd_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9
!    CHARACTER (LEN=72) :: str

    call ddout (u, q1%re(1))

    if (present(q2)) then
      call ddout (u, q2%re(1))
    end if

    if (present(q3)) then
      call ddout (u, q3%re(1))
    end if

    if (present(q4)) then
      call ddout (u, q4%re(1))
    end if

    if (present(q5)) then
      call ddout (u, q5%re(1))
    end if

    if (present(q6)) then
      call ddout (u, q6%re(1))
    end if

    if (present(q7)) then
      call ddout (u, q7%re(1))
    end if

    if (present(q8)) then
      call ddout (u, q8%re(1))
    end if

    if (present(q9)) then
      call ddout (u, q9%re(1))
   end if

  end subroutine ddoutq

  subroutine ddcoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (dd_complex), intent(in) :: q1
    type (dd_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call ddout (u, q1%cmp(1))
    call ddout (u, q1%cmp(3))

    if (present(q2)) then
      call ddout (u, q2%cmp(1))
      call ddout (u, q2%cmp(3))
    end if

    if (present(q3)) then
      call ddout (u, q3%cmp(1))
      call ddout (u, q3%cmp(3))
    end if

    if (present(q4)) then
      call ddout (u, q4%cmp(1))
      call ddout (u, q4%cmp(3))
    end if

    if (present(q5)) then
      call ddout (u, q5%cmp(1))
      call ddout (u, q5%cmp(3))
    end if

    if (present(q6)) then
      call ddout (u, q6%cmp(1))
      call ddout (u, q6%cmp(3))
    end if

    if (present(q7)) then
      call ddout (u, q7%cmp(1))
      call ddout (u, q7%cmp(3))
    end if

    if (present(q8)) then
      call ddout (u, q8%cmp(1))
      call ddout (u, q8%cmp(3))
    end if

    if (present(q9)) then
      call ddout (u, q9%cmp(1))
      call ddout (u, q9%cmp(3))
    end if

  end subroutine ddcoutq

  type (dd_real) function ddmin2(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r == 1) then
      ddmin2 = b
    else
      ddmin2 = a
    end if
  end function ddmin2

  type (dd_real) function ddmin(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (dd_real), intent(in) :: a1, a2, a3
    type (dd_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    ddmin = ddmin2(ddmin2(a1, a2), a3)
    if (present(a4)) ddmin = ddmin2(ddmin, a4)
    if (present(a5)) ddmin = ddmin2(ddmin, a5)
    if (present(a6)) ddmin = ddmin2(ddmin, a6)
    if (present(a7)) ddmin = ddmin2(ddmin, a7)
    if (present(a8)) ddmin = ddmin2(ddmin, a8)
    if (present(a9)) ddmin = ddmin2(ddmin, a9)
  end function ddmin

  type (dd_real) function ddmax2(a, b)
    type (dd_real), intent(in) :: a, b
    integer :: r
    call f_dd_comp(a, b, r)
    if (r == -1) then
      ddmax2 = b
    else
      ddmax2 = a
    end if
  end function ddmax2

  type (dd_real) function ddmax(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (dd_real), intent(in) :: a1, a2, a3
    type (dd_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    ddmax = ddmax2(ddmax2(a1, a2), a3)
    if (present(a4)) ddmax = ddmax2(ddmax, a4)
    if (present(a5)) ddmax = ddmax2(ddmax, a5)
    if (present(a6)) ddmax = ddmax2(ddmax, a6)
    if (present(a7)) ddmax = ddmax2(ddmax, a7)
    if (present(a8)) ddmax = ddmax2(ddmax, a8)
    if (present(a9)) ddmax = ddmax2(ddmax, a9)
  end function ddmax

  type (dd_real) function dd_pi()
    call f_dd_pi(dd_pi)
  end function dd_pi

subroutine ddinp (iu, a)

!   This routine reads the DD number A from logical unit IU.  The input
!   value must be placed on a single line of not more than 80 characters.

implicit none
integer iu, ln
parameter (ln = 80)
character*80 cs
real*8 a(2)

read (iu, '(a)', end = 100) cs
call ddinpc (cs, a)
goto 110

100 write (6, 1)
1  format ('*** ddinp: End-of-file encountered.')
! call ddabrt
stop

110 return
end subroutine

subroutine ddinpc (a, b)

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
real*8 b(2), f(2), s0(2), s1(2), s2(2)

id = 0
ip = -1
is = 0
inz = 0
s1(1) = 0.d0
s1(2) = 0.d0

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
!    call ddmuld (s1, 10.d0, s0)
      call f_dd_mul_dd_d (s1, 10.d0, s0)
      f(1) = bi
      f(2) = 0.d0
!    call dddqc (bi, f)
!    call ddadd (s0, f, s1)
      call f_dd_add (s0, f, s1)
    endif
  endif
enddo

100   continue
if (is .eq. -1) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
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
! call ddnpwr (s0, ie, s2)
call f_dd_npwr (s0, ie, s2)
! call ddmul (s1, s2, b)
call f_dd_mul (s1, s2, b)
goto 220

210  write (6, 1) a
1 format ('*** ddinpc: Syntax error in literal string: ', a)
! call ddabrt
stop

220  return
end subroutine

subroutine ddout (iu, a)

!   This routine writes the DD number A on logical unit iu using a standard
!   E format, with lines 40 characters long.

implicit none
integer iu, ln
parameter (ln = 40)
character cs(40)
real*8 a(2)

call ddoutc (a, cs)
write (iu, '(40a)') cs

return
end subroutine

subroutine ddoutc (a, b)
  implicit none
  real*8 a(2)
  character b(40)

  call f_dd_swrite(a, b)
end subroutine

  real*8 function dddigin (ca, n)
    implicit none
    real*8 d1
    character*(*), ca
    character*16 digits
    integer i, k, n
    parameter (digits = '0123456789')

    d1 = 0.d0

    do i = 1, n
      k = index (digits, ca(i:i)) - 1
      if (k < 0) then
        write (6, *) 'dddigin: non-digit in character string'
      elseif (k <= 9) then
        d1 = 10.d0 * d1 + k
      endif
    enddo

    dddigin = d1
  end function

  character*16 function dddigout (a, n)
    implicit none
    real*8 a, d1, d2
    character*16 ca, digits
    parameter (digits = '0123456789')
    integer i, is, k, n
    real*8 abs, aint

    ca = ' '
    is = sign (1.d0, a)
    d1 = abs (a)

    do i = n, 1, -1
      d2 = aint (d1 / 10.d0)
      k = 1.d0 + (d1 - 10.d0 * d2)
      d1 = d2
      ca(i:i) = digits(k:k)
      if (d1 == 0.d0) goto 100
    enddo

    i = 0

100 continue

    if (is < 0 .and. i > 1) then
      ca(i-1:i-1) = '-'
    elseif (i == 0 .or. is < 0 .and. i == 1) then
      ca = '****************'
    endif

    dddigout = ca
    return
  end function

type (dd_real) function ddhuge(a) 
  type (dd_real), intent(in) :: a
  ddhuge = dd_huge
end function ddhuge

type (dd_real) function dd_safe_huge(a)
  type (dd_real), intent(in) :: a
  dd_safe_huge = dd_real((/1.7976931080746007281d+308, &
                           9.97920154767359795037d+291/));
end function dd_safe_huge

type (dd_real) function ddtiny(a) 
  type (dd_real), intent(in) :: a
  ddtiny = dd_tiny
end function ddtiny

type (dd_real) function ddepsilon(a) 
  type (dd_real), intent(in) :: a
  ddepsilon = dd_eps
end function ddepsilon

integer function dd_radix(a)
  type (dd_real), intent(in) :: a
  dd_radix = 2
end function dd_radix

integer function dd_digits(a)
  type (dd_real), intent(in) :: a
  dd_digits = 104
end function dd_digits

integer function dd_max_expn(a)
  type (dd_real), intent(in) :: a
  dd_max_expn = 1023
end function dd_max_expn

integer function dd_min_expn(a)
  type (dd_real), intent(in) :: a
  dd_min_expn = -969
end function dd_min_expn

integer function dd_precision(a)
  type (dd_real), intent(in) :: a
  dd_precision = 31
end function dd_precision

integer function dd_range(a)
  type (dd_real), intent(in) :: a
  dd_range = 291
end function dd_range

type (dd_real) function dd_nan(a)
  type (dd_real), intent(in) :: a
  call f_dd_nan(dd_nan)
end function dd_nan

end module ddmodule


