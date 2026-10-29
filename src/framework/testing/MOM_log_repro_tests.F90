! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Unit tests for log_repro
module MOM_log_repro_tests

use, intrinsic :: ieee_arithmetic, only : ieee_value
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_signaling_nan
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

use MOM_error_handler, only : assert
use MOM_unit_testing, only : TestSuite
use MOM_intrinsic_functions, only : exp_repro, log_repro

implicit none ; private

public :: add_log_repro_tests

integer, parameter :: realquad = selected_real_kind(p=30, r=300)
  !< Potential real128 precision kind.  If unavailable, this will be negative.
integer, parameter :: realq = merge(realquad, kind(1.), realquad >= 0.)
  !< Placeholder real128 for declarations.  Unused if real128 is unavailable.
integer, parameter :: int_kind = selected_int_kind(18)
  !< Integer kind large enough to hold a binary64 bit pattern.

real, parameter :: rmold = 0.
  !< Mold value for default real [nondim]

! Mathematical constants
real, parameter :: ln2 = 0.69314718055994531
  !< ln(2) [nondim]
real, parameter :: half_ln2 = 0.34657359027997265
  !< ln(2)/2 [nondim]
real, parameter :: sqrt2 = 1.4142135623730950
  !< sqrt(2) [nondim]
real, parameter :: log_huge = 709.78271289338397
  !< log(huge()) [nondim]
real, parameter :: log_tiny = -708.39641853226408
  !< log(tiny()) [nondim]
real, parameter :: log_033 = -1.1086626245216111
  !< log(0.33) [nondim]
real, parameter :: log_2 = 0.69314718055994531
  !< log(2) [nondim]
real, parameter :: log_10 = 2.3025850929940457
  !< log(10) [nondim]
real, parameter :: log_smallest_subnormal = -744.44007192138122
  !< log(nearest(0, 1)) [nondim]
real, parameter :: log_one_plus_eps = 2.2204460492503128e-16
  !< log(1 + epsilon(1)) [nondim]
real, parameter :: log_one_minus_eps = -2.2204460492503136e-16
  !< log(1 - epsilon(1)) [nondim]
integer, parameter :: log_dump_npts = 100000
  !< Number of points in each optional log comparison dump.
integer(kind=int_kind), parameter :: log_dump_wide_start_bits = 4457293557087583675_int_kind
  !< Bit pattern for the lower end of the optional dump range, 1.e-10 (Z'3DDB7CDFD9D7BDBB').
integer(kind=int_kind), parameter :: log_dump_wide_end_bits = 4756540486875873280_int_kind
  !< Bit pattern for the upper end of the optional dump range, 1.e10 (Z'4202A05F20000000').

contains

!> Test that log_repro is elemental (works on arrays)
subroutine test_log_elemental
  real :: x(4), val(4), ref(4), err
    ! Inputs, results, references, and error [nondim]
  real, parameter :: tol = 1.e-14
    ! Error tolerance [nondim]
  integer :: i

  x(:) = [0.5, 1., 2., 10.]
  ref(:) = [-ln2, 0., ln2, log_10]

  val = log_repro(x)

  do i = 1, 4
    if (ref(i) /= 0.) then
      err = abs(val(i) - ref(i)) / abs(ref(i))
    else
      err = abs(val(i) - ref(i))
    endif
    call assert(err < tol, "log_repro elemental test failed")
  enddo
end subroutine test_log_elemental

!> Test that log_repro(1) = 0 exactly
subroutine test_log_one
  real :: val
    ! Result [nondim]

  val = log_repro(1.)

  print '(2x, "log_repro(1) = ", ES22.15)', val

  call assert(val == 0., "log_repro(1) should equal 0 exactly")
end subroutine test_log_one

!> Test log_repro at a general value (0.33)
subroutine test_log_general
  real :: x, val, ref, err
    ! Input, result, reference, and relative error [nondim]
  real, parameter :: tol = 1.e-12
    ! Relative tolerance [nondim]

  x = 0.33
  ref = log_033
  val = log_repro(x)
  err = abs(val - ref) / abs(ref)

  print '(2x, "log_repro(0.33) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(0.33) relative error exceeds tolerance")
end subroutine test_log_general

!> Test log_repro(2) against ln(2)
subroutine test_log_two
  real :: val, err
    ! Result and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  val = log_repro(2.)
  err = abs(val - log_2) / log_2

  print '(2x, "log_repro(2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(2) relative error exceeds tolerance")
end subroutine test_log_two

!> Test log_repro at range reduction boundary: 0.5
subroutine test_log_half
  real :: val, err
    ! Result and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  val = log_repro(0.5)
  err = abs(val + ln2) / ln2

  print '(2x, "log_repro(0.5) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(0.5) should be close to -ln(2)")
end subroutine test_log_half

!> Test log_repro at half the range reduction boundary: sqrt(2)
subroutine test_log_sqrt2
  real :: val, err
    ! Result and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  val = log_repro(sqrt2)
  err = abs(val - half_ln2) / half_ln2

  print '(2x, "log_repro(sqrt2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(sqrt(2)) should be close to ln(2)/2")
end subroutine test_log_sqrt2

!> Test log_repro at half the range reduction boundary: 1 / sqrt(2)
subroutine test_log_inv_sqrt2
  real :: val, err
    ! Result and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  val = log_repro(1. / sqrt2)
  err = abs(val + half_ln2) / half_ln2

  print '(2x, "log_repro(1/sqrt2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(1/sqrt(2)) should be close to -ln(2)/2")
end subroutine test_log_inv_sqrt2

!> Test log_repro(10) against reference value
subroutine test_log_ten
  real :: val, err
    ! Result and relative error [nondim]
  real, parameter :: tol = 1.e-12
    ! Relative tolerance [nondim]

  val = log_repro(10.)
  err = abs(val - log_10) / log_10

  print '(2x, "log_repro(10) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(10) relative error exceeds tolerance")
end subroutine test_log_ten

!> Test log_repro(+Inf) = +Inf
subroutine test_log_pos_inf
  real :: x, val, ref
    ! Input, result, and reference [nondim]

  x = ieee_value(0., ieee_positive_inf)
  ref = ieee_value(0., ieee_positive_inf)
  val = log_repro(x)

  print '(2x, "log_repro(+Inf) = ", ES22.15)', val

  call assert(val == ref, "log_repro(+Inf) should equal +Inf")
end subroutine test_log_pos_inf

!> Test log_repro(0) = -Inf
subroutine test_log_zero
  real :: val, ref
    ! Result and reference [nondim]

  ref = ieee_value(0., ieee_negative_inf)
  val = log_repro(0.)

  print '(2x, "log_repro(0) = ", ES22.15)', val

  call assert(val == ref, "log_repro(0) should equal -Inf")
end subroutine test_log_zero

!> Test log_repro(-0) = -Inf
subroutine test_log_neg_zero
  real :: val, ref
    ! Result and reference [nondim]

  ref = ieee_value(0., ieee_negative_inf)
  val = log_repro(-0.)

  print '(2x, "log_repro(-0) = ", ES22.15)', val

  call assert(val == ref, "log_repro(-0) should equal -Inf")
end subroutine test_log_neg_zero

!> Test log_repro(-Inf) = NaN
subroutine test_log_neg_inf
  real :: x, val
    ! Input and result [nondim]

  x = ieee_value(0., ieee_negative_inf)
  val = log_repro(x)

  print '(2x, "log_repro(-Inf) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-Inf) should be NaN")
end subroutine test_log_neg_inf

!> Test log_repro(negative) = NaN
subroutine test_log_negative
  real :: val
    ! Result [nondim]

  val = log_repro(-1.)

  print '(2x, "log_repro(-1) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-1) should be NaN")
end subroutine test_log_negative

!> Test log_repro(NaN) = NaN
subroutine test_log_nan
  real :: x, val
    ! Input and result [nondim]

  x = ieee_value(0., ieee_quiet_nan)
  val = log_repro(x)

  print '(2x, "log_repro(NaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(NaN) should be NaN")
end subroutine test_log_nan

!> Test log_repro(-NaN) = NaN
subroutine test_log_neg_nan
  real :: x, val
    ! Input and result [nondim]

  x = -ieee_value(0., ieee_quiet_nan)
  val = log_repro(x)

  print '(2x, "log_repro(-NaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-NaN) should be NaN")
end subroutine test_log_neg_nan

!> Test log_repro(sNaN) = NaN
subroutine test_log_signaling_nan
  real :: x, val
    ! Input and result [nondim]

  x = ieee_value(0., ieee_signaling_nan)
  val = log_repro(x)

  print '(2x, "log_repro(sNaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(sNaN) should be NaN")
end subroutine test_log_signaling_nan

!> Test log_repro(-sNaN) = NaN
subroutine test_log_neg_signaling_nan
  real :: x, val
    ! Input and result [nondim]

  x = -ieee_value(0., ieee_signaling_nan)
  val = log_repro(x)

  print '(2x, "log_repro(-sNaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-sNaN) should be NaN")
end subroutine test_log_neg_signaling_nan

!> Test log_repro at largest representable float boundary
subroutine test_log_largest_float
  real :: val, ref, err
    ! Result, reference, and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  ref = log_huge
  val = log_repro(huge(1.))
  err = abs(val - ref) / ref

  print '(2x, "log_repro(huge) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(huge) relative error exceeds tolerance")
end subroutine test_log_largest_float

!> Test log_repro at smallest normal float boundary
subroutine test_log_smallest_float
  real :: val, ref, err
    ! Result, reference, and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  ref = log_tiny
  val = log_repro(tiny(1.))
  err = abs(val - ref) / abs(ref)

  print '(2x, "log_repro(tiny) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(tiny) relative error exceeds tolerance")
end subroutine test_log_smallest_float

!> Test log_repro at smallest positive subnormal float boundary
subroutine test_log_smallest_subnormal
  real :: val, ref, err
    ! Result, reference, and relative error [nondim]
  real, parameter :: tol = 1.e-15
    ! Relative tolerance [nondim]

  ref = log_smallest_subnormal
  val = log_repro(nearest(0., 1.))
  err = abs(val - ref) / abs(ref)

  print '(2x, "log_repro(nearest(0,1)) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(nearest(0,1)) relative error exceeds tolerance")
end subroutine test_log_smallest_subnormal

!> Test log_repro near 1 from above
subroutine test_log_near_one_above
  real :: x, val, ref, err
    ! Input, result, reference, and absolute error [nondim]
  real, parameter :: tol = 1.e-30
    ! Absolute tolerance [nondim]

  x = 1. + epsilon(1.)
  ref = log_one_plus_eps
  val = log_repro(x)
  err = abs(val - ref)

  print '(2x, "log_repro(1+eps) = ", ES22.15, ", abs err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(1+eps) absolute error exceeds tolerance")
end subroutine test_log_near_one_above

!> Test log_repro near 1 from below
subroutine test_log_near_one_below
  real :: x, val, ref, err
    ! Input, result, reference, and absolute error [nondim]
  real, parameter :: tol = 1.e-30
    ! Absolute tolerance [nondim]

  x = 1. - epsilon(1.)
  ref = log_one_minus_eps
  val = log_repro(x)
  err = abs(val - ref)

  print '(2x, "log_repro(1-eps) = ", ES22.15, ", abs err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(1-eps) absolute error exceeds tolerance")
end subroutine test_log_near_one_below

!> Test a near-one input that exposed a one-bit GCC/NVHPC difference.
subroutine test_log_near_one_sensitive_bit
  integer(kind=int_kind), parameter :: x_bits = 4607434661109602316_int_kind
    !< Bit pattern Z'3FF0E569BD42980C', a sensitive near-one input.
  integer(kind=int_kind), parameter :: result_bits_1 = 4588014630356761006_int_kind
    !< Observed result bit pattern Z'3FABE6FF7AB5B9AE'.
  integer(kind=int_kind), parameter :: result_bits_2 = 4588014630356761007_int_kind
    !< Observed result bit pattern Z'3FABE6FF7AB5B9AF'.

  real :: x, val, ref, ulp_err
    ! Input, result, reference, and ULP error [nondim]
  integer(kind=int_kind) :: val_bits
    !< Bit pattern of the log_repro() result.

  x = transfer(x_bits, 0.)
  val = log_repro(x)
  ref = real(log(real(x, realq)))
  val_bits = transfer(val, 0_int_kind)
  ulp_err = abs(val - ref) / abs(spacing(ref))

  print '(2x,"log_repro sensitive near-1 x_bits = ",a)', real_bits_hex(x)
  print '(2x,"log_repro sensitive near-1 result = ",ES22.15,", bits = ",a,", ULP err = ",f12.10)', &
      val, real_bits_hex(val), ulp_err

  call assert((val_bits == result_bits_1) .or. (val_bits == result_bits_2), &
      "log_repro sensitive near-one result is not one of the observed adjacent bit patterns")
  call assert(ulp_err <= 1., "log_repro sensitive near-one ULP error exceeds tolerance")
end subroutine test_log_near_one_sensitive_bit

!> Print an informational accuracy sweep for log_repro near 1.
!!
!! The broad bit-pattern sweep does not sample values close to 1 densely.
!! This diagnostic sweep extends beyond both the original and symmetric
!! near-one cutoffs, comparing both paths with the same quad-precision reference.
!! Integer bit coordinates keep the input grid independent of FMA contraction.
subroutine test_log_near_one_sweep
  integer, parameter :: npts = 100000
  real, parameter :: umin = -0.075
    ! Lower sweep offset [nondim]
  real, parameter :: umax = 0.075
    ! Upper sweep offset [nondim]

  real :: x, val, val_log
    ! Input, log_repro result, and log result [nondim]
  real(kind=realq) :: ref
    ! Quad reference [nondim]
  real(kind=realq) :: abs_err, rel_err
    ! Absolute and relative log_repro errors [nondim]
  real(kind=realq) :: abs_err_log, rel_err_log
    ! Absolute and relative log errors [nondim]
  real :: ulp_err, ulp_err_log
    ! ULP errors for log_repro and log [nondim]
  real :: max_abs_err, max_rel_err, max_ulp
    ! Maximum absolute, relative, and ULP log_repro errors [nondim]
  real :: max_abs_err_log, max_rel_err_log, max_ulp_log
    ! Maximum absolute, relative, and ULP log errors [nondim]
  real :: x_max_abs, x_max_rel, x_max_ulp
    ! Inputs at maximum log_repro errors [nondim]
  real :: x_max_abs_log, x_max_rel_log, x_max_ulp_log
    ! Inputs at maximum log errors [nondim]
  real :: ulp_val
    ! Reference spacing [nondim]
  integer(kind=int_kind) :: start_bits, end_bits, span, step, rem, offset, x_bits
    !< Input bit coordinates and integer-grid spacing
  integer(kind=int_kind) :: denom
    !< Number of intervals in the input grid

  integer :: i, n_half_ulp, n_half_ulp_log

  max_abs_err = 0.
  max_rel_err = 0.
  max_ulp = 0.
  max_abs_err_log = 0.
  max_rel_err_log = 0.
  max_ulp_log = 0.
  x_max_abs = 1.
  x_max_rel = 1.
  x_max_ulp = 1.
  x_max_abs_log = 1.
  x_max_rel_log = 1.
  x_max_ulp_log = 1.
  n_half_ulp = 0
  n_half_ulp_log = 0

  start_bits = transfer(1. + umin, 0_int_kind)
  end_bits = transfer(1. + umax, 0_int_kind)
  span = end_bits - start_bits
  denom = int(npts - 1, int_kind)
  step = span / denom
  rem = modulo(span, denom)
  do i = 1, npts
    offset = int(i - 1, int_kind)
    x_bits = start_bits + offset * step + (offset * rem) / denom
    x = transfer(x_bits, 0.)
    val = log_repro(x)
    val_log = log(x)
    ref = log(real(x, realq))
    ulp_val = abs(spacing(real(ref, kind(rmold))))

    abs_err = abs(real(val, realq) - ref)
    if (abs_err > max_abs_err) then
      max_abs_err = abs_err
      x_max_abs = x
    endif

    abs_err_log = abs(real(val_log, realq) - ref)
    if (abs_err_log > max_abs_err_log) then
      max_abs_err_log = abs_err_log
      x_max_abs_log = x
    endif

    if (ref /= 0.) then
      rel_err = abs_err / abs(ref)
      if (rel_err > max_rel_err) then
        max_rel_err = rel_err
        x_max_rel = x
      endif

      rel_err_log = abs_err_log / abs(ref)
      if (rel_err_log > max_rel_err_log) then
        max_rel_err_log = rel_err_log
        x_max_rel_log = x
      endif
    endif

    ulp_err = real(abs_err, kind(rmold)) / ulp_val
    if (ulp_err > max_ulp) then
      max_ulp = ulp_err
      x_max_ulp = x
    endif
    if (ulp_err > 0.5) n_half_ulp = n_half_ulp + 1

    ulp_err_log = real(abs_err_log, kind(rmold)) / ulp_val
    if (ulp_err_log > max_ulp_log) then
      max_ulp_log = ulp_err_log
      x_max_ulp_log = x
    endif
    if (ulp_err_log > 0.5) n_half_ulp_log = n_half_ulp_log + 1

  enddo

  print '(2x,"Tested ", i0, " values across near-1 interval [", ES12.5, ",", ES12.5, "]")', &
      npts, 1. + umin, 1. + umax
  print '(2x,"log_repro max abs err:", t34, ES12.5, " at input = ", ES22.15)', &
      max_abs_err, x_max_abs
  print '(2x,"log_repro max rel err:", t34, ES12.5, " at input = ", ES22.15)', &
      max_rel_err, x_max_rel
  print '(2x,"log_repro max ULP err:", t34, f12.10, " at input = ", ES22.15)', &
      max_ulp, x_max_ulp
  print '(2x,"log_repro points > 0.5 ULP:", t34, i12)', n_half_ulp
  print '(2x,"log() max abs err:", t34, ES12.5, " at input = ", ES22.15)', &
      max_abs_err_log, x_max_abs_log
  print '(2x,"log() max rel err:", t34, ES12.5, " at input = ", ES22.15)', &
      max_rel_err_log, x_max_rel_log
  print '(2x,"log() max ULP err:", t34, f12.10, " at input = ", ES22.15)', &
      max_ulp_log, x_max_ulp_log
  print '(2x,"log() points > 0.5 ULP:", t34, i12)', n_half_ulp_log
end subroutine test_log_near_one_sweep


!> Check the immediate neighbors of the original and candidate near-one cutoffs.
subroutine test_log_near_one_boundaries
  real, parameter :: boundaries(3) = [1. - 1./16., 1. + 1./16., 1. + 265./4096.]
    ! Original and candidate cutoffs [nondim]
  real :: x(3*size(boundaries)), val(3*size(boundaries))
    ! Boundary inputs and results [nondim]
  real(kind=realq) :: ref(3*size(boundaries))
    ! Quad references [nondim]
  integer :: i, j

  do i = 1, size(boundaries)
    j = 3 * i - 2
    x(j) = nearest(boundaries(i), -1.)
    x(j + 1) = boundaries(i)
    x(j + 2) = nearest(boundaries(i), 1.)
  enddo
  val(:) = log_repro(x(:))
  ref(:) = log(real(x(:), realq))
  call check_ulp_accuracy(x, val, ref, max_ulp_tol=1.)
end subroutine test_log_near_one_boundaries


!> Test ULP accuracy over a wide range of positive values.
!!
!! log_repro should be within a few ULP of the true value.
subroutine test_log_ulp_accuracy
  integer, parameter :: npts = 100000
  real, parameter :: xmin = 1.e-10
    ! Lower input bound [nondim]
  real, parameter :: xmax = 1.e10
    ! Upper input bound [nondim]

  real :: x(npts)
    ! Input grid [nondim]
  integer(kind=int_kind) :: x_bits, span, step, rem, offset
    !< Input bit pattern and integer-grid spacing values.
  integer(kind=int_kind) :: denom
    !< Denominator for the uniformly spaced integer grid.

  real :: val(npts), val_vec(npts)
    ! Scalar and vector log_repro results [nondim]
  real :: val_log(npts), val_log_vec(npts)
    ! Scalar and vector log results [nondim]
  real(kind=realq) :: val_quad(npts), val_quad_vec(npts)
    ! Scalar and vector quad references [nondim]

  integer :: i

  ! Generate reproducible positive test points with uniformly spaced binary64
  ! bit coordinates.  This gives broad scale-aware coverage without relying on
  ! exp() or exp_repro() to define the test inputs.
  span = log_dump_wide_end_bits - log_dump_wide_start_bits
  denom = int(npts - 1, int_kind)
  step = span / denom
  rem = modulo(span, denom)

  do i = 1, npts
    offset = int(i - 1, int_kind)
    x_bits = log_dump_wide_start_bits + offset * step + (offset * rem) / denom
    x(i) = transfer(x_bits, 0.)
  enddo

  ! Several libraries have scalar and vector implementations, chosen at the
  ! discretion of the compiler.  The following attempts to test each case.

  ! Scalar evaluations
  do i = 1, npts
    val_log(i) = log(x(i))
    val(i) = log_repro(x(i))
    val_quad(i) = log(real(x(i), realq))

    ! Impossible branch to prevent vectorization
    if (val(i) > huge(val(i))) exit
  enddo

  ! Vector-favorable evaluation
  val_log_vec = log(x)
  val_vec = log_repro(x)
  val_quad_vec = log(real(x, realq))

  ! Assert that log_repro() is within 8 ULP and compare with intrinsic log().
  print '(1x,a)', '=== scalar log_repro() and log() accuracy'
  call check_log_accuracy_comparison(x, val, val_log, val_quad, xmin, xmax, max_ulp_tol=8.)

  ! We expect scalar and vector implementations to agree.
  call assert(all(val == val_vec), 'Scalar and vector log_repro() do not agree')
  print '(1x,a)', '=== vector log_repro() and log() accuracy'
  call check_log_accuracy_comparison(x, val_vec, val_log_vec, val_quad_vec, xmin, xmax, max_ulp_tol=8.)
end subroutine test_log_ulp_accuracy


!> Test bitwise agreement between scalar log_repro() and a do-concurrent loop.
!! When compiled with NVIDIA standard parallelism offload, this exercises
!! log_repro() in a GPU kernel while keeping reference generation and bit
!! comparisons on the host.
subroutine test_log_do_concurrent_bitwise
  integer, parameter :: npts = 10000

  real :: x(npts), val_scalar(npts), val_concurrent(npts)
    ! Inputs, scalar results, and do-concurrent results [nondim]
  integer(kind=int_kind) :: x_bits, start_bits, end_bits, span, step, rem, offset
    !< Input bit pattern and integer-grid spacing values.
  integer(kind=int_kind) :: bits_scalar, bits_concurrent
    !< Binary64 bit patterns for scalar and do-concurrent results.
  integer(kind=int_kind) :: denom
    !< Denominator for the uniformly spaced integer grid.
  integer :: i

  start_bits = log_dump_wide_start_bits
  end_bits = log_dump_wide_end_bits
  span = end_bits - start_bits
  denom = int(npts - 1, int_kind)
  step = span / denom
  rem = modulo(span, denom)

  do i = 1, npts
    offset = int(i - 1, int_kind)
    x_bits = start_bits + offset * step + (offset * rem) / denom
    x(i) = transfer(x_bits, 0.)
  enddo

  ! Add cases that stress special handling without relying on quad precision or
  ! intrinsic log() inside the offloaded loop.
  x(1) = 0.
  x(2) = tiny(0.)
  x(3) = nearest(0., 1.)
  x(4) = huge(0.)
  x(5) = 1. - epsilon(1.)
  x(6) = 1.
  x(7) = 1. + epsilon(1.)
  x(8) = ieee_value(0., ieee_positive_inf)

  do i = 1, npts
    val_scalar(i) = log_repro(x(i))
  enddo

  call eval_log_repro_do_concurrent(x, val_concurrent)

  do i = 1, npts
    bits_scalar = transfer(val_scalar(i), 0_int_kind)
    bits_concurrent = transfer(val_concurrent(i), 0_int_kind)

    if (bits_scalar /= bits_concurrent) then
      print '(2x,"log_repro do-concurrent mismatch at index ",i0)', i
      print '(4x,"x_bits             = ",a)', real_bits_hex(x(i))
      print '(4x,"scalar bits        = ",a)', real_bits_hex(val_scalar(i))
      print '(4x,"do-concurrent bits = ",a)', real_bits_hex(val_concurrent(i))
      call assert(.false., "do-concurrent log_repro() differs from scalar log_repro()")
    endif
  enddo

  print '(2x,"do-concurrent log_repro() matched scalar results for ",i0," inputs")', npts
end subroutine test_log_do_concurrent_bitwise


!> Evaluate log_repro() over an array using a do-concurrent loop.
subroutine eval_log_repro_do_concurrent(x, val)
  real, intent(in) :: x(:)
    !< Inputs [nondim]
  real, intent(out) :: val(:)
    !< Results [nondim]

  integer :: i

  !$omp target data map(to: x) map(tofrom: val)
  do concurrent (i=1:size(x))
    val(i) = log_repro(x(i))
  enddo
  !$omp end target data
end subroutine eval_log_repro_do_concurrent


!> Compare scalar log_repro() and intrinsic log() accuracy against a real128 reference.
subroutine check_log_accuracy_comparison(x, val, val_log, ref, xmin, xmax, max_ulp_tol)
  real, intent(in) :: x(:)
    !< Input grid [nondim]
  real, intent(in) :: val(:)
    !< log_repro results [nondim]
  real, intent(in) :: val_log(:)
    !< Intrinsic log results [nondim]
  real(kind=realq), intent(in) :: ref(:)
    !< Quad reference [nondim]
  real, intent(in) :: xmin
    !< Lower input bound [nondim]
  real, intent(in) :: xmax
    !< Upper input bound [nondim]
  real, optional, intent(in) :: max_ulp_tol
    !< ULP tolerance [nondim]

  real(kind=realq) :: err, err_log, rel_err, rel_err_log
    ! Absolute and relative errors for log_repro and log [nondim]
  real :: max_abs_err, max_rel_err, max_ulp
    ! Maximum absolute, relative, and ULP log_repro errors [nondim]
  real :: max_abs_err_log, max_rel_err_log, max_ulp_log
    ! Maximum absolute, relative, and ULP log errors [nondim]
  real :: x_max_abs, x_max_rel, x_max_ulp
    ! Inputs at maximum log_repro errors [nondim]
  real :: x_max_abs_log, x_max_rel_log, x_max_ulp_log
    ! Inputs at maximum log errors [nondim]
  real :: ulp_val, ulp_err, ulp_err_log
    ! Reference spacing and ULP errors for log_repro and log [nondim]

  integer :: count_half_ulp, count_half_ulp_log
    !< Number of points with errors above 0.5 ULP
  integer :: i, npts

  npts = size(x)

  max_abs_err = 0.
  max_rel_err = 0.
  max_ulp = 0.
  max_abs_err_log = 0.
  max_rel_err_log = 0.
  max_ulp_log = 0.
  x_max_abs = x(1)
  x_max_rel = x(1)
  x_max_ulp = x(1)
  x_max_abs_log = x(1)
  x_max_rel_log = x(1)
  x_max_ulp_log = x(1)
  count_half_ulp = 0
  count_half_ulp_log = 0

  do i = 1, npts
    ulp_val = abs(spacing(real(ref(i), kind(rmold))))

    err = abs(real(val(i), realq) - ref(i))
    if (err > max_abs_err) then
      max_abs_err = err
      x_max_abs = x(i)
    endif

    err_log = abs(real(val_log(i), realq) - ref(i))
    if (err_log > max_abs_err_log) then
      max_abs_err_log = err_log
      x_max_abs_log = x(i)
    endif

    if (ref(i) /= 0.) then
      rel_err = err / abs(ref(i))
      if (rel_err > max_rel_err) then
        max_rel_err = rel_err
        x_max_rel = x(i)
      endif

      rel_err_log = err_log / abs(ref(i))
      if (rel_err_log > max_rel_err_log) then
        max_rel_err_log = rel_err_log
        x_max_rel_log = x(i)
      endif
    endif

    ulp_err = real(err, kind(rmold)) / ulp_val
    if (ulp_err > max_ulp) then
      max_ulp = ulp_err
      x_max_ulp = x(i)
    endif
    if (ulp_err > 0.5) count_half_ulp = count_half_ulp + 1

    ulp_err_log = real(err_log, kind(rmold)) / ulp_val
    if (ulp_err_log > max_ulp_log) then
      max_ulp_log = ulp_err_log
      x_max_ulp_log = x(i)
    endif
    if (ulp_err_log > 0.5) count_half_ulp_log = count_half_ulp_log + 1

  enddo

  print '(2x,"Tested ", i0, " values across x interval [", ES12.5, ",", ES12.5, "]")', &
      npts, xmin, xmax
  print '(2x,"log_repro max abs err:", t34, ES12.5, " at input = ", ES14.5E3)', &
      max_abs_err, x_max_abs
  print '(2x,"log_repro max rel err:", t34, ES12.5, " at input = ", ES14.5E3)', &
      max_rel_err, x_max_rel
  print '(2x,"log_repro max ULP err:", t34, f12.10, " at input = ", ES14.5E3)', &
      max_ulp, x_max_ulp
  print '(2x,"log_repro points > 0.5 ULP:", t34, i12)', count_half_ulp
  print '(2x,"log() max abs err:", t34, ES12.5, " at input = ", ES14.5E3)', &
      max_abs_err_log, x_max_abs_log
  print '(2x,"log() max rel err:", t34, ES12.5, " at input = ", ES14.5E3)', &
      max_rel_err_log, x_max_rel_log
  print '(2x,"log() max ULP err:", t34, f12.10, " at input = ", ES14.5E3)', &
      max_ulp_log, x_max_ulp_log
  print '(2x,"log() points > 0.5 ULP:", t34, i12)', count_half_ulp_log

  if (present(max_ulp_tol)) then
    call assert(max_ulp < max_ulp_tol, "Max ULP error exceeds tolerance")
  endif
end subroutine check_log_accuracy_comparison

!> Compute the function accuracy relative to a real128-precision reference.
!! Absolute, relative, and ULP error is computed, as well as the number of
!! points above 0.5 and 1 ULP.  An error is raised if any point exceeds the
!! optionally prescribed maximum ULP tolerance.
subroutine check_ulp_accuracy(x, val, ref, max_ulp_tol)
  real, intent(in) :: x(:)
    !< Input grid [nondim]
  real, intent(in) :: val(:)
    !< Results [nondim]
  real(kind=realq), intent(in) :: ref(:)
    !< Quad reference [nondim]
  real, optional, intent(in) :: max_ulp_tol
    !< ULP tolerance [nondim]

  real(kind=realq) :: err, rel_err, ulp_err
    ! Absolute, relative, and ULP errors [nondim]

  real :: max_abs_err, max_rel_err, max_ulp
    ! Maximum absolute, relative, and ULP errors [nondim]
  real :: sum_abs_err, sum_rel_err, sum_sq_err
    ! Sums of absolute, relative, and squared errors [nondim]
  real :: x_max_abs, x_max_rel, x_max_ulp
    ! Inputs at maximum absolute, relative, and ULP errors [nondim]
  real :: ulp_val
    ! Reference spacing [nondim]

  integer :: count_exact, count_half_ulp, count_one_ulp
  integer :: i, npts

  npts = size(x)

  ! Initialize statistics
  max_abs_err = 0.
  max_rel_err = 0.
  max_ulp = 0.
  sum_abs_err = 0.
  sum_rel_err = 0.
  sum_sq_err = 0.
  count_exact = 0
  count_half_ulp = 0
  count_one_ulp = 0

  do i = 1, npts
    ! Absolute error
    err = abs(real(val(i), realq) - ref(i))

    sum_abs_err = sum_abs_err + err
    sum_sq_err = sum_sq_err + err * err

    ! Update max error
    if (err > max_abs_err) then
      max_abs_err = err
      x_max_abs = x(i)
    endif

    ! Relative error

    if (ref(i) /= 0.) then
      rel_err = err / abs(ref(i))

      sum_rel_err = sum_rel_err + rel_err

      if (rel_err > max_rel_err) then
        max_rel_err = rel_err
        x_max_rel = x(i)
      endif
    endif

    ! Report error relative to ULP

    ulp_val = spacing(real(ref(i), kind(rmold)))
    ulp_err = abs(val(i) - ref(i)) / ulp_val

    if (ulp_err > max_ulp) then
      max_ulp = ulp_err
      x_max_ulp = x(i)
    endif

    if (ulp_err < 0.5) count_exact = count_exact + 1
    if (ulp_err >= 0.5) count_half_ulp = count_half_ulp + 1
    if (ulp_err >= 1.0) count_one_ulp = count_one_ulp + 1
  enddo

  print '(2x,"Tested ", i0, " points")', npts

  ! NOTE: Floats use t25 assuming a positive sign as blank
  print '(2x,"max abs err:", t25, ES12.5, " at input = ", ES14.5E3)', &
      max_abs_err, x_max_abs
  print '(2x,"max rel err:", t25, ES12.5, " at input = ", ES14.5E3)', &
      max_rel_err, x_max_rel
  print '(2x,"max ULP err (vs quad):", t26, f12.10, " at input = ", ES14.5E3)', &
      max_ulp, x_max_ulp
  print '(2x,"mean abs err:", t25, ES12.5)', sum_abs_err / npts
  print '(2x,"mean rel err:", t25, ES12.5)', sum_rel_err / npts
  print '(2x,"RMS err:", t25, ES12.5)', sqrt(sum_sq_err / npts)

  print '(2x,"correct (<0.5 ULP):", t25, i0, 1x, "(", f6.2, "%)")', &
      count_exact, 100. * count_exact / npts
  print '(2x,"above 0.5 ULP:", t26, i0, 1x, "(", f6.2, "%)")', &
      count_half_ulp, 100. * count_half_ulp / npts
  print '(2x,"above 1 ULP:", t26, i0, 1x, "(", f6.2, "%)")', &
      count_one_ulp, 100. * count_one_ulp / npts

  ! exp_repro should be within 2 ULP of the quad-precision reference.
  if (present(max_ulp_tol)) then
    call assert(max_ulp < max_ulp_tol, "Max ULP error exceeds tolerance")
  endif
end subroutine check_ulp_accuracy

!> Test the logarithmic property: log(x*y) = log(x) + log(y)
!!
!! This property test verifies that log_repro produces consistent results
!! by checking the logarithm of a product against the sum of logarithms,
!! within floating-point tolerance.
subroutine test_log_product_property
  integer, parameter :: npts = 1000
  real, parameter :: tol = 1.e-13
    ! Error tolerance [nondim]

  real :: r1, r2
    ! Input logarithms [nondim]
  real :: x1, x2
    ! Test inputs [nondim]
  real :: log_product, log_sum
    ! Product and sum-form log_repro results [nondim]
  real :: log_product_std, log_sum_std
    ! Product and sum-form log results [nondim]
  real :: err, err_std, max_err, max_err_std
    ! Errors and maxima for log_repro and log [nondim]
  real :: r1_max, r2_max, r1_max_std, r2_max_std
    ! Input logarithms at maximum errors [nondim]
  integer :: i
  real :: seed
    ! Pseudo-random seed [nondim]

  max_err = 0.
  max_err_std = 0.

  ! Use a simple deterministic sequence for reproducibility.
  seed = 0.246813579

  do i = 1, npts
    ! Generate pseudo-random logarithms in a range where the product is finite.
    ! Keep r1, r2 in [-5, 5] so x1*x2 is in exp([-10, 10]).
    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r1 = (seed / 2.**31) * 10. - 5.

    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r2 = (seed / 2.**31) * 10. - 5.

    x1 = exp_repro(r1)
    x2 = exp_repro(r2)

    log_product = log_repro(x1 * x2)
    log_sum = log_repro(x1) + log_repro(x2)
    log_product_std = log(x1 * x2)
    log_sum_std = log(x1) + log(x2)

    if (log_product /= 0.) then
      err = abs(log_sum - log_product) / abs(log_product)
    else
      err = abs(log_sum - log_product)
    endif

    if (err > max_err) then
      max_err = err
      r1_max = r1
      r2_max = r2
    endif

    if (log_product_std /= 0.) then
      err_std = abs(log_sum_std - log_product_std) / abs(log_product_std)
    else
      err_std = abs(log_sum_std - log_product_std)
    endif

    if (err_std > max_err_std) then
      max_err_std = err_std
      r1_max_std = r1
      r2_max_std = r2
    endif
  enddo

  print '("Tested ", i0, " random positive (x1, x2) pairs")', npts
  print '("log_repro max rel err in log_repro(x1*x2) vs log_repro(x1)+log_repro(x2):", t82, ES12.5)', &
      max_err
  print '("  at log(x1) = ", f10.6, ", log(x2) = ", f10.6)', r1_max, r2_max
  print '("log() max rel err in log(x1*x2) vs log(x1)+log(x2):", t82, ES12.5)', max_err_std
  print '("  at log(x1) = ", f10.6, ", log(x2) = ", f10.6)', r1_max_std, r2_max_std

  call assert(max_err < tol, "log_repro product property test failed")
end subroutine test_log_product_property

!> Test the logarithmic property: log(1/x) = -log(x)
!!
!! This property test verifies that log_repro(1/x) equals the negation of
!! log_repro(x), within floating-point tolerance.
subroutine test_log_reciprocal_property
  integer, parameter :: npts = 1000
  real, parameter :: tol = 1.e-13
    ! Error tolerance [nondim]

  real :: r
    ! Input logarithm [nondim]
  real :: x
    ! Test input [nondim]
  real :: log_inv_x, neg_log_x
    ! Reciprocal-input and negated log_repro results [nondim]
  real :: log_inv_x_std, neg_log_x_std
    ! Reciprocal-input and negated log results [nondim]
  real :: err, err_std, max_err, max_err_std
    ! Errors and maxima for log_repro and log [nondim]
  real :: r_max, r_max_std
    ! Input logarithms at maximum errors [nondim]
  integer :: i
  real :: seed
    ! Pseudo-random seed [nondim]

  max_err = 0.
  max_err_std = 0.

  ! Use a simple deterministic sequence for reproducibility.
  seed = 0.135792468

  do i = 1, npts
    ! Generate pseudo-random logarithms that avoid overflow and underflow.
    ! Keep r in [-300, 300] so x and 1/x are both representable.
    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r = (seed / 2.**31) * 600. - 300.

    x = exp_repro(r)

    log_inv_x = log_repro(1. / x)
    neg_log_x = -log_repro(x)
    log_inv_x_std = log(1. / x)
    neg_log_x_std = -log(x)

    if (neg_log_x /= 0.) then
      err = abs(log_inv_x - neg_log_x) / abs(neg_log_x)
    else
      err = abs(log_inv_x - neg_log_x)
    endif

    if (err > max_err) then
      max_err = err
      r_max = r
    endif

    if (neg_log_x_std /= 0.) then
      err_std = abs(log_inv_x_std - neg_log_x_std) / abs(neg_log_x_std)
    else
      err_std = abs(log_inv_x_std - neg_log_x_std)
    endif

    if (err_std > max_err_std) then
      max_err_std = err_std
      r_max_std = r
    endif
  enddo

  print '("Tested ", i0, " random positive x values")', npts
  print '("log_repro max rel err in log_repro(1/x) vs -log_repro(x):", t72, ES12.5)', max_err
  print '("  at log(x) = ", f10.6)', r_max
  print '("log() max rel err in log(1/x) vs -log(x):", t72, ES12.5)', max_err_std
  print '("  at log(x) = ", f10.6)', r_max_std

  call assert(max_err < tol, "log_repro reciprocal property test failed")
end subroutine test_log_reciprocal_property


!> Optionally print deterministic log comparison values for platform comparisons.
!!
!! Set MOM_DUMP_LOG in the environment to print two CSV tables: one with
!! points equally spaced in log(x) over [-700, 700], and one with points equally
!! spaced in x over the near-one special-case interval.
subroutine test_log_value_dump
  character(len=8) :: concurrent_mode
  character(len=1) :: dump_enabled
  integer :: status
  logical :: use_do_concurrent

  call get_environment_variable("MOM_DUMP_LOG", dump_enabled, status=status)
  if (status /= 0) return

  use_do_concurrent = .false.
  call get_environment_variable("MOM_DUMP_LOG_DO_CONCURRENT", concurrent_mode, status=status)
  if (status == 0) then
    if ((concurrent_mode(1:1) == "1") .or. (concurrent_mode(1:1) == "T") .or. &
        (concurrent_mode(1:1) == "t") .or. (concurrent_mode(1:1) == "Y") .or. &
        (concurrent_mode(1:1) == "y")) use_do_concurrent = .true.
  endif

  print '(1x,a)', '=== optional log comparison dump'
  print '(a)', 'case,index,x_bits,log_repro_bits,log_bits,quad_ref_bits,x,log_repro,log,quad_ref'
  call dump_log_bit_grid("wide", log_dump_npts, log_dump_wide_start_bits, log_dump_wide_end_bits, use_do_concurrent)
  call dump_log_bit_grid("near1", log_dump_npts, &
      transfer(1. - (1. / 16.), 0_int_kind), transfer(1. + 0.064697265625, 0_int_kind), &
      use_do_concurrent)
end subroutine test_log_value_dump


!> Print one deterministic binary64 bit-grid of log comparison values.
subroutine dump_log_bit_grid(grid_name, npts, start_bits, end_bits, use_do_concurrent)
  character(len=*), intent(in) :: grid_name
    !< Label for the dumped grid.
  integer, intent(in) :: npts
    !< Number of points to print.
  integer(kind=int_kind), intent(in) :: start_bits
    !< Binary64 bit pattern for the first input value.
  integer(kind=int_kind), intent(in) :: end_bits
    !< Binary64 bit pattern for the last input value.
  logical, intent(in) :: use_do_concurrent
    !< If true, evaluate log_repro() using a do-concurrent loop.

  real :: x_vals(npts), val_repro_vals(npts), val_log_vals(npts), ref_vals(npts)
    ! Inputs, log_repro results, log results, and rounded references [nondim]
  real(kind=realq) :: ref_quad_vals(npts)
    ! Quad references [nondim]
  integer(kind=int_kind) :: x_bits_vals(npts)
    !< Input bit patterns.
  integer(kind=int_kind) :: span, step, rem
    !< Integer-grid spacing values.
  integer(kind=int_kind) :: denom
    !< Denominator for the uniformly spaced integer grid.
  integer :: i

  span = end_bits - start_bits
  denom = int(npts - 1, int_kind)
  step = span / denom
  rem = modulo(span, denom)

  if (use_do_concurrent) then
    !$omp target data map(from: x_bits_vals)
    do concurrent (i=1:npts)
      x_bits_vals(i) = start_bits + int(i - 1, int_kind) * step + &
          (int(i - 1, int_kind) * rem) / denom
    enddo
    !$omp end target data

    !$omp target data map(to: x_bits_vals) map(from: x_vals)
    do concurrent (i=1:npts)
      ! Decode positive normal binary64 inputs exactly without device TRANSFER.
      x_vals(i) = scale(1. + scale(real(ibits(x_bits_vals(i), 0, 52)), -52), &
          int(ibits(x_bits_vals(i), 52, 11)) - 1023)
    enddo
    !$omp end target data

    ! Check device input generation independently of the logarithm calls.
    call check_log_dump_grid(x_vals, x_bits_vals, start_bits, end_bits)
    !$omp target data map(to: x_vals) map(from: val_repro_vals)
    do concurrent (i=1:npts)
      val_repro_vals(i) = log_repro(x_vals(i))
    enddo
    !$omp end target data

    !call check_log_dump_grid(x_vals, x_bits_vals, start_bits, end_bits)
  else
    do i = 1, npts
      x_bits_vals(i) = start_bits + int(i - 1, int_kind) * step + &
          (int(i - 1, int_kind) * rem) / denom
      x_vals(i) = scale(1. + scale(real(ibits(x_bits_vals(i), 0, 52)), -52), &
          int(ibits(x_bits_vals(i), 52, 11)) - 1023)
      val_repro_vals(i) = log_repro(x_vals(i))
    enddo
  endif

  call check_log_dump_grid(x_vals, x_bits_vals, start_bits, end_bits)

  do i = 1, npts
    val_log_vals(i) = log(x_vals(i))
    ref_quad_vals(i) = log(real(x_vals(i), realq))
    ref_vals(i) = real(ref_quad_vals(i))

    print '(a,",",i0,",",a,",",a,",",a,",",a,3(",",ES26.17E3),",",ES46.37E3)', &
        trim(grid_name), i, real_bits_hex(x_vals(i)), real_bits_hex(val_repro_vals(i)), &
        real_bits_hex(val_log_vals(i)), real_bits_hex(ref_vals(i)), &
        x_vals(i), val_repro_vals(i), val_log_vals(i), ref_quad_vals(i)
  enddo
end subroutine dump_log_bit_grid


!> Check dump inputs on the host before formatting or reference evaluation.
subroutine check_log_dump_grid(x, x_bits, start_bits, end_bits)
  real, intent(in) :: x(:)
    !< Generated inputs [nondim]
  integer(kind=int_kind), intent(in) :: x_bits(:)
    !< Generated input bit patterns
  integer(kind=int_kind), intent(in) :: start_bits
    !< Expected start bit pattern
  integer(kind=int_kind), intent(in) :: end_bits
    !< Expected end bit pattern

  integer(kind=int_kind) :: denom, step, rem, offset, expected_bits, actual_bits
    !< Host grid spacing, coordinate, and expected and observed bit patterns.
  integer :: i

  denom = int(size(x) - 1, int_kind)
  step = (end_bits - start_bits) / denom
  rem = modulo(end_bits - start_bits, denom)
  do i = 1, size(x)
    offset = int(i - 1, int_kind)
    expected_bits = (start_bits + offset * step) + (offset * rem) / denom
    actual_bits = transfer(x(i), 0_int_kind)
    if ((x_bits(i) /= expected_bits) .or. (actual_bits /= expected_bits)) then
      print '(2x,"log dump input mismatch at index ",i0)', i
      print '(4x,"expected bits = ",Z16.16)', expected_bits
      print '(4x,"grid bits     = ",Z16.16)', x_bits(i)
      print '(4x,"x bits        = ",Z16.16)', actual_bits
      call assert(.false., "MOM_log_repro_tests: check_log_dump_grid: input grid differs from host formula")
      exit
    endif
  enddo
end subroutine check_log_dump_grid


!> Return the binary64 bit pattern of x formatted as 16 hexadecimal digits.
function real_bits_hex(x) result(hex_string)
  real, intent(in) :: x
    !< Value to format [nondim]
  character(len=16) :: hex_string
    !< Binary representation of x as hexadecimal digits.

  write(hex_string, '(Z16.16)') transfer(x, 0_int_kind)
end function real_bits_hex


!> Add log_repro tests to the intrinsic function test suite.
subroutine add_log_repro_tests(suite)
  type(TestSuite), intent(inout) :: suite

  ! Elemental test
  call suite%add(test_log_elemental, "test_log_elemental")

  ! Basic value tests
  call suite%add(test_log_one, "test_log_one")
  call suite%add(test_log_general, "test_log_general")
  call suite%add(test_log_two, "test_log_two")
  call suite%add(test_log_half, "test_log_half")
  call suite%add(test_log_sqrt2, "test_log_sqrt2")
  call suite%add(test_log_inv_sqrt2, "test_log_inv_sqrt2")
  call suite%add(test_log_ten, "test_log_ten")
  call suite%add(test_log_pos_inf, "test_log_pos_inf")
  call suite%add(test_log_zero, "test_log_zero")
  call suite%add(test_log_neg_zero, "test_log_neg_zero")
  call suite%add(test_log_neg_inf, "test_log_neg_inf")
  call suite%add(test_log_negative, "test_log_negative")
  call suite%add(test_log_nan, "test_log_nan")
  call suite%add(test_log_neg_nan, "test_log_neg_nan")
  call suite%add(test_log_signaling_nan, "test_log_signaling_nan")
  call suite%add(test_log_neg_signaling_nan, "test_log_neg_signaling_nan")
  call suite%add(test_log_largest_float, "test_log_largest_float")
  call suite%add(test_log_smallest_float, "test_log_smallest_float")
  call suite%add(test_log_smallest_subnormal, "test_log_smallest_subnormal")
  call suite%add(test_log_near_one_above, "test_log_near_one_above")
  call suite%add(test_log_near_one_below, "test_log_near_one_below")
  call suite%add(test_log_near_one_sensitive_bit, "test_log_near_one_sensitive_bit")

  ! Log property tests
  call suite%add(test_log_product_property, "test_log_product_property")
  call suite%add(test_log_reciprocal_property, "test_log_reciprocal_property")
  call suite%add(test_log_do_concurrent_bitwise, "test_log_do_concurrent_bitwise")
  call suite%add(test_log_value_dump, "test_log_value_dump")

  ! Evaluate log_repro error if quad precision is available
  if (realquad >= 0) then
    call suite%add(test_log_ulp_accuracy, "test_log_ulp_accuracy")
    call suite%add(test_log_near_one_sweep, "test_log_near_one_sweep")
    call suite%add(test_log_near_one_boundaries, "test_log_near_one_boundaries")
  else
    print '(1x,a)', 'Skipping log_repro near-1 sweep and ULP accuracy test: quad precision is unavailable.'
  endif
end subroutine add_log_repro_tests

end module MOM_log_repro_tests
