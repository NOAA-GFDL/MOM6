! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Unit tests for MOM_intrinsic_functions
module MOM_intrinsic_functions_tests

use, intrinsic :: ieee_arithmetic, only : ieee_value
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_signaling_nan
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use, intrinsic :: ieee_exceptions, only : ieee_set_flag, ieee_get_flag
use, intrinsic :: ieee_exceptions, only : ieee_support_flag
use, intrinsic :: ieee_exceptions, only : ieee_all
use, intrinsic :: ieee_exceptions, only : ieee_invalid, ieee_overflow
use, intrinsic :: ieee_exceptions, only : ieee_underflow, ieee_inexact
use, intrinsic :: ieee_exceptions, only : ieee_divide_by_zero

use MOM_error_handler, only : assert
use MOM_unit_testing, only : TestSuite
use MOM_intrinsic_functions, only : exp_repro

implicit none ; private

public :: run_intrinsic_functions_tests

integer, parameter :: realquad = selected_real_kind(p=30, r=300)
  !< Potential real128 precision kind.  If unavailable, this will be negative.
integer, parameter :: realq = merge(realquad, kind(1.), realquad >= 0.)
  !< Placeholder real128 for declarations.  Unused if real128 is unavailable.

real, parameter :: rmold = 0.
  !< Mold value for default real

! Mathematical constants
real, parameter :: ln2 = 0.69314718055994531
  !< ln(2)
real, parameter :: half_ln2 = 0.34657359027997265
  !< ln(2)/2
real, parameter :: sqrt2 = 1.4142135623730950
  !< sqrt(2)
real, parameter :: exp1 = 2.7182818284590452
  !< e
real, parameter :: inv_e = 0.36787944117144232
  !< 1/e
real, parameter :: exp_033 = 1.3909681284637803
  !< exp(0.33)
real, parameter :: exp_tiny_pos = 1.0000000000000010
  !< exp(1e-15)
real, parameter :: exp_tiny_neg = 0.99999999999999900
  !< exp(-1e-15)
real, parameter :: exp_709_78 = 1.7928227943945156e+308
  !< exp(709.78)
real, parameter :: exp_neg_708_39 = 2.2394014988804678e-308
  !< exp(-708.39)
real, parameter :: exp_neg2 = 0.13533528323661269
  !< exp(-2)
real, parameter :: exp2 = 7.3890560989306502
  !< exp(2)
real, parameter :: log_huge = 709.78271289338397
  !< log(huge())
real, parameter :: exp_log_huge = 1.7976931348622732e+308
  !< exp(log_huge)
real, parameter :: log_tiny = -708.39641853226408
  !< log(tiny())
real, parameter :: exp_log_tiny = 2.2250738585072625e-308
  !< exp(log_tiny)

! Module-level flag to enable/disable IEEE exception tests
logical :: ieee_flags_supported = .false.
  !< True if the platform supports IEEE exception flag queries

contains

!> Check if IEEE exception flags are supported on this platform
subroutine check_ieee_support
  ieee_flags_supported = ieee_support_flag(ieee_overflow) &
      .and. ieee_support_flag(ieee_underflow) &
      .and. ieee_support_flag(ieee_inexact) &
      .and. ieee_support_flag(ieee_invalid) &
      .and. ieee_support_flag(ieee_divide_by_zero)
end subroutine check_ieee_support


!> Test that exp_repro(0) = 1 exactly
subroutine test_exp_zero
  real :: val

  val = exp_repro(0.)

  print '(2x, "exp_repro(0) = ", ES22.15)', val

  call assert(val == 1., "exp_repro(0) should equal 1 exactly")
end subroutine test_exp_zero


!> Test that exp_repro(-0) = 1 exactly
subroutine test_exp_neg_zero
  real :: val

  val = exp_repro(-0.)

  print '(2x, "exp_repro(-0) = ", ES22.15)', val

  call assert(val == 1., "exp_repro(-0) should equal 1 exactly")
end subroutine test_exp_neg_zero


!> Test exp_repro(1) against reference value
subroutine test_exp_one
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(1.)
  err = abs(val - exp1) / exp1

  print '(2x, "exp_repro(1) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(1) relative error exceeds tolerance")
end subroutine test_exp_one


!> Test exp_repro(-1) against reference value
subroutine test_exp_neg_one
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(-1.)
  err = abs(val - inv_e) / inv_e

  print '(2x, "exp_repro(-1) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(-1) relative error exceeds tolerance")
end subroutine test_exp_neg_one


!> Test exp_repro at a general value (0.33)
subroutine test_exp_general
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-15

  x = 0.33
  ref = exp_033
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(0.33) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(0.33) relative error exceeds tolerance")
end subroutine test_exp_general


!> Test exp_repro at range reduction boundary: ln(2)
subroutine test_exp_ln2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(ln2)
  err = abs(val - 2.0) / 2.0

  print '(2x, "exp_repro(ln2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(ln(2)) should be close to 2")
end subroutine test_exp_ln2


!> Test exp_repro at range reduction boundary: -ln(2)
subroutine test_exp_neg_ln2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(-ln2)
  err = abs(val - 0.5) / 0.5

  print '(2x, "exp_repro(-ln2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(-ln(2)) should be close to 0.5")
end subroutine test_exp_neg_ln2


!> Test exp_repro at half the range reduction boundary: 0.5*ln(2)
subroutine test_exp_half_ln2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(half_ln2)
  err = abs(val - sqrt2) / sqrt2

  print '(2x, "exp_repro(ln2/2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(0.5*ln(2)) should be close to sqrt(2)")
end subroutine test_exp_half_ln2


!> Test exp_repro with tiny argument (tests 1 + x accuracy)
subroutine test_exp_tiny_pos
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-15

  x = 1.0e-15
  ref = exp_tiny_pos
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(1e-15) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(1e-15) relative error exceeds tolerance")
end subroutine test_exp_tiny_pos


!> Test exp_repro with tiny negative argument
subroutine test_exp_tiny_neg
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-15

  x = -1.0e-15
  ref = exp_tiny_neg
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(-1e-15) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(-1e-15) relative error exceeds tolerance")
end subroutine test_exp_tiny_neg


!> Test exp_repro overflow behavior: exp(1000)
subroutine test_exp_overflow
  real :: x, val, ref

  x = 1000.
  ref = ieee_value(0., ieee_positive_inf)
  val = exp_repro(x)

  print '(2x, "exp_repro(1000) = ", ES22.15)', val

  call assert(val == ref, "exp_repro(1000) should overflow to +Inf")
end subroutine test_exp_overflow


!> Test exp_repro extreme overflow behavior: exp(10000)
subroutine test_exp_extreme_overflow
  real :: x, val, ref

  x = 10000.
  ref = ieee_value(0., ieee_positive_inf)
  val = exp_repro(x)

  print '(2x, "exp_repro(10000) = ", ES22.15)', val

  call assert(val == ref, "exp_repro(10000) should overflow to +Inf")
end subroutine test_exp_extreme_overflow


!> Test exp_repro underflow behavior: exp(-1000)
subroutine test_exp_underflow
  real :: x, val, ref

  x = -1000.
  ref = 0.
  val = exp_repro(x)

  print '(2x, "exp_repro(-1000) = ", ES22.15)', val

  call assert(val == ref, "exp_repro(-1000) should underflow to 0")
end subroutine test_exp_underflow


!> Test exp_repro extreme underflow behavior: exp(-10000)
subroutine test_exp_extreme_underflow
  real :: x, val, ref

  x = -10000.
  ref = 0.
  val = exp_repro(x)

  print '(2x, "exp_repro(-10000) = ", ES22.15)', val

  call assert(val == ref, "exp_repro(-10000) should underflow to 0")
end subroutine test_exp_extreme_underflow


!> Test exp_repro near overflow boundary (709.78)
subroutine test_exp_near_overflow
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-14

  x = 709.78
  ref = exp_709_78
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(709.78) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(709.78) relative error exceeds tolerance")
end subroutine test_exp_near_overflow


!> Test exp_repro near underflow boundary (-708.39)
subroutine test_exp_near_underflow
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-14

  x = -708.39
  ref = exp_neg_708_39
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(-708.39) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(-708.39) relative error exceeds tolerance")
end subroutine test_exp_near_underflow


!> Test exp_repro in subnormal region (-745.13)
subroutine test_exp_subnormal
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-10  ! Relaxed for subnormals

  x = -745.13
  ref = nearest(0., 1.)
  val = exp_repro(x)

  print '(2x, "exp_repro(-745.13) = ", ES22.15)', val

  ! For subnormals, check absolute error or that both are very small
  if (ref > 0.) then
    err = abs(val - ref) / ref
    call assert(err < tol, "exp_repro(-745.13) relative error exceeds tolerance")
  else
    call assert(val == 0., "exp_repro(-745.13) should be zero or subnormal")
  endif
end subroutine test_exp_subnormal


!> Test exp_repro at largest representable float boundary
subroutine test_exp_largest_float
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-13  ! Relaxed due to log/exp round-trip error

  x = log_huge
  ref = exp_log_huge
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(log(huge)) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(log(huge)) relative error exceeds tolerance")
end subroutine test_exp_largest_float


!> Test exp_repro at smallest normal float boundary
subroutine test_exp_smallest_float
  real :: val, ref, err
  real, parameter :: tol = 5.e-14

  ref = exp_log_tiny
  val = exp_repro(log_tiny)
  err = abs(val - ref) / ref

  print '(2x, "exp_repro(log(tiny)) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "exp_repro(log(tiny)) relative error exceeds tolerance")
end subroutine test_exp_smallest_float


!> Test exp_repro(+Inf) = +Inf
subroutine test_exp_pos_inf
  real :: x, val, ref

  x = ieee_value(0., ieee_positive_inf)
  ref = ieee_value(0., ieee_positive_inf)
  val = exp_repro(x)

  print '(2x, "exp_repro(+Inf) = ", ES22.15)', val

  call assert(val == ref, "exp_repro(+Inf) should equal +Inf")
end subroutine test_exp_pos_inf


!> Test exp_repro(-Inf) = 0
subroutine test_exp_neg_inf
  real :: x, val

  x = ieee_value(0., ieee_negative_inf)
  val = exp_repro(x)

  print '(2x, "exp_repro(-Inf) = ", ES22.15)', val

  call assert(val == 0., "exp_repro(-Inf) should equal 0")
end subroutine test_exp_neg_inf


!> Test exp_repro(NaN) = NaN
subroutine test_exp_nan
  real :: x, val

  x = ieee_value(0., ieee_quiet_nan)
  val = exp_repro(x)

  print '(2x, "exp_repro(NaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "exp_repro(NaN) should be NaN")
end subroutine test_exp_nan


!> Test exp_repro(-NaN) = NaN
subroutine test_exp_neg_nan
  real :: x, val

  x = -ieee_value(0., ieee_quiet_nan)
  val = exp_repro(x)

  print '(2x, "exp_repro(-NaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "exp_repro(-NaN) should be NaN")
end subroutine test_exp_neg_nan


!> Test ULP accuracy over a range of values [-10, 10]
!!
!! exp_repro should be within 1.5 ULP of the true value (1 ULP with FMA).
!! We test against 2 ULP to allow for combined error from exp_repro and the
!! quad-precision reference.
subroutine test_exp_ulp_accuracy
  integer, parameter :: npts = 100000
  real, parameter :: xmin = -10.
  real, parameter :: xmax = 10.

  ! Input axis
  real :: x(npts)
  real :: I_npts

  real :: val(npts), val_vec(npts)
  real :: val_exp(npts), val_exp_vec(npts)
  real(kind=realq) :: val_quad(npts), val_quad_vec(npts)

  integer :: i

  ! Generate test points
  I_npts = 1. / (npts - 1)
  do i = 1, npts
    x(i) = xmin + (i - 1) * ((xmax - xmin) * I_npts)
  enddo

  ! Several libraries have scalar and vector implementations, chosen at the
  ! discretion of the compiler.  The following attempts to test each case.

  ! Scalar evaluations
  do i = 1, npts
    val_exp(i) = exp(x(i))
    val(i) = exp_repro(x(i))
    val_quad(i) = exp(real(x(i), realq))

    ! Impossible branch to prevent vectorization
    if (val(i) < 0.) exit
  enddo

  ! Vector-favorable evaluation
  val_exp_vec = exp(x)
  val_vec = exp_repro(x)
  val_quad_vec = exp(real(x, realq))

  ! Assert that exp_repro() is within 2 ULP.
  print '(1x,a)', '=== scalar exp_repro() accuracy'
  call check_ulp_accuracy(x, val, val_quad, max_ulp_tol=2.)

  ! We expect scalar and vector implementations to agree.
  call assert(all(val == val_vec), 'Scalar and vector exp_repro() do not agree')
  print '(1x,a)', '=== vector exp_repro() matches scalar'

  ! exp() accuracy is provided for comparison.
  print '(1x,a)', '=== scalar exp() accuracy'
  call check_ulp_accuracy(x, val_exp, val_quad)

  if (all(val_exp == val_exp_vec)) then
    print '(1x,a)', '=== vector exp() matches scalar'
  else
    print '(1x,a)', '=== vector exp() accuracy'
    call check_ulp_accuracy(x, val_exp_vec, val_quad_vec)
  endif
end subroutine test_exp_ulp_accuracy


!> Compute the function accuracy relative to a real128-precision reference.
!! Absolute, relative, and ULP error is computed, as well as the number of
!! points above 0.5 and 1 ULP.  An error is raised if any point exceeds the
!! optionally prescribed maximum ULP tolerance.
subroutine check_ulp_accuracy(x, val, ref, max_ulp_tol)
  real, intent(in) :: x(:)
    !< Input grid
  real, intent(in) :: val(:)
    !< Output estimates
  real(kind=realq), intent(in) :: ref(:)
    !< Reference estimates in real128 precision
  real, optional, intent(in) :: max_ulp_tol
    !< Maximum ULP tolerance

  real(kind=realq) :: err, rel_err, ulp_err
    ! Absolute, relative, and ULP error of val relative to ref.

  real :: max_abs_err, max_rel_err, max_ulp
  real :: sum_abs_err, sum_rel_err, sum_sq_err
  real :: x_max_abs, x_max_rel, x_max_ulp
  real :: ulp_val

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

  print '(2x,"Tested ", i0, " points in [-10, 10]")', npts

  ! NOTE: Floats use t25 assuming a positive sign as blank
  print '(2x,"max abs err:", t25, ES12.5, " at x = ", f10.4)', &
      max_abs_err, x_max_abs
  print '(2x,"max rel err:", t25, ES12.5, " at x = ", f10.4)', &
      max_rel_err, x_max_rel
  print '(2x,"max ULP err (vs quad):", t26, f12.10, " at x = ", f10.4)', &
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
    call assert(max_ulp < max_ulp_tol, &
        "exp_repro max ULP error exceeds 2 over [-10, 10]")
  endif
end subroutine check_ulp_accuracy


!> Test the exponential property: exp(r1) * exp(r2) = exp(r1 + r2)
!!
!! This property test verifies that exp_repro produces consistent results
!! by checking that the product of two exponentials equals the exponential
!! of the sum, within floating-point tolerance.
subroutine test_exp_product_property
  integer, parameter :: npts = 1000
  real, parameter :: tol = 1.e-14

  real :: r1, r2
  real :: exp_r1, exp_r2, exp_sum
  real :: product_result, sum_result
  real :: err, max_err
  real :: r1_max, r2_max
  integer :: i
  real :: seed

  max_err = 0.

  ! Use a simple deterministic sequence for reproducibility
  seed = 0.123456789

  do i = 1, npts
    ! Generate pseudo-random values in a range where sum won't overflow
    ! Keep r1, r2 in [-5, 5] so r1+r2 in [-10, 10]
    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r1 = (seed / 2.**31) * 10. - 5.

    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r2 = (seed / 2.**31) * 10. - 5.

    exp_r1 = exp_repro(r1)
    exp_r2 = exp_repro(r2)
    exp_sum = exp_repro(r1 + r2)

    product_result = exp_r1 * exp_r2
    sum_result = exp_sum

    if (sum_result /= 0.) then
      err = abs(product_result - sum_result) / abs(sum_result)
      if (err > max_err) then
        max_err = err
        r1_max = r1
        r2_max = r2
      endif
    endif
  enddo

  print '("Tested ", i0, " random (r1, r2) pairs")', npts
  print '("max rel err in exp(r1)*exp(r2) vs exp(r1+r2):", t48, ES12.5)', max_err
  print '("  at r1 = ", f10.6, ", r2 = ", f10.6)', r1_max, r2_max

  call assert(max_err < tol, "exp_repro product property test failed")
end subroutine test_exp_product_property


!> Test the exponential property: exp(-x) = 1/exp(x)
!!
!! This property test verifies that exp_repro(-x) equals the reciprocal
!! of exp_repro(x), within floating-point tolerance.
subroutine test_exp_negation_property
  integer, parameter :: npts = 1000
  real, parameter :: tol = 1.e-14

  real :: x
  real :: exp_neg_x, inv_exp_x
  real :: err, max_err
  real :: x_max
  integer :: i
  real :: seed

  max_err = 0.

  ! Use a simple deterministic sequence for reproducibility
  seed = 0.987654321

  do i = 1, npts
    ! Generate pseudo-random values in a range that avoids overflow/underflow
    ! Keep x in [-300, 300] so both exp(x) and exp(-x) are representable
    seed = mod(seed * 1103515245. + 12345., 2.**31)
    x = (seed / 2.**31) * 600. - 300.

    exp_neg_x = exp_repro(-x)
    inv_exp_x = 1. / exp_repro(x)

    if (inv_exp_x /= 0. .and. exp_neg_x /= 0.) then
      err = abs(exp_neg_x - inv_exp_x) / abs(inv_exp_x)
      if (err > max_err) then
        max_err = err
        x_max = x
      endif
    endif
  enddo

  print '("Tested ", i0, " random x values")', npts
  print '("max rel err in exp(-x) vs 1/exp(x):", t48, ES12.5)', max_err
  print '("  at x = ", f10.6)', x_max

  call assert(max_err < tol, "exp_repro negation property test failed")
end subroutine test_exp_negation_property


!> Test that exp_repro is elemental (works on arrays)
subroutine test_exp_elemental
  real :: x(5), val(5), ref(5), err
  real, parameter :: tol = 1.e-14
  integer :: i

  x(:) = [-2., -1., 0., 1., 2.]
  ref(:) = [exp_neg2, inv_e, 1., exp1, exp2]

  val= exp_repro(x)

  do i=1,5
    if (ref(i) /= 0.) then
      err = abs(val(i) - ref(i)) / abs(ref(i))
      call assert(err < tol, "exp_repro elemental test failed")
    endif
  enddo
end subroutine test_exp_elemental


!> Test IEEE exception flags for normal input (should raise inexact only)
subroutine test_exp_flags_normal
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = 0.33

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(0.33) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(0.33) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(0.33) should not raise underflow")
  call assert(flag_inexact, "exp_repro(0.33) should raise inexact")
end subroutine test_exp_flags_normal


!> Test IEEE exception flags for overflow (should raise overflow + inexact)
subroutine test_exp_flags_overflow
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = 1000.

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(1000) should not raise invalid")
  call assert(flag_overflow, "exp_repro(1000) should raise overflow")
  call assert(.not. flag_underflow, "exp_repro(1000) should not raise underflow")
  call assert(flag_inexact, "exp_repro(1000) should raise inexact")
end subroutine test_exp_flags_overflow


!> Test IEEE exception flags for underflow (should raise underflow + inexact)
subroutine test_exp_flags_underflow
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = -1000.

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(-1000) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(-1000) should not raise overflow")
  call assert(flag_underflow, "exp_repro(-1000) should raise underflow")
  call assert(flag_inexact, "exp_repro(-1000) should raise inexact")
end subroutine test_exp_flags_underflow


!> Test IEEE exception flags for +Inf input (should raise no flags)
subroutine test_exp_flags_pos_inf
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = ieee_value(0., ieee_positive_inf)

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(+Inf) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(+Inf) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(+Inf) should not raise underflow")
  call assert(.not. flag_inexact, "exp_repro(+Inf) should not raise inexact")
end subroutine test_exp_flags_pos_inf


!> Test IEEE exception flags for -Inf input (should raise no flags)
subroutine test_exp_flags_neg_inf
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = ieee_value(0., ieee_negative_inf)

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(-Inf) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(-Inf) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(-Inf) should not raise underflow")
  call assert(.not. flag_inexact, "exp_repro(-Inf) should not raise inexact")
end subroutine test_exp_flags_neg_inf


!> Test IEEE exception flags for NaN input (should raise no flags for quiet NaN)
subroutine test_exp_flags_nan
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = ieee_value(0., ieee_quiet_nan)

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(NaN) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(NaN) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(NaN) should not raise underflow")
  call assert(.not. flag_inexact, "exp_repro(NaN) should not raise inexact")
end subroutine test_exp_flags_nan


!> Print a summary table of IEEE exception flags for various inputs
!!
!! This is informational only - it shows which flags are raised by exp() and
!! exp_repro() for different input categories. Not enforced by assertions.
subroutine print_ieee_flags_summary
  real, volatile :: x, val_intrinsic, val_repro
  logical :: flags_intrinsic(5), flags_repro(5)
  character(len=5) :: flag_str_intrinsic, flag_str_repro
  integer :: i

  ! Test cases: name, input value
  character(len=20) :: test_names(12)
  real :: test_values(12)

  if (.not. ieee_flags_supported) then
    print '(2x, "IEEE flags not supported on this platform")'
    return
  endif

  ! Define test cases
  test_names(1) = "exact (0)"
  test_values(1) = 0.

  test_names(2) = "normal"
  test_values(2) = 0.33

  test_names(3) = "overflow"
  test_values(3) = 1000.

  test_names(4) = "underflow"
  test_values(4) = -1000.

  test_names(5) = "near overflow"
  test_values(5) = 709.78

  test_names(6) = "near underflow"
  test_values(6) = -708.39

  test_names(7) = "largest float"
  test_values(7) = log_huge

  test_names(8) = "smallest normal"
  test_values(8) = log_tiny

  test_names(9) = "+Inf"
  test_values(9) = ieee_value(0., ieee_positive_inf)

  test_names(10) = "-Inf"
  test_values(10) = ieee_value(0., ieee_negative_inf)

  test_names(11) = "NaN"
  test_values(11) = ieee_value(0., ieee_quiet_nan)

  test_names(12) = "sNaN"
  test_values(12) = ieee_value(0., ieee_signaling_nan)

  print '(/, "IEEE flags summary:", 1x, "exp()", 3x, "exp_repro")'
  print '(a18, 2x, "IOUXZ", 3x, "IOUXZ")', ""

  do i = 1, 12
    x = test_values(i)

    ! Test intrinsic exp()
    call ieee_set_flag(ieee_all, .false.)
    val_intrinsic = exp(x)
    call get_flag_string(flags_intrinsic, flag_str_intrinsic)

    ! Test exp_repro()
    call ieee_set_flag(ieee_all, .false.)
    val_repro = exp_repro(x)
    call get_flag_string(flags_repro, flag_str_repro)

    print '(a18, ": ", a5, 3x, a5)', &
        trim(test_names(i)), flag_str_intrinsic, flag_str_repro
  enddo

contains

  subroutine get_flag_string(flags, str)
    logical, intent(out) :: flags(5)
    character(len=5), intent(out) :: str

    call ieee_get_flag(ieee_invalid, flags(1))
    call ieee_get_flag(ieee_overflow, flags(2))
    call ieee_get_flag(ieee_underflow, flags(3))
    call ieee_get_flag(ieee_inexact, flags(4))
    call ieee_get_flag(ieee_divide_by_zero, flags(5))

    str(1:1) = merge('I', '.', flags(1))
    str(2:2) = merge('O', '.', flags(2))
    str(3:3) = merge('U', '.', flags(3))
    str(4:4) = merge('X', '.', flags(4))
    str(5:5) = merge('Z', '.', flags(5))
  end subroutine get_flag_string

end subroutine print_ieee_flags_summary


!> Run all intrinsic function tests
subroutine run_intrinsic_functions_tests
  type(TestSuite) :: suite

  ! Check IEEE exception support before running tests
  call check_ieee_support

  suite = TestSuite()

  ! Basic value tests
  call suite%add(test_exp_zero, "test_exp_zero")
  call suite%add(test_exp_neg_zero, "test_exp_neg_zero")
  call suite%add(test_exp_one, "test_exp_one")
  call suite%add(test_exp_neg_one, "test_exp_neg_one")
  call suite%add(test_exp_general, "test_exp_general")

  ! Range reduction boundary tests
  call suite%add(test_exp_ln2, "test_exp_ln2")
  call suite%add(test_exp_neg_ln2, "test_exp_neg_ln2")
  call suite%add(test_exp_half_ln2, "test_exp_half_ln2")

  ! Tiny argument tests
  call suite%add(test_exp_tiny_pos, "test_exp_tiny_pos")
  call suite%add(test_exp_tiny_neg, "test_exp_tiny_neg")

  ! Overflow/underflow tests
  call suite%add(test_exp_overflow, "test_exp_overflow")
  call suite%add(test_exp_extreme_overflow, "test_exp_extreme_overflow")
  call suite%add(test_exp_underflow, "test_exp_underflow")
  call suite%add(test_exp_extreme_underflow, "test_exp_extreme_underflow")
  call suite%add(test_exp_near_overflow, "test_exp_near_overflow")
  call suite%add(test_exp_near_underflow, "test_exp_near_underflow")
  call suite%add(test_exp_subnormal, "test_exp_subnormal")

  ! Boundary tests
  call suite%add(test_exp_largest_float, "test_exp_largest_float")
  call suite%add(test_exp_smallest_float, "test_exp_smallest_float")

  ! Special value tests
  call suite%add(test_exp_pos_inf, "test_exp_pos_inf")
  call suite%add(test_exp_neg_inf, "test_exp_neg_inf")
  call suite%add(test_exp_nan, "test_exp_nan")
  call suite%add(test_exp_neg_nan, "test_exp_neg_nan")

  ! Evaluate error if quad precision is available
  if (realq >= 0) then
    call suite%add(test_exp_ulp_accuracy, "test_exp_ulp_accuracy")
  endif

  ! Property tests
  call suite%add(test_exp_product_property, "test_exp_product_property")
  call suite%add(test_exp_negation_property, "test_exp_negation_property")

  ! Elemental test
  call suite%add(test_exp_elemental, "test_exp_elemental")

  ! IEEE exception flag tests (skipped if not supported)
  call suite%add(test_exp_flags_normal, "test_exp_flags_normal")
  call suite%add(test_exp_flags_overflow, "test_exp_flags_overflow")
  call suite%add(test_exp_flags_underflow, "test_exp_flags_underflow")
  call suite%add(test_exp_flags_pos_inf, "test_exp_flags_pos_inf")
  call suite%add(test_exp_flags_neg_inf, "test_exp_flags_neg_inf")
  call suite%add(test_exp_flags_nan, "test_exp_flags_nan")

  call suite%run()

  ! Print IEEE flags summary table (informational, not enforced)
  call print_ieee_flags_summary()
end subroutine run_intrinsic_functions_tests

end module MOM_intrinsic_functions_tests
