! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Unit test driver for MOM_intrinsic_functions
module MOM_intrinsic_functions_tests

use, intrinsic :: ieee_arithmetic, only : ieee_value
use, intrinsic :: ieee_arithmetic, only : ieee_signaling_nan
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan
use, intrinsic :: ieee_exceptions, only : ieee_set_flag, ieee_get_flag
use, intrinsic :: ieee_exceptions, only : ieee_support_flag
use, intrinsic :: ieee_exceptions, only : ieee_all
use, intrinsic :: ieee_exceptions, only : ieee_invalid, ieee_overflow
use, intrinsic :: ieee_exceptions, only : ieee_underflow, ieee_inexact
use, intrinsic :: ieee_exceptions, only : ieee_divide_by_zero

use MOM_unit_testing, only : TestSuite
use MOM_intrinsic_functions, only : exp_repro
use MOM_exp_repro_tests, only : add_exp_repro_tests

implicit none ; private

public :: run_intrinsic_functions_tests

real, parameter :: log_huge = 709.78271289338397
  !< log(huge()) [nondim]
real, parameter :: log_tiny = -708.39641853226408
  !< log(tiny()) [nondim]

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

!> Print a summary table of IEEE exception flags for various inputs
!!
!! This is informational only - it shows which flags are raised by exp() and
!! exp_repro() for different input categories. Not enforced by assertions.
subroutine print_ieee_flags_summary
  real, volatile :: x
    ! Input [nondim]
  real, volatile :: val_intrinsic, val_repro
    ! Intrinsic and reproducible results [nondim]
  logical :: flags_intrinsic(5), flags_repro(5)
  character(len=5) :: flag_str_intrinsic, flag_str_repro
  integer :: i

  ! Test cases: name, input value
  character(len=20) :: test_names(12)
  real :: test_values(12)
    ! Test inputs [nondim]

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

  ! Check IEEE exception support before printing the informational summary.
  call check_ieee_support

  suite = TestSuite()

  call add_exp_repro_tests(suite)

  call suite%run()

  ! Print IEEE flags summary table (informational, not enforced)
  call print_ieee_flags_summary()
end subroutine run_intrinsic_functions_tests

end module MOM_intrinsic_functions_tests
