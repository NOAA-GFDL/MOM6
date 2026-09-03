! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Timing tests for MOM_intrinsic_functions (exp_repro)
program time_MOM_intrinsic_functions

use, intrinsic :: iso_fortran_env, only : int64, real64
use MOM_intrinsic_functions, only : exp_repro

implicit none

integer, parameter :: npts = 100000
  !< Number of test points
integer, parameter :: niter = 200
  !< Number of timing iterations
real, parameter :: xmin = -10.
  !< Minimum x value [nondim]
real, parameter :: xmax = 10.
  !< Maximum x value [nondim]

real, allocatable :: x(:), val(:)
integer(kind=int64) :: count_rate, c1, c2
real(kind=real64) :: clock_rate, time_intrinsic, time_repro
real(kind=real64) :: time_scalar_baseline, time_scalar_intrinsic, time_scalar_repro
real :: I_npts
real, volatile :: scalar_x, scalar_val
integer :: i, j

call system_clock(count_rate=count_rate)
clock_rate = real(count_rate, real64)

allocate(x(npts), val(npts))
I_npts = 1. / (npts - 1)
do i = 1, npts
  x(i) = xmin + (i - 1) * ((xmax - xmin) * I_npts)
enddo

print '("=== MOM_intrinsic_functions timing ===")'
print '("npts = ", i0, ", niter = ", i0)', npts, niter
print '("x range: [", f6.1, ", ", f6.1, "]")', xmin, xmax
print *

! Warm-up and time intrinsic exp()
do j = 1, 3
  val = exp(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c1)
do j = 1, niter
  val = exp(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c2)
time_intrinsic = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("exp() time/elem:", t30, f8.2, " ns")', time_intrinsic
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

! Warm-up and time exp_repro()
do j = 1, 3
  val = exp_repro(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c1)
do j = 1, niter
  val = exp_repro(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c2)
time_repro = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("exp_repro() time/elem:", t30, f8.2, " ns")', time_repro
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

print *
print '("slowdown factor:", t30, f8.2, "x")', time_repro / time_intrinsic

print *
print '("=== scalar loop-carried timing ===")'

! Warm-up and time a scalar dependency chain without exp().
! This gives an approximate loop overhead for interpreting the scalar timings.
do j = 1, 3
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = 1. + scalar_x
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = 1. + scalar_x
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo
call system_clock(count=c2)
time_scalar_baseline = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("baseline scalar time/call:", t30, f8.2, " ns")', time_scalar_baseline
print '("  final scalar value:", t30, ES12.5)', scalar_val

! Warm-up and time intrinsic exp() in a scalar loop with a carried dependency.
! This prevents the loop from being converted into independent vector lanes.
do j = 1, 3
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo
call system_clock(count=c2)
time_scalar_intrinsic = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("exp() scalar time/call:", t30, f8.2, " ns")', time_scalar_intrinsic
print '("  final scalar value:", t30, ES12.5)', scalar_val

! Warm-up and time exp_repro() in the same scalar loop pattern.
do j = 1, 3
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp_repro(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp_repro(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo
call system_clock(count=c2)
time_scalar_repro = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("exp_repro() scalar time/call:", t30, f8.2, " ns")', time_scalar_repro
print '("  final scalar value:", t30, ES12.5)', scalar_val

print *
print '("scalar slowdown factor:", t30, f8.2, "x")', &
    time_scalar_repro / time_scalar_intrinsic
print '("exp() minus baseline:", t30, f8.2, " ns")', &
    time_scalar_intrinsic - time_scalar_baseline
print '("exp_repro() minus baseline:", t30, f8.2, " ns")', &
    time_scalar_repro - time_scalar_baseline
print '("adjusted slowdown factor:", t30, f8.2, "x")', &
    (time_scalar_repro - time_scalar_baseline) &
      / (time_scalar_intrinsic - time_scalar_baseline)

deallocate(x, val)

end program time_MOM_intrinsic_functions
