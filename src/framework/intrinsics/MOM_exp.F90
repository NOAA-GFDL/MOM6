! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of exp().

#include "MOM_exp.h"

submodule (MOM_intrinsic_functions) MOM_exp

use, intrinsic :: iso_fortran_env, only : int32, int64
use, intrinsic :: ieee_arithmetic, only : ieee_rint
use MOM_exp_data_n128, only : ndiv, idiv_scale_lookup, idiv_residual_lookup

implicit none

real, parameter :: real_mold = 0.
  !< Real mold for transfers and numerical format queries

integer, parameter :: int_kind &
    = merge(int64, int32, storage_size(real_mold) > storage_size(0_int32))
  !< Integer kind with the same storage size as default real
integer(kind=int_kind), parameter :: int_mold = 0
  !< Integer mold value

! Floating point layout
integer, parameter :: expbit = digits(real_mold) - 1
  !< Position of lowest exponent bit
  ! NOTE: digits() includes the implicit leading digit
integer, parameter :: signbit = storage_size(real_mold) - 1
  !< Position of sign bit
integer, parameter :: expwidth = signbit - expbit
  !< Number of exponent bits
integer, parameter :: expbias = maxexponent(real_mold) - 1
  !< Exponent bias

! IEEE 754 special values
integer(kind=int_kind), parameter :: pos_inf_bits &
    = ishft(2_int_kind**expwidth - 1_int_kind, expbit)
  !< IEEE +Inf bit pattern
integer(kind=int_kind), parameter :: neg_inf_bits &
    = ior(pos_inf_bits, ishft(-1_int_kind, signbit))
  !< IEEE -Inf bit pattern

! Fast integer rounding offset
real, parameter :: round_bias = 1.5 * 2_int_kind**(digits(real_mold) - 1)
  !< Binary offset used to trigger rounding of fractional values
integer(int_kind), parameter :: round_bias_bits &
    = transfer(round_bias, int_mold)
  !< Bit representation of rounding bias

contains

!> Reproducible exponential function
!!
!! Compute exp(x) with bitwise reproducibility across platforms.
module procedure exp_repro
  ! ln2 estimates
  real, parameter :: ln2 = 0.693147180559945309417232121458176568
    !< ln2: 0.693147180559945309417232... [nondim]
  real, parameter :: I_ln2 = 1.44269504088896340735992468100189214
    !< 1 / ln2: 1.4426950408889634073599... [nondim]

  ! The max and min x values between which exp(x) remains finite.
  real, parameter :: xmax &
      = real(maxexponent(real_mold) + 1) * ln2
    !< Largest x value before exp(x) overflow [nondim]
  real, parameter :: xmin &
      = real(minexponent(real_mold) - digits(real_mold) - 1) * ln2
    !< Smallest x value before exp(x) underflow [nondim]

  ! Double-real precision of ln2 used in Cody-Waite range reduction
  ! NOTE: This split assumes real64 precision.
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Upper 32 bits of ln2: 6.93147180369123816490e-01 [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Lower precision bits of ln2: 1.90821492927058770002e-10 [nondim]

  ! Subdivide [-ln2/2, ln2/2] into ndiv subintervals to reduce approximation range.
  ! This allows for a smaller polynomial at the cost of a lookup table.
  ! ndiv and the lookup tables are imported from MOM_exp_data_r64.
  real, parameter :: I_ndiv = 1. / real(ndiv)
    !< 1 / ndiv [nondim]
  real, parameter :: n_ln2 = ndiv * I_ln2
    !< ndiv / ln2 [nondim]
  real, parameter :: ln2_ndiv_hi = I_ndiv * ln2_hi
    !< Upper 32 bits of ln2 / ndiv [nondim]
  real, parameter :: ln2_ndiv_lo = I_ndiv * ln2_lo
    !< Lower precision bits of ln2 / ndiv [nondim]
  integer(int_kind), parameter :: idiv_mask = int(ndiv - 1, int_kind)
    !< Used for fast modulo of ndiv

  ! Range of K = nint(x / ln2) for which direct exponent scaling is safe.
  ! Beyond this range, a bias is applied to handle subnormals and overflow.
  ! NOTE: Fortran exponent is defined as one less than IEEE exponent.
  integer, parameter :: Kmin = minexponent(real_mold) + 1
    !< Minimum K before subnormal scaling is needed
    !! Kmin = (minexponent() - 1) + 1 (min exp(r)) (+1 safety buffer)
  integer, parameter :: Kmax = maxexponent(real_mold) - 2
    !< Maximum K before overflow scaling is needed
    !! Kmax = (maxexponent() - 1) - 0 (max exp(r)) (-1 safety buffer)
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used for overflow and subnormal scaling
    !! Any bias which rescales 2**K exp(r) to O(1) works here.

  ! Nonfinite testing
  logical :: nonfinite
    ! True if input is a nonfinite float (+/-Inf, NaN)
  integer(kind=int_kind) :: xb
    ! Bit representation of x

  ! Range-reduction variables

  real :: xc
    ! x clamped between xmin and xmax [nondim]
  integer(kind=int_kind) :: K
    ! Nearest IEEE-rounded integer to x/ln2 [nondim]
  real :: Z
    ! Nearest integer to ndiv * K + idiv for ndiv subdivisions
    ! NOTE: Z is stored as real to avoid int-real type conversions.
  integer(kind=int_kind) :: Zi
    ! Integer representation of Z
  integer :: idiv
    ! Subdivision index
  real :: r
    ! Range-reduced input, r = x - Z ln2 / ndiv [nondim]

  ! Polynomial estimation variables

  real :: e
    ! Exponent of range-reduced input, e = 2**(idiv/ndiv) exp(r) [nondim]
  real :: expm1_r
    ! Approximation to exp(r) - 1 [nondim]
  real :: idiv_scale
    ! Estimate of 2**(idiv/ndiv) [nondim]
  real :: idiv_residual
    ! Relative residual, 2**(j/N) = idiv_scale * (1 + idiv_residual) [nondim]

  ! Descaling and subnormal handling

  integer(kind=int_kind) :: eb
    ! Bit representations of e
  integer(kind=int_kind) :: j
    ! Bias added to K to compensate for exponent K beyond {-1022,..,+1023}.
  integer(kind=int_kind) :: fb
    ! Bit representation 2**j, the K exponent rescale

  ! 1. Nonfinite handling
  ! ---------------------
  ! Nonfinites must be handled first to prevent their appearance in
  ! calculations, which may raise unwanted floating point signals.

  xb = transfer(x, int_mold)
  nonfinite = iand(xb, pos_inf_bits) == pos_inf_bits

  if (nonfinite) then
    ! exp(-Inf) = 0, otherwise pass-through +Inf and +/-NaN values
    ! Compute x + x to trigger `Invalid` for signaled NaNs.
    a = merge(0., x + x, xb == neg_inf_bits)
    return
  endif

  ! 2. Range Reduction
  ! ------------------
  ! Apply a range reduction of r = x - K ln2 - (idiv / ndiv) ln2, so that
  !     exp(x) = 2**K 2**(idiv / ndiv) exp(r).
  ! If K = nint(x / ln2) then r is in [-ln2/(2*ndiv), ln2/(2*ndiv)] and exp(r)
  ! can be estimated by a sufficiently accurate polynomial.

  ! Clamp x to [xmin,xmax] to avoid extreme exponents in subnormal handler.
  xc = min(max(x, xmin), xmax)

  ! Compute Z = ndiv K + idiv, where r = x - Z (ln 2 / ndiv)
  Z = NEAREST_INT(xc * n_ln2)
  Zi = transfer(Z + round_bias, int_mold) - round_bias_bits

  ! Compute the subdivision index and integer offset K
  idiv = iand(Zi, idiv_mask)
  K = (Zi - int(idiv, int_kind)) / ndiv

  ! Since Z ~ x N / ln2, the terms in r will nearly cancel and there is some
  ! expected loss of precision.  To compensate, we use a Cody-Waite correction.
  r = (xc - Z * ln2_ndiv_hi) - Z * ln2_ndiv_lo

  ! 3. Polynomial approximation
  ! ---------------------------
  ! a = 2**K e where e = 2**(idiv/ndiv) exp(r), and exp(r) = 1 + r * P(r) where
  ! P(r) is an order-5 Remez minimax polynomial of expm1(r) / r.

  idiv_scale = idiv_scale_lookup(idiv)
  idiv_residual = idiv_residual_lookup(idiv)
  expm1_r = exp_remez_expm1_estrin_4(r)

  ! Evaluate the small correction before the final addition to idiv_scale
  e = idiv_scale + idiv_scale * (idiv_residual + r * expm1_r)

  ! 4. Unscaling
  ! ------------
  ! Compute a = 2**K e, an exact power-of-2 calculation.
  ! Adjust scaling to compensate for subnormal output.

  ! exp(r) has range [0.707/ndiv, 1.414/ndiv], so K shifts by either 0 or -1.
  ! Resolved exponent are in the range {-1022..1023}, so for a to be resolved,
  ! K must be in {-1021,1023}.  (For safety, we actually do {-1020,1022}.)

  ! Determine if K is outside the supported exponent range.
  ! If so, then apply a bias j to normalize the exponent.
  ! Kbias is chosen so that the exponent is "something near 1".
  j = merge(Kbias, 0_int_kind, K < Kmin) + merge(-Kbias, 0_int_kind, K > Kmax)

  ! Get the bit representation of e
  eb = transfer(e, int_mold)

  ! Rescale to e to exp(x), possibly including the j bias.
  eb = eb + ishft(k + j, expbit)
  a = transfer(eb, real_mold)

  ! Undo the 2**j bias as floating point multiplication.
  ! - For "normals", this has no effect.
  ! - For subnormals, this will force subnormal estimation (if enabled).
  ! - For resolvable K beyond this range, it triggers an over/underflow.
  ! - Extreme values of K have already been filtered out by the min/max step.
  fb = ishft(int(expbias, int_kind) - j, expbit)
  a = a * transfer(fb, real_mold)
end procedure exp_repro


!> Remez polynomial estimate of (exp(x) - 1) / x over [-ln2/256, ln2/256].
!! Coefficients are generated by Sollya 8, and evaluation is in Estrin form.
pure function exp_remez_expm1_estrin_4(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-ln2/256, ln2/256] [nondim]
  real :: e
    !< Approximation of (exp(x) - 1) / x [nondim]

  real :: x2, x4
    !< Powers of x [nondim]
  real :: p01, p23, p
    !< Polynomial partial sums [nondim]

  !> Remez coefficients for (exp(x) - 1) / x on [-ln2/256, ln2/256] [nondim]
  real, parameter :: c(0:4) = [ &
      1.0, &
      0.4999999999999766853164828717126511037349700927734375, &
      0.166666666666670015839457619222230277955532073974609375, &
      4.1666679392304360740606483659576042555272579193115234375e-2, &
      8.3333340579158539374038383584775147028267383575439453125e-3 &
  ]

  x2 = x * x
  x4 = x2 * x2

  p01 = c(0) + c(1) * x
  p23 = c(2) + c(3) * x

  ! Final assembly: (exp(x) - 1) / x
  e = (p01 + x2 * p23) + x4 * c(4)
end function exp_remez_expm1_estrin_4


!> An optimized nearest-integer function for floating point reals.
!!
!! The value x is shifted from 2**K (1 + a) to 2**(digits-1+K) (1.5 + [a]),
!! causing the fractional part to be rounded according to the current IEEE
!! settings.  In almost all cases, this is nearest ties-to-even.
!!
!! The +0.5 ensures that the biased exponent of negative numbers does not drop
!! by one, which can cause half-value rounding.
!!
!! The behavior of this function does not match nint() or anint().  The nint()
!! function always ties away from zero, e.g. nint(2.5) = 3.
!!
!! It is essential that compilers not reduce (x+b)-b to x.  This can typically
!! be ensured as long as parentheses are respected.  This is managed by the
!! ENABLE_FAST_RINT macro in MOM_exp.h and assigned to NEAREST_INT().  If
!! unset, then ieee_rint() is used.
pure function fast_rint(x) result(n)
  real, intent(in) :: x
    !< Real value to be rounded to the nearest integer
  real :: n
    !< Nearest integer to x, stored as a real

  n = (x + round_bias) - round_bias
end function fast_rint

end submodule MOM_exp
