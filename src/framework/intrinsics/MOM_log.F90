! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of log().

submodule (MOM_intrinsic_functions) MOM_log

use MOM_log_data_n128, only : binbits
use MOM_log_data_n128, only : inv_c, log_c
use MOM_log_data_n128, only : c_hi, c_lo

implicit none

! The near-one bounds follow glibc.
! Both are exact binary fractions enclosing the inputs with |log(x)| <= 1/16.

real, parameter :: near_one_lo = 1. - 1. / 16.
  !< Lower endpoint of the near-one input interval [nondim]
real, parameter :: near_one_hi = 1. + 265. / 4096.
  !< Upper endpoint of the near-one input interval [nondim]

contains

!> Reproducible natural logarithm function
!!
!! Compute log(x) with bitwise reproducibility across platforms.
module procedure log_repro
  integer(kind=int_rsize), parameter :: near_lo_bits = transfer(near_one_lo, int_mold)
    !< Lower endpoint bit pattern
  integer(kind=int_rsize), parameter :: near_hi_bits = transfer(near_one_hi, int_mold)
    !< Upper endpoint bit pattern
  integer(kind=int_rsize) :: xb
    ! Bit representation of x

  xb = transfer(x, int_mold)

  ! Use log_table() for positive normal inputs away from 1.
  ! Otherwise use log_fallback():
  !   xb < expstep:
  !     zeros, negatives (including -Inf), positive subnormals, and -NaN
  !   xb >= expmask:
  !     +Inf and +NaN
  !   near_lo_bits <= xb <= near_hi_bits:
  !     near-one inputs needing the compensated polynomial

  if (xb < expstep .or. xb >= expmask &
      .or. (xb >= near_lo_bits .and. xb <= near_hi_bits)) then
    a = log_fallback(x)
  else
    a = log_table(xb, 0_int_rsize)
  endif
end procedure log_repro


!> Scalar fallback for exceptional, near-one, and subnormal inputs.
pure function log_fallback(x) result(a)
  !DIR$ ATTRIBUTES NOINLINE :: log_fallback
  real, intent(in) :: x
    !< Input [nondim]
  real :: a
    !< Natural logarithm [nondim]
  integer(kind=int_rsize), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used to normalize subnormal inputs
  integer(kind=int_rsize), parameter :: near_lo_bits = transfer(near_one_lo, int_mold)
    !< Lower endpoint bit pattern
  integer(kind=int_rsize), parameter :: near_hi_bits = transfer(near_one_hi, int_mold)
    !< Upper endpoint bit pattern
  integer(kind=int_rsize) :: xb
    ! Bit representation of x
  integer(kind=int_rsize) :: exponent_adjustment
    ! Correction to the normalized exponent for subnormal inputs
  integer(kind=int_rsize) :: raw_exp
    ! Biased IEEE exponent field
  integer(kind=int_rsize) :: frac_bits
    ! Fraction bits of a subnormal input
  integer :: top_frac_bit
    ! Highest set fraction bit of a subnormal input
  real :: u
    ! Difference from 1 for near-one log1p path [nondim]

  xb = transfer(x, int_mold)

  ! Handle exceptional values before arithmetic range reduction.  This gives
  ! log(0) = -Inf with divide-by-zero, log(negative) = NaN with invalid,
  ! log(+Inf) = +Inf, and NaNs pass through with the usual signaling behavior.
  if (iand(xb, not(signmask)) == 0_int_rsize) then
    a = -1. / abs(x)
    return
  endif

  if (iand(xb, signmask) /= 0_int_rsize) then
    a = (x - x) / (x - x)
    return
  endif

  raw_exp = iand(ishft(xb, -expbit), expmax)

  if (raw_exp == expmax) then
    a = x + x
    return
  endif

  ! Avoid table-reduction cancellation for values very close to 1.  The
  ! x - 1 is exact throughout the near-one interval by Sterbenz's lemma.
  u = x - 1.
  if (xb >= near_lo_bits .and. xb <= near_hi_bits) then
    a = log1p_compensated_near(u)
    return
  endif

  ! Normalize subnormals before entering the shared table solver.
  exponent_adjustment = 0_int_rsize
  if (raw_exp == 0_int_rsize) then
    ! Normalize subnormal inputs using integer bit operations rather than
    ! arithmetic scaling, which can be defeated by flush-to-zero modes.
    frac_bits = iand(xb, expstep - 1_int_rsize)
    do top_frac_bit = expbit - 1, 0, -1
      if (btest(frac_bits, top_frac_bit)) exit
    enddo
    xb = ishft(int(expbias - expbit + top_frac_bit, int_rsize), expbit) + &
        ishft(frac_bits - ishft(1_int_rsize, top_frac_bit), expbit - top_frac_bit)
    exponent_adjustment = -Kbias
  endif

  a = log_table(xb, exponent_adjustment)
end function log_fallback


!> Shared table solver for positive normal or normalized subnormal inputs.
pure function log_table(xb, exponent_adjustment) result(a)
  integer(kind=int_rsize), intent(in) :: xb
    !< Bit representation of the positive normal input
  integer(kind=int_rsize), intent(in) :: exponent_adjustment
    !< Subnormal exponent correction; zero for normal inputs
  real :: a
    !< Natural logarithm of the original input [nondim]

  integer, parameter :: binshift = expbit - binbits
    !< Shift to move fraction subdivision index (idiv) to lowest bits
  integer(kind=int_rsize) :: K, mb, frac
    !< Binary exponent, normalized significand bits, and fraction bits
  integer :: idiv
    !< Lookup table subdivision index
  real :: z, r
    !< Reduced significand and reduced log1p argument [nondim]
  real :: w, hi, lo
    !< Compensated partial sums for log(x) [nondim]

  ! For encoding x = 2**K z, z ∈ [1,2),
  !   log(x) = K ln2 + log(z)

  ! Using z = 1 + r permits a Taylor-like Remez approximation of log(1 + r).
  ! This can be improved by decomposing [1,2) over N subintervals:
  !   z ∈ [ 1 + idiv/ndiv, 1 + (idiv + 1)/ndiv )

  ! Extract the fractional part of x, and use its top bits to determine idiv.
  frac = iand(xb, expstep - 1)
  idiv = int(shiftr(frac, binshift))

  ! Extract the exponent and significand from xb
  K = shiftr(xb, expbit) - expbias
  mb = ior(frac, transfer(1., int_mold))
  z = transfer(mb, real_mold)

  ! Apply the subnormal rescaling if present.
  K = K + exponent_adjustment

  ! Use idiv to rescale z to z' = z/c = 1 + r so that 1 + r is close to 1.
  ! After scaling,
  !   log(x) = K ln2 + log(c) + log(1 + r).

  ! NOTE: See log data table module for c(idiv) derivation.

  ! Compute r as ((z - c_hi) - c_lo) * (1/c) to avoid precision loss near c.
  r = ((z - c_hi(idiv)) - c_lo(idiv)) * inv_c(idiv)

  ! NOTE: The following can be more accurate, faster, and does not need c_hi
  !   and c_lo lookups, but must be evaluated by FMA, which cannot be ensured.
  ! r = z * inv_c(idiv) - 1

  ! Reconstruct log(x) with compensated high and low parts.

  ! Compute the upper part of log(x)
  w = real(K) * ln2_hi + log_c(idiv)
  hi = w + r

  ! XXX: Prevent (w - (w + r)) ~> r optimization by adding a subroundoff term.
  hi = hi + epsilon(real_mold) / 256.

  ! Recover the rounding residual of hi, and append K ln2_lo.
  ! Cancel w first, since it is the larger of w + r.
  lo = ((w - hi) + r) + real(K) * ln2_lo

  ! Apply the complete residual to hi to recover log(x).
  a = hi + (lo + log1p_tail_remez_4(r))
end function log_table


!> Remez estimate of log1p(x) - x over the log_repro() table-reduced range.
!! Uses a degree-4 polynomial for (log1p(x) - x) / x**2, giving a total
!! degree-6 approximation to log1p(x) - x.
pure function log1p_tail_remez_4(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-0.004, 0.004] [nondim]
  real :: a
    !< Approximation of log1p(x) - x [nondim]

  real, parameter :: c(0:4) = [ &
      -0.50000000000000003200073731434623651190475459492601993106505, &
       0.33333333332190456250105515593864671177817759161003655071438, &
      -0.24999999997399946327580391068339139441490813642200722246253, &
       0.20000285718021696112769866984902598077845774074847461743975, &
      -0.16667016672293213609527712180397636993208899588796181673368]
    !< Remez coefficients for (log1p(x) - x) / x**2 on [-0.004, 0.004] [nondim]
    !! Generated by generators/generate_MOM_log_remez.sollya.

  real :: x2, x4
    !< Powers of x [nondim]
  real :: p01, p23
    !< Grouped polynomial partial sums [nondim]

  x2 = x * x
  x4 = x2 * x2

  p01 = c(0) + x * c(1)
  p23 = c(2) + x * c(3)

  a = x2 * ((p01 + x2 * p23) + x4 * c(4))
end function log1p_tail_remez_4


!> Compensated near-one estimate of log1p(x).
pure function log1p_compensated_near(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-1/16, 0.065] [nondim]
  real :: a
    !< Approximation of log1p(x) [nondim]

  real, parameter :: split_scale = 134217728.
    !< Exact power-of-two scale used to split x into high and low parts [nondim]
  real, parameter :: c(3:13) = [ &
       0.3333333333333333495663409158244431039492812795260412782820753513055676037, &
      -0.2500000000000048410898553018121243145959173159979113851695191401132202101, &
       0.1999999999999401954600182759053457391449548132751979049871662291688667322, &
      -0.1666666666544826424331902956193816320602454785603791935540999385211449934, &
       0.142857142926292331379290143400260129042455919821774060063839433541603189, &
      -0.1250000111081452467282831140237236744260745509124580241506815300101707851, &
       0.1111110780004715236340979900688580293441238517810620127033976097471253938, &
      -9.999529690174345780528722225743902151690897998496979148198959705548315274e-2, &
       9.09152833243368478159650834402877235721130565569984927237483052376511638e-2, &
      -8.42729367861804930413103634851785932895197979844461832654488007987301762e-2, &
       7.684686137510498609385153891204008164960977636361018003750735920585730455e-2]
    !< Remez coefficients for log1p(x) - x + x*x/2 on [-1/16, 265/4096] [nondim]
  integer(kind=int_rsize), parameter :: split_mask = &
      -(2_int_rsize**((digits(real_mold) + 1) / 2))
    !< Retain at most half the significant bits, making x_hi*x_hi exact.

  real :: x2, x3
    !< Powers of the reduced argument [nondim]
  real :: x_hi, x_lo
    !< High and low parts of x [nondim]
  real :: w, hi, lo, tail
    !< Compensated partial sums and polynomial tail [nondim]
  real :: p1, p2, p3, p4
    !< Grouped polynomial partial sums [nondim]

  x2 = x * x
  x3 = x * x2

  p1 = (c(3) + (x * c(4))) + (x2 * c(5))
  p2 = (c(6) + (x * c(7))) + (x2 * c(8))
  p3 = (c(9) + (x * c(10))) + (x2 * c(11))
  p4 = c(12) + (x * c(13))
  tail = x3 * (p1 + (x3 * (p2 + (x3 * (p3 + (x3 * p4))))))

  ! Split x and compensate the dominant x - x*x/2 terms.  This follows the
  ! structure used by high-quality libm log1p kernels and reduces final-rounding
  ! error near x = 0, where the leading terms dominate the result.
  w = x * split_scale

  x_hi = transfer(iand(transfer(x, int_mold), split_mask), real_mold)
  x_lo = x - x_hi

  w = -0.5 * (x_hi * x_hi)
  hi = x + w
  ! Clear the lowest bit to prevent (x - (x + w)) => -w optimization
  hi = transfer(iand(transfer(hi, int_mold), not(1_int_rsize)), real_mold)
  lo = ((x - hi) + w) + ((-0.5 * x_lo) * (x_hi + x))

  a = (tail + lo) + hi
end function log1p_compensated_near

end submodule MOM_log
