! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Parabolic Method 1D reconstruction ported from "hybgen" module in Hycom, reproducing
!! the PPM_HYBGEN scheme
!!
!! The "_v1" indicates this is the first version. It extends recon1d_ppm_hybgen and differs only
!! in f() and average(), which retain the expressions used by the PPM_HYBGEN scheme (in
!! mom_remapping.average_value_ppoly()) so as to reproduce it bitwise. These take the form
!! based on the right edge value at x=1/2, and for intervals centered on x=1/2, so they are not
!! reflection symmetric. This implementation is reached with "C_PPM_HYBGEN_V1".
module Recon1d_PPM_hybgen_v1

use Recon1d_type, only : testing
use Recon1d_PPM_hybgen, only : PPM_hybgen

implicit none ; private

public PPM_hybgen_v1, testing

!> PPM reconstruction following "hybgen", reproducing the PPM_HYBGEN scheme
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_cw.init()
!! - reconstruct()          -> recon1d_ppm_hybgen.reconstruct()
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                 -> recon1d_ppm_cw.dfdx()
!! - check_reconstruction() -> recon1d_ppm_hybgen.check_reconstruction()
!! - unit_tests()           -> recon1d_ppm_hybgen.unit_tests()
!! - destroy()              -> recon1d_ppm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_ppm_cw.init()
!! - reconstruct_parent()   -> recon1d_ppm_cw.reconstruct()
type, extends (PPM_hybgen) :: PPM_hybgen_v1

contains
  !> Implementation of the PPM_hybgen_v1 average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PPM_hybgen_v1 reconstruction at a point [A]
  procedure :: f => f

end type PPM_hybgen_v1

contains

!> Average between xa and xb for cell k of a PPM_hybgen_v1 reconstruction [A]
!!
!! Calls unbounded_average(), or the point value f() when xa=xb, and then clamps the result to
!! [min(ul,ur), max(ul,ur)]. This replicates the evaluation and the force_bounds_in_subcell
!! behavior of the equivalent old-style PPM_HYBGEN scheme.
real function average(this, k, xa, xb)
  class(PPM_hybgen_v1), intent(in) :: this !< This reconstruction
  integer,              intent(in) :: k    !< Cell number
  real,                 intent(in) :: xa   !< Start of averaging interval on element (0 to 1) [nondim]
  real,                 intent(in) :: xb   !< End of averaging interval on element (0 to 1) [nondim]
  real :: u_lo, u_hi ! Bounds on the sub-cell average given by the edge values [A]

  if (xb > xa) then
    average = unbounded_average(this, k, xa, xb)
  else ! For a vanishing interval, use the point value. The expressions in f() are the same as
       ! those used for this case in the OM4-era average_value_ppoly().
    average = this%f(k, xa)
  endif
  u_lo = min(this%ul(k), this%ur(k))
  u_hi = max(this%ul(k), this%ur(k))
  average = max(u_lo, min(u_hi, average))

end function average

!> Average between xa and xb for cell k of a 1D PPM reconstruction, before bounding [A]
!!
!! This retains the form used by the PPM_HYBGEN scheme (in average_value_ppoly()), which
!! takes the form based on the right edge value for intervals centered on x=1/2.
real function unbounded_average(this, k, xa, xb)
  class(PPM_hybgen_v1), intent(in) :: this !< This reconstruction
  integer,              intent(in) :: k    !< Cell number
  real,                 intent(in) :: xa   !< Start of averaging interval on element (0 to 1) [nondim]
  real,                 intent(in) :: xb   !< End of averaging interval on element (0 to 1) [nondim]
  real :: xapxb                      ! A sum of fracional positions [nondim]
  real :: mx, Ya, Yb, my             ! Various fractional positions [nondim]
  real :: u_a, u_b                   ! Values at xa and xb [A]
  real :: xa2pxb2,  xa2b2ab, Ya2b2ab ! Sums of squared fractional positions [nondim]
  real :: a_L, a_R, u_c, a_c         ! Values of the polynomial at various locations [A]

  mx = 0.5 * ( xa + xb )
  a_L = this%ul(k)
  a_R = this%ur(k)
  u_c = this%u_mean(k)
  a_c = 0.5 * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6 / 6
  if (mx<0.5) then
    ! This integration of the PPM reconstruction is expressed in distances from the left edge
    xa2b2ab = (xa * xa + xb * xb) + xa * xb
    unbounded_average = a_L + ( ( a_R - a_L ) * mx &
                    + a_c * ( 3. * ( xb + xa ) - 2. * xa2b2ab ) )
  else
    ! This integration of the PPM reconstruction is expressed in distances from the right edge
    Ya = 1. - xa
    Yb = 1. - xb
    my = 0.5 * ( Ya + Yb )
    Ya2b2ab = (Ya * Ya + Yb * Yb) + Ya * Yb
    unbounded_average = a_R  + ( ( a_L - a_R ) * my &
                     + a_c * ( 3. * ( Yb + Ya ) - 2. * Ya2b2ab ) )
  endif

end function unbounded_average

!> Value of PPM_hybgen_v1 reconstruction at a point in cell k [A]
!!
!! This retains the form used by the PPM_HYBGEN scheme, which takes the form based on the
!! right edge value at x=1/2.
real function f(this, k, x)
  class(PPM_hybgen_v1), intent(in) :: this !< This reconstruction
  integer,              intent(in) :: k    !< Cell number
  real,                 intent(in) :: x    !< Non-dimensional position within element [nondim]
  real :: xc ! Bounded version of x [nondim]
  real :: du ! Difference across cell [A]
  real :: a6 ! Collela and Woordward curvature parameter [A]
  real :: u_a, u_b ! Two estimate of f [A]
  real :: lmx ! 1 - x [nondim]
  real :: wb ! Weight based on x [nondim]

  du = this%ur(k) - this%ul(k)
  a6 = 3.0 * ( ( this%u_mean(k) - this%ul(k) ) + ( this%u_mean(k) - this%ur(k) ) )
  xc = max( 0., min( 1., x ) )
  lmx = 1.0 - xc

  ! This expression for u_a can overshoot u_r but is good for x<<1
  u_a = this%ul(k) + xc * ( du + a6 * lmx )
  ! This expression for u_b can overshoot u_l but is good for 1-x<<1
  u_b = this%ur(k) + lmx * ( - du + a6 * xc )

  ! Since u_a and u_b are both side-bounded, using weights=0 or 1 will preserve uniformity
  wb = 0.5 + sign(0.5, xc - 0.5 ) ! = 1 @ x=0, = 0 @ x=1
  f = ( ( 1. - wb ) * u_a ) + ( wb * u_b )

end function f

!> \namespace recon1d_ppm_hybgen_v1
!!
!! \section section_recon1d_ppm_hybgen_v1 PPM_hybgen_v1 reconstruction
!!
!! Identical to recon1d_ppm_hybgen except that f() and average() reproduce the expressions of the
!! PPM_HYBGEN scheme.

end module Recon1d_PPM_hybgen_v1
