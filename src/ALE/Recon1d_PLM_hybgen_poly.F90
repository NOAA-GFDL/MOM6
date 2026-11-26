! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise Linear Method 1D reconstruction ported from "hybgen" module in Hycom, reproducing
!! the OM4-era evaluation of the reconstruction.
!!
!! This extends recon1d_plm_hybgen and differs only in average(), which uses the same expressions
!! as the OM4-era function mom_remapping::average_value_ppoly() so that the answers are bitwise
!! identical to REMAPPING_SCHEME="PLM_HYBGEN". In particular, for a vanishing interval (xa=xb)
!! the point value is evaluated from the nearer edge, which introduces an asymmetry that is not
!! present in recon1d_plm_hybgen.
module Recon1d_PLM_hybgen_poly

use Recon1d_PLM_hybgen, only : PLM_hybgen, testing

implicit none ; private

public PLM_hybgen_poly, testing

!> PLM reconstruction following "hybgen", reproducing the OM4-era evaluation.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_hybgen.init()
!! - reconstruct()          -> recon1d_plm_hybgen.reconstruct()
!! - average()                 *locally defined
!! - f()                    -> recon1d_plm_hybgen.f()
!! - dfdx()                 -> recon1d_plm_hybgen.dfdx()
!! - check_reconstruction() -> recon1d_plm_hybgen.check_reconstruction()
!! - unit_tests()           -> recon1d_plm_hybgen.unit_tests()
!! - destroy()              -> recon1d_plm_hybgen.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_plm_hybgen.init()
!! - reconstruct_parent()   -> recon1d_plm_hybgen.reconstruct()
type, extends (PLM_hybgen) :: PLM_hybgen_poly

contains
  !> Implementation of the PLM_hybgen_poly average over an interval [A]
  procedure :: average => average

end type PLM_hybgen_poly

contains

!> Average between xa and xb for cell k of a 1D PLM reconstruction [A]
!!
!! For xb>xa this is the same as recon1d_plm_hybgen.average(). For xa=xb the point value is
!! evaluated from the nearer edge, as in mom_remapping::average_value_ppoly().
real function average(this, k, xa, xb)
  class(PLM_hybgen_poly), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real,          intent(in) :: xa   !< Start of averaging interval on element (0 to 1) [nondim]
  real,          intent(in) :: xb   !< End of averaging interval on element (0 to 1) [nondim]

  if (xb > xa) then
    average = this%PLM_hybgen%average(k, xa, xb)
  else ! Point value, using the same expressions as the OM4-era average_value_ppoly()
    if (xa < 0.5) then
      average = this%ul(k) + xa * ( this%ur(k) - this%ul(k) )
    else
      average = this%ur(k) + ( 1. - xa ) * ( this%ul(k) - this%ur(k) )
    endif
  endif

end function average

!> \namespace recon1d_plm_hybgen_poly
!!
!! \section section_recon1d_plm_hybgen_poly PLM (hybgen) reproducing OM4-era answers
!!
!! This is recon1d_plm_hybgen with the evaluation of averages over vanishing intervals changed to
!! reproduce the OM4-era scheme REMAPPING_SCHEME="PLM_HYBGEN" bitwise. Use recon1d_plm_hybgen
!! for new work.

end module Recon1d_PLM_hybgen_poly
