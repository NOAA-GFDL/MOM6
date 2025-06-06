!> Provides a transparent vertical ocean grid type and supporting routines
module MOM_verticalGrid

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public verticalGridInit, verticalGridEnd
public setVerticalGridAxes
public get_flux_units, get_thickness_units, get_tr_flux_units

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Describes the vertical ocean grid, including unit conversion factors
type, public :: verticalGrid_type

  ! Commonly used parameters
  integer :: ke     !< The number of layers/levels in the vertical
  real :: max_depth !< The maximum depth of the ocean [Z ~> m].
!  real :: mks_g_Earth !< The gravitational acceleration in unscaled MKS units [m s-2].  This might not be used.
  real :: g_Earth   !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: g_Earth_Z_T2 !< The gravitational acceleration in alternatively rescaled units [Z T-2 ~> m s-2]
  real :: Rho0      !< The density used in the Boussinesq approximation or nominal
                    !! density used to convert depths into mass units [R ~> kg m-3].

  ! Vertical coordinate descriptions for diagnostics and I/O
  character(len=40) :: zAxisUnits !< The units that vertical coordinates are written in
  character(len=40) :: zAxisLongName !< Coordinate name to appear in files,
                                  !! e.g. "Target Potential Density" or "Height"
  real, allocatable, dimension(:) :: sLayer !< Coordinate values of layer centers, in unscaled
                        !! units that depend on the vertical coordinate, such as [kg m-3] for an
                        !! isopycnal or some hybrid coordinates, [m] for a Z* coordinate,
                        !! or [nondim] for a sigma coordinate.
  real, allocatable, dimension(:) :: sInterface !< Coordinate values on interfaces, in the same
                        !! unscale units as sLayer [various].
  integer :: direction = 1 !< Direction defaults to 1, positive up.

  ! The following variables give information about the vertical grid.
  logical :: Boussinesq !< If true, make the Boussinesq approximation.
  logical :: semi_Boussinesq !< If true, do non-Boussinesq pressure force calculations and
                        !! use mass-based "thicknesses, but use Rho0 to convert layer thicknesses
                        !! into certain height changes.  This only applies if BOUSSINESQ is false.
  real :: Angstrom_H    !< A one-Angstrom thickness in the model thickness units [H ~> m or kg m-2].
  real :: Angstrom_Z    !< A one-Angstrom thickness in the model depth units [Z ~> m].
  real :: Angstrom_m    !< A one-Angstrom thickness [m].
  real :: H_subroundoff !< A thickness that is so small that it can be added to a thickness of
                        !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
                        !! If Angstrom is 0 or exceedingly small, this is negligible compared to 1e-17 m.
  real :: dZ_subroundoff !< A thickness in height units that is so small that it can be added to a
                        !! vertical distance of Angstrom_Z or 1e-17 m without changing it at the bit
                        !! level [Z ~> m].  This is the height equivalent of H_subroundoff.
  real, allocatable, dimension(:) :: &
    g_prime, &          !< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
    Rlay                !< The target coordinate value (potential density) in each layer [R ~> kg m-3].
  integer :: nkml = 0   !< The number of layers at the top that should be treated
                        !! as parts of a homogeneous region.
  integer :: nk_rho_varies = 0 !< The number of layers at the top where the
                        !! density does not track any target density.
  real :: H_to_kg_m2    !< A constant that translates thicknesses from the units of thickness
                        !! to kg m-2 [kg m-2 H-1 ~> kg m-3 or 1].
  real :: kg_m2_to_H    !< A constant that translates thicknesses from kg m-2 to the units
                        !! of thickness [H m2 kg-1 ~> m3 kg-1 or 1].
  real :: m_to_H        !< A constant that translates distances in m to the units of
                        !! thickness [H m-1 ~> 1 or kg m-3].
  real :: H_to_m        !< A constant that translates distances in the units of thickness
                        !! to m [m H-1 ~> 1 or m3 kg-1].
  real :: H_to_Pa       !< A constant that translates the units of thickness to pressure
                        !! [Pa H-1 = kg m-1 s-2 H-1 ~> kg m-2 s-2 or m s-2].
  real :: H_to_Z        !< A constant that translates thickness units to the units of
                        !! depth [Z H-1 ~> 1 or m3 kg-1].
  real :: Z_to_H        !< A constant that translates depth units to thickness units
                        !! depth [H Z-1 ~> 1 or kg m-3].
  real :: H_to_RZ       !< A constant that translates thickness units to the units of
                        !! mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real :: RZ_to_H       !< A constant that translates mass per unit area units to
                        !! thickness units [H R-1 Z-1 ~> m3 kg-2 or 1].
  real :: H_to_MKS      !< A constant that translates thickness units to its MKS unit
                        !! (m or kg m-2) based on GV%Boussinesq [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
  real :: m2_s_to_HZ_T  !< The combination of conversion factors that converts kinematic viscosities
                        !! in m2 s-1 to the internal units of the kinematic (in Boussinesq mode)
                        !! or dynamic viscosity [H Z s T-1 m-2 ~> 1 or kg m-3]
  real :: HZ_T_to_m2_s  !< The combination of conversion factors that converts the viscosities from
                        !! their internal representation into a kinematic viscosity in m2 s-1
                        !! [T m2 H-1 Z-1 s-1 ~> 1 or m3 kg-1]
  real :: HZ_T_to_MKS   !< The combination of conversion factors that converts the viscosities from
                        !! their internal representation into their unnscaled MKS units
                        !! (m2 s-1 or Pa s), depending on whether the model is Boussinesq
                        !! [T m2 H-1 Z-1 s-1 ~> 1] or [T Pa s H-1 Z-1 ~> 1]

end type verticalGrid_type

contains

!> Allocates and initializes the ocean model vertical grid structure.
subroutine verticalGridInit( param_file, GV, US )
  type(param_file_type),   intent(in) :: param_file !< Parameter file handle/type
  type(verticalGrid_type), pointer    :: GV         !< The container for vertical grid data
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type
  ! This routine initializes the verticalGrid_type structure (GV).
  ! All memory is allocated but not necessarily set to meaningful values until later.

  ! Local variables
  integer :: nk, H_power
  real    :: H_rescale_factor ! The integer power of 2 by which thicknesses are rescaled [nondim]
  real    :: rho_Kv  ! The density used convert input kinematic viscosities into dynamic viscosities
                     ! when in non-Boussinesq mode [R ~> kg m-3]
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=16) :: mdl = 'MOM_verticalGrid'

  if (associated(GV)) call MOM_error(FATAL, &
     'verticalGridInit: called with an associated GV pointer.')
  allocate(GV)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
                   "Parameters providing information about the vertical grid.", &
                   log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "G_EARTH", GV%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%Z_to_m*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "RHO_0", GV%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "BOUSSINESQ", GV%Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true.)
  call get_param(param_file, mdl, "SEMI_BOUSSINESQ", GV%semi_Boussinesq, &
                 "If true, do non-Boussinesq pressure force calculations and use mass-based "//&
                 "thicknesses, but use RHO_0 to convert layer thicknesses into certain "//&
                 "height changes.  This only applies if BOUSSINESQ is false.", &
                 default=.true., do_not_log=GV%Boussinesq)
  if (GV%Boussinesq) GV%semi_Boussinesq = .true.
  call get_param(param_file, mdl, "RHO_KV_CONVERT", Rho_Kv, &
                 "The density used to convert input vertical distances into thickesses in "//&
                 "non-BOUSSINESQ mode, and to convert kinematic viscosities into dynamic "//&
                 "viscosities and similarly for vertical diffusivities.  GV%m_to_H is set "//&
                 "using this value, whereas GV%Z_to_H is set using RHO_0.  The default is "//&
                 "RHO_0, but this can be set separately to demonstrate the independence of the "//&
                 "non-Boussinesq solutions of the value of RHO_0.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=GV%Boussinesq)
  call get_param(param_file, mdl, "ANGSTROM", GV%Angstrom_Z, &
                 "The minimum layer thickness, usually one-Angstrom.", &
                 units="m", default=1.0e-10, scale=US%m_to_Z)
  call get_param(param_file, mdl, "H_RESCALE_POWER", H_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of thickness.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  if (abs(H_power) > 300) call MOM_error(FATAL, "verticalGridInit: "//&
                 "H_RESCALE_POWER is outside of the valid range of -300 to 300.")
  H_rescale_factor = 1.0
  if (H_power /= 0) H_rescale_factor = 2.0**H_power
  if (.not.GV%Boussinesq) then
    call get_param(param_file, mdl, "H_TO_KG_M2", GV%H_to_kg_m2,&
                 "A constant that translates thicknesses from the model's "//&
                 "internal units of thickness to kg m-2.", units="kg m-2 H-1", &
                 default=1.0)
    GV%H_to_kg_m2 = GV%H_to_kg_m2 * H_rescale_factor
  else
    call get_param(param_file, mdl, "H_TO_M", GV%H_to_m, &
                 "A constant that translates the model's internal "//&
                 "units of thickness into m.", units="m H-1", default=1.0)
    GV%H_to_m = GV%H_to_m * H_rescale_factor
  endif
  ! This is not used:  GV%mks_g_Earth = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth
  GV%g_Earth_Z_T2 = US%L_to_Z**2 * GV%g_Earth  ! This would result from scale=US%m_to_Z*US%T_to_s**2.
#ifdef STATIC_MEMORY_
  ! Here NK_ is a macro, while nk is a variable.
  call get_param(param_file, mdl, "NK", nk, &
                 "The number of model layers.", units="nondim", &
                 default=NK_)
  if (nk /= NK_) call MOM_error(FATAL, "verticalGridInit: " // &
       "Mismatched number of layers NK_ between MOM_memory.h and param_file")

#else
  call get_param(param_file, mdl, "NK", nk, &
                 "The number of model layers.", units="nondim", fail_if_missing=.true.)
#endif
  GV%ke = nk

  if (GV%Boussinesq) then
    GV%H_to_kg_m2 = US%R_to_kg_m3*GV%Rho0 * GV%H_to_m
    GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
    GV%m_to_H = 1.0 / GV%H_to_m
    GV%H_to_MKS = GV%H_to_m
    GV%m2_s_to_HZ_T = GV%m_to_H * US%m_to_Z * US%T_to_s

    GV%H_to_Z = GV%H_to_m * US%m_to_Z
    GV%Z_to_H = US%Z_to_m * GV%m_to_H
  else
    GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
    !  GV%m_to_H = US%R_to_kg_m3*GV%Rho0 * GV%kg_m2_to_H
    GV%m_to_H = US%R_to_kg_m3*rho_Kv * GV%kg_m2_to_H
    GV%H_to_MKS = GV%H_to_kg_m2
    GV%m2_s_to_HZ_T = US%R_to_kg_m3*rho_Kv * GV%kg_m2_to_H * US%m_to_Z * US%T_to_s
    GV%H_to_m = 1.0 / GV%m_to_H

    GV%H_to_Z = US%m_to_Z * ( GV%H_to_kg_m2 / (US%R_to_kg_m3*GV%Rho0) )
    GV%Z_to_H = US%Z_to_m * ( US%R_to_kg_m3*GV%Rho0 * GV%kg_m2_to_H )
  endif

  GV%Angstrom_H = (US%Z_to_m * GV%m_to_H) * GV%Angstrom_Z
  GV%Angstrom_m = US%Z_to_m * GV%Angstrom_Z

  GV%H_subroundoff = 1e-20 * max(GV%Angstrom_H, GV%m_to_H*1e-17)
  GV%dZ_subroundoff = 1e-20 * max(GV%Angstrom_Z, US%m_to_Z*1e-17)

  GV%H_to_Pa = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth * GV%H_to_kg_m2

  GV%H_to_RZ = GV%H_to_kg_m2 * US%kg_m3_to_R * US%m_to_Z
  GV%RZ_to_H = GV%kg_m2_to_H * US%R_to_kg_m3 * US%Z_to_m

  GV%HZ_T_to_m2_s = 1.0 / GV%m2_s_to_HZ_T
  GV%HZ_T_to_MKS = GV%H_to_MKS * US%Z_to_m * US%s_to_T

  ! Note based on the above that for both Boussinsq and non-Boussinesq cases that:
  !     GV%Rho0 = GV%Z_to_H * GV%H_to_RZ
  !     1.0/GV%Rho0 = GV%H_to_Z * GV%RZ_to_H
  ! This is exact for power-of-2 scaling of the units, regardless of the value of Rho0, but
  ! the first term on the right hand side is invertable in Boussinesq mode, but the second
  ! is invertable when non-Boussinesq.

  ! Log derivative values.
  call log_param(param_file, mdl, "M to THICKNESS", GV%m_to_H*H_rescale_factor, units="H m-1")
  call log_param(param_file, mdl, "M to THICKNESS rescaled by 2^-n", GV%m_to_H, units="2^n H m-1")
  call log_param(param_file, mdl, "THICKNESS to M rescaled by 2^n", GV%H_to_m, units="2^-n m H-1")

  allocate( GV%sInterface(nk+1) )
  allocate( GV%sLayer(nk) )
  allocate( GV%g_prime(nk+1), source=0.0 )
  allocate( GV%Rlay(nk), source=0.0 )

end subroutine verticalGridInit

!> Returns the model's thickness units, usually m or kg/m^2.
function get_thickness_units(GV)
  character(len=48)                   :: get_thickness_units !< The vertical thickness units
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  !   This subroutine returns the appropriate units for thicknesses,
  ! depending on whether the model is Boussinesq or not and the scaling for
  ! the vertical thickness.

  if (GV%Boussinesq) then
    get_thickness_units = "m"
  else
    get_thickness_units = "kg m-2"
  endif
end function get_thickness_units

!> Returns the model's thickness flux units, usually m^3/s or kg/s.
function get_flux_units(GV)
  character(len=48)                   :: get_flux_units !< The thickness flux units
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  !   This subroutine returns the appropriate units for thickness fluxes,
  ! depending on whether the model is Boussinesq or not and the scaling for
  ! the vertical thickness.

  if (GV%Boussinesq) then
    get_flux_units = "m3 s-1"
  else
    get_flux_units = "kg s-1"
  endif
end function get_flux_units

!> Returns the model's tracer flux units.
function get_tr_flux_units(GV, tr_units, tr_vol_conc_units, tr_mass_conc_units)
  character(len=48)                      :: get_tr_flux_units !< The model's flux units
                                                              !! for a tracer.
  type(verticalGrid_type),    intent(in) :: GV                !< The ocean's vertical
                                                              !! grid structure.
  character(len=*), optional, intent(in) :: tr_units          !< Units for a tracer, for example
                                                              !! Celsius or PSU.
  character(len=*), optional, intent(in) :: tr_vol_conc_units !< The concentration units per unit
                                                              !! volume, for example if the units are
                                                              !! umol m-3, tr_vol_conc_units would
                                                              !! be umol.
  character(len=*), optional, intent(in) :: tr_mass_conc_units !< The concentration units per unit
                                                              !! mass of sea water, for example if
                                                              !! the units are mol kg-1,
                                                              !! tr_vol_conc_units would be mol.

  !   This subroutine returns the appropriate units for thicknesses and fluxes,
  ! depending on whether the model is Boussinesq or not and the scaling for
  ! the vertical thickness.
  integer :: cnt

  cnt = 0
  if (present(tr_units)) cnt = cnt+1
  if (present(tr_vol_conc_units)) cnt = cnt+1
  if (present(tr_mass_conc_units)) cnt = cnt+1

  if (cnt == 0) call MOM_error(FATAL, "get_tr_flux_units: One of the three "//&
    "arguments tr_units, tr_vol_conc_units, or tr_mass_conc_units "//&
    "must be present.")
  if (cnt > 1) call MOM_error(FATAL, "get_tr_flux_units: Only one of "//&
    "tr_units, tr_vol_conc_units, and tr_mass_conc_units may be present.")
  if (present(tr_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_units)//" m3 s-1"
    else
      get_tr_flux_units = trim(tr_units)//" kg s-1"
    endif
  endif
  if (present(tr_vol_conc_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_vol_conc_units)//" s-1"
    else
      get_tr_flux_units = trim(tr_vol_conc_units)//" m-3 kg s-1"
    endif
  endif
  if (present(tr_mass_conc_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_mass_conc_units)//" kg-1 m3 s-1"
    else
      get_tr_flux_units = trim(tr_mass_conc_units)//" s-1"
    endif
  endif

end function get_tr_flux_units

!> This sets the coordinate data for the "layer mode" of the isopycnal model.
subroutine setVerticalGridAxes( Rlay, GV, scale )
  type(verticalGrid_type), intent(inout) :: GV    !< The container for vertical grid data
  real, dimension(GV%ke),  intent(in)    :: Rlay  !< The layer target density [R ~> kg m-3]
  real,                    intent(in)    :: scale !< A unit scaling factor for Rlay to convert
                                                  !! it into the units of sInterface, usually
                                                  !! [kg m-3 R-1 ~> 1] when used in layer mode.
  ! Local variables
  integer :: k, nk

  nk = GV%ke

  GV%zAxisLongName = 'Target Potential Density'
  GV%zAxisUnits = 'kg m-3'
  do k=1,nk ; GV%sLayer(k) = scale*Rlay(k) ; enddo
  if (nk > 1) then
    GV%sInterface(1) = scale * (1.5*Rlay(1) - 0.5*Rlay(2))
    do K=2,nk ; GV%sInterface(K) = scale * 0.5*( Rlay(k-1) + Rlay(k) ) ; enddo
    GV%sInterface(nk+1) = scale * (1.5*Rlay(nk) - 0.5*Rlay(nk-1))
  else
    GV%sInterface(1) = 0.0 ; GV%sInterface(nk+1) = 2.0*scale*Rlay(nk)
  endif

end subroutine setVerticalGridAxes

!> Deallocates the model's vertical grid structure.
subroutine verticalGridEnd( GV )
  type(verticalGrid_type), pointer :: GV !< The ocean's vertical grid structure

  deallocate( GV%g_prime, GV%Rlay )
  deallocate( GV%sInterface , GV%sLayer )
  deallocate( GV )

end subroutine verticalGridEnd

end module MOM_verticalGrid
