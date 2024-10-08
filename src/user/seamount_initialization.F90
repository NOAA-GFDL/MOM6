!> Configures the model for the idealized seamount test case.
module seamount_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : sum_across_PEs
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mdl = "seamount_initialization" !< This module's name.

! The following routines are visible to the outside world
public seamount_initialize_topography
public seamount_initialize_thickness
public seamount_initialize_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Initialization of topography.
subroutine seamount_initialize_topography( D, G, param_file, max_depth )
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth [Z ~> m]

  ! Local variables
  real :: delta     ! Height of the seamount as a fraction of the maximum ocean depth [nondim]
  real :: x, y      ! Normalized positions relative to the domain center [nondim]
  real :: Lx, Ly    ! Seamount length scales normalized by the relevant domain sizes [nondim]
  real :: rLx, rLy  ! The Adcroft reciprocals of Lx and Ly [nondim]
  integer   :: i, j

  call get_param(param_file, mdl,"SEAMOUNT_DELTA", delta, &
                 "Non-dimensional height of seamount.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mdl,"SEAMOUNT_X_LENGTH_SCALE", Lx, &
                 "Length scale of seamount in x-direction. "//&
                 "Set to zero make topography uniform in the x-direction.", &
                 units=G%x_ax_unit_short, default=20.)
  call get_param(param_file, mdl,"SEAMOUNT_Y_LENGTH_SCALE", Ly, &
                 "Length scale of seamount in y-direction. "//&
                 "Set to zero make topography uniform in the y-direction.", &
                 units=G%y_ax_unit_short, default=0.)

  Lx = Lx / G%len_lon
  Ly = Ly / G%len_lat
  rLx = 0. ; if (Lx>0.) rLx = 1. / Lx
  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Compute normalized zonal coordinates (x,y=0 at center of domain)
    x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
    y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
    D(i,j) = G%max_depth * ( 1.0 - delta * exp(-((rLx*x)**2) - ((rLy*y)**2)) )
  enddo ; enddo

end subroutine seamount_initialize_topography

!> Initialization of thicknesses.
!! This subroutine initializes the layer thicknesses to be uniform.
subroutine seamount_initialize_thickness (h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  real :: e0(SZK_(GV)+1)  ! The resting interface heights [Z ~> m], usually
                          ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1) ! Interface height relative to the sea surface, positive upward [Z ~> m]
  real :: min_thickness   ! The minimum layer thicknesses [Z ~> m].
  real :: S_ref           ! A default value for salinities [S ~> ppt].
  real :: S_surf, S_range, S_light, S_dense ! Various salinities [S ~> ppt].
  real :: eta_IC_quanta   ! The granularity of quantization of intial interface heights [Z-1 ~> m-1].
  character(len=20) :: verticalCoordinate
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("seamount_initialization.F90, seamount_initialize_thickness: setting thickness")

  call get_param(param_file, mdl,"MIN_THICKNESS",min_thickness, &
                'Minimum thickness for layer', &
                 units='m', default=1.0e-3, do_not_log=just_read, scale=US%m_to_Z)
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)

  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)
  !do k=1,nz+1
  !  e0(k) = -G%max_depth * real(k-1) / real(nz)
  !enddo

  select case ( coordinateMode(verticalCoordinate) )

  case ( REGRIDDING_LAYER, REGRIDDING_RHO ) ! Initial thicknesses for isopycnal coordinates
    call get_param(param_file, mdl,"INITIAL_SSS", S_surf, &
                   units="ppt", default=34., scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl,"INITIAL_S_RANGE", S_range, &
                   units="ppt", default=2., scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "S_REF", S_ref, &
                   units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "TS_RANGE_S_LIGHT", S_light, &
                   units="ppt", default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "TS_RANGE_S_DENSE", S_dense, &
                   units="ppt", default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "INTERFACE_IC_QUANTA", eta_IC_quanta, &
                   "The granularity of initial interface height values "//&
                   "per meter, to avoid sensivity to order-of-arithmetic changes.", &
                   default=2048.0, units="m-1", scale=US%Z_to_m, do_not_log=just_read)
    if (just_read) return ! All run-time parameters have been read, so return.

    do K=1,nz+1
      ! Salinity of layer k is S_light + (k-1)/(nz-1) * (S_dense - S_light)
      ! Salinity of interface K is S_light + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Salinity at depth z should be S(z) = S_surf - S_range * z/max_depth
      ! Equating: S_surf - S_range * z/max_depth = S_light + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Equating: - S_range * z/max_depth = S_light - S_surf + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Equating: z/max_depth = - ( S_light - S_surf + (K-3/2)/(nz-1) * (S_dense - S_light) ) / S_range
      e0(K) = - G%max_depth * ( ( S_light  - S_surf ) + ( S_dense - S_light ) * &
                              ( (real(K)-1.5) / real(nz-1) ) ) / S_range
      ! Force round numbers ... the above expression has irrational factors ...
      if (eta_IC_quanta > 0.0) &
        e0(K) = nint(eta_IC_quanta*e0(K)) / eta_IC_quanta
      e0(K) = min(real(1-K)*GV%Angstrom_Z, e0(K)) ! Bound by surface
      e0(K) = max(-G%max_depth, e0(K)) ! Bound by bottom
    enddo
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_Z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_ZSTAR )                       ! Initial thicknesses for z coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = min_thickness
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      h(i,j,:) = depth_tot(i,j) / real(nz)
    enddo ; enddo

end select

end subroutine seamount_initialize_thickness

!> Initial values for temperature and salinity
subroutine seamount_initialize_temperature_salinity(T, S, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h !< Layer thickness [Z ~> m]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.

  ! Local variables
  real :: xi0, xi1  ! Fractional positions within the depth range [nondim]
  real :: r         ! A nondimensional sharpness parameter with an exponetial profile [nondim]
  real :: S_Ref     ! Default salinity range parameters [S ~> ppt].
  real :: T_Ref     ! Default temperature range parameters [C ~> degC].
  real :: S_Light, S_Dense, S_surf, S_range ! Salinity range parameters [S ~> ppt].
  real :: T_Light, T_Dense, T_surf, T_range ! Temperature range parameters [C ~> degC].
  real :: res_rat   ! The ratio of density space resolution in the denser part
                    ! of the range to that in the lighter part of the range.
                    ! Setting this greater than 1 increases the resolution for
                    ! the denser water [nondim].
  real :: a1, frac_dense, k_frac  ! Nondimensional temporary variables [nondim]
  integer :: i, j, k, is, ie, js, je, nz, k_light

  character(len=20) :: verticalCoordinate, density_profile

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl,"INITIAL_DENSITY_PROFILE", density_profile, &
                 'Initial profile shape. Valid values are "linear", "parabolic" '//&
                 'and "exponential".', default='linear', do_not_log=just_read)
  call get_param(param_file, mdl,"INITIAL_SSS", S_surf, &
                 'Initial surface salinity', &
                 units="ppt", default=34., scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl,"INITIAL_SST", T_surf, &
                 'Initial surface temperature', &
                 units="degC", default=0., scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl,"INITIAL_S_RANGE", S_range, &
                 'Initial salinity range (bottom - surface)', &
                 units="ppt", default=2., scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl,"INITIAL_T_RANGE", T_range, &
                 'Initial temperature range (bottom - surface)', &
                 units="degC", default=0., scale=US%degC_to_C, do_not_log=just_read)

  select case ( coordinateMode(verticalCoordinate) )
    case ( REGRIDDING_LAYER ) ! Initial thicknesses for layer isopycnal coordinates
      ! These parameters are used in MOM_fixed_initialization.F90 when CONFIG_COORD="ts_range"
      call get_param(param_file, mdl, "T_REF", T_ref, &
                 units="degC", default=10.0, scale=US%degC_to_C, do_not_log=.true.)
      call get_param(param_file, mdl, "TS_RANGE_T_LIGHT", T_light, &
                 units="degC", default=US%C_to_degC*T_Ref, scale=US%degC_to_C, do_not_log=.true.)
      call get_param(param_file, mdl, "TS_RANGE_T_DENSE", T_dense, &
                 units="degC", default=US%C_to_degC*T_Ref, scale=US%degC_to_C, do_not_log=.true.)
      call get_param(param_file, mdl, "S_REF", S_ref, &
                 units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=.true.)
      call get_param(param_file, mdl, "TS_RANGE_S_LIGHT", S_light, &
                 units="ppt", default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S, do_not_log=.true.)
      call get_param(param_file, mdl, "TS_RANGE_S_DENSE", S_dense, &
                 units="ppt", default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S, do_not_log=.true.)
      call get_param(param_file, mdl, "TS_RANGE_RESOLN_RATIO", res_rat, &
                 units="nondim", default=1.0, do_not_log=.true.)
      if (just_read) return ! All run-time parameters have been read, so return.

      ! Emulate the T,S used in the "ts_range" coordinate configuration code
      k_light = GV%nk_rho_varies + 1
      do j=js,je ; do i=is,ie
        T(i,j,k_light) = T_light ; S(i,j,k_light) = S_light
      enddo ; enddo
      a1 = 2.0 * res_rat / (1.0 + res_rat)
      do k=k_light+1,nz
        k_frac = real(k-k_light)/real(nz-k_light)
        frac_dense = a1 * k_frac + (1.0 - a1) * k_frac**2
        do j=js,je ; do i=is,ie
          T(i,j,k) = frac_dense * (T_Dense - T_Light) + T_Light
          S(i,j,k) = frac_dense * (S_Dense - S_Light) + S_Light
        enddo ; enddo
      enddo
    case ( REGRIDDING_SIGMA, REGRIDDING_ZSTAR, REGRIDDING_RHO ) ! All other coordinate use FV initialization
      if (just_read) return ! All run-time parameters have been read, so return.
      do j=js,je ; do i=is,ie
        xi0 = 0.0
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth
          select case ( trim(density_profile) )
            case ('linear')
             !S(i,j,k) = S_surf + S_range * 0.5 * (xi0 + xi1)
              S(i,j,k) = S_surf + ( 0.5 * S_range ) * (xi0 + xi1) ! Coded this way to reproduce old hard-coded answers
              T(i,j,k) = T_surf + T_range * 0.5 * (xi0 + xi1)
            case ('parabolic')
              S(i,j,k) = S_surf + S_range * (2.0 / 3.0) * (xi1**3 - xi0**3) / (xi1 - xi0)
              T(i,j,k) = T_surf + T_range * (2.0 / 3.0) * (xi1**3 - xi0**3) / (xi1 - xi0)
            case ('exponential')
              r = 0.8 ! small values give sharp profiles
              S(i,j,k) = S_surf + S_range * (exp(xi1/r)-exp(xi0/r)) / (xi1 - xi0)
              T(i,j,k) = T_surf + T_range * (exp(xi1/r)-exp(xi0/r)) / (xi1 - xi0)
            case default
              call MOM_error(FATAL, 'Unknown value for "INITIAL_DENSITY_PROFILE"')
          end select
          xi0 = xi1
        enddo
      enddo ; enddo
  end select

end subroutine seamount_initialize_temperature_salinity

end module seamount_initialization
