!> Interface to CVMix double diffusion scheme.
module MOM_CVMix_ddiff

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,  only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,  only : post_data
use MOM_EOS,            only : calculate_density_derivs
use MOM_error_handler,  only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_debugging,      only : hchksum
use MOM_grid,           only : ocean_grid_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use cvmix_ddiff,        only : cvmix_init_ddiff, CVMix_coeffs_ddiff
use cvmix_kpp,          only : CVmix_kpp_compute_kOBL_depth
implicit none ; private

#include <MOM_memory.h>

public CVMix_ddiff_init, CVMix_ddiff_end, CVMix_ddiff_is_used, compute_ddiff_coeffs

!> Control structure including parameters for CVMix double diffusion.
type, public :: CVMix_ddiff_cs ; private

  ! Parameters
  real    :: strat_param_max !< maximum value for the stratification parameter [nondim]
  real    :: kappa_ddiff_s   !< leading coefficient in formula for salt-fingering regime
                             !! for salinity diffusion [Z2 T-1 ~> m2 s-1]
  real    :: ddiff_exp1      !< interior exponent in salt-fingering regime formula [nondim]
  real    :: ddiff_exp2      !< exterior exponent in salt-fingering regime formula [nondim]
  real    :: mol_diff        !< molecular diffusivity [Z2 T-1 ~> m2 s-1]
  real    :: kappa_ddiff_param1 !< exterior coefficient in diffusive convection regime [nondim]
  real    :: kappa_ddiff_param2 !< middle coefficient in diffusive convection regime [nondim]
  real    :: kappa_ddiff_param3 !< interior coefficient in diffusive convection regime [nondim]
  real    :: min_thickness      !< Minimum thickness allowed [H ~> m or kg m-2]
  character(len=4) :: diff_conv_type !< type of diffusive convection to use. Options are Marmorino &
                                !! Caldwell 1976 ("MC76"; default) and Kelley 1988, 1990 ("K90")
  logical :: debug              !< If true, turn on debugging

end type CVMix_ddiff_cs

character(len=40)  :: mdl = "MOM_CVMix_ddiff"     !< This module's name.

contains

!> Initialized the CVMix double diffusion module.
logical function CVMix_ddiff_init(Time, G, GV, US, param_file, diag, CS)

  type(time_type),         intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),   intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(CVMix_ddiff_cs),    pointer       :: CS         !< This module's control structure.

  ! This include declares and sets the variable "version".
# include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(WARNING, "CVMix_ddiff_init called with an associated "// &
                            "control structure.")
    return
  endif

  ! Read parameters
  call get_param(param_file, mdl, "USE_CVMIX_DDIFF", CVMix_ddiff_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
           "Parameterization of mixing due to double diffusion processes via CVMix", &
           all_default=.not.CVMix_ddiff_init)
  call get_param(param_file, mdl, "USE_CVMIX_DDIFF", CVMix_ddiff_init, &
                 "If true, turns on double diffusive processes via CVMix. "//&
                 "Note that double diffusive processes on viscosity are ignored "//&
                 "in CVMix, see http://cvmix.github.io/ for justification.", &
                 default=.false.)

  if (.not. CVMix_ddiff_init) return
  allocate(CS)

  call get_param(param_file, mdl, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

  call get_param(param_file, mdl, 'MIN_THICKNESS', CS%min_thickness, &
                 units="m", scale=GV%m_to_H, default=0.001, do_not_log=.True.)

  call openParameterBlock(param_file,'CVMIX_DDIFF')

  call get_param(param_file, mdl, "STRAT_PARAM_MAX", CS%strat_param_max, &
                 "The maximum value for the double dissusion stratification parameter", &
                 units="nondim", default=2.55)

  call get_param(param_file, mdl, "KAPPA_DDIFF_S", CS%kappa_ddiff_s, &
                 "Leading coefficient in formula for salt-fingering regime for salinity diffusion.", &
                 units="m2 s-1", default=1.0e-4, scale=US%m2_s_to_Z2_T)

  call get_param(param_file, mdl, "DDIFF_EXP1", CS%ddiff_exp1, &
                 "Interior exponent in salt-fingering regime formula.", &
                 units="nondim", default=1.0)

  call get_param(param_file, mdl, "DDIFF_EXP2", CS%ddiff_exp2, &
                 "Exterior exponent in salt-fingering regime formula.", &
                 units="nondim", default=3.0)

  call get_param(param_file, mdl, "KAPPA_DDIFF_PARAM1", CS%kappa_ddiff_param1, &
                "Exterior coefficient in diffusive convection regime.", &
                 units="nondim", default=0.909)

  call get_param(param_file, mdl, "KAPPA_DDIFF_PARAM2", CS%kappa_ddiff_param2, &
                "Middle coefficient in diffusive convection regime.", &
                 units="nondim", default=4.6)

  call get_param(param_file, mdl, "KAPPA_DDIFF_PARAM3", CS%kappa_ddiff_param3, &
                "Interior coefficient in diffusive convection regime.", &
                 units="nondim", default=-0.54)

  call get_param(param_file, mdl, "MOL_DIFF", CS%mol_diff, &
                 "Molecular diffusivity used in CVMix double diffusion.", &
                 units="m2 s-1", default=1.5e-6,  scale=US%m2_s_to_Z2_T)

  call get_param(param_file, mdl, "DIFF_CONV_TYPE", CS%diff_conv_type, &
                 "type of diffusive convection to use. Options are Marmorino \n" //&
                 "and Caldwell 1976 (MC76) and Kelley 1988, 1990 (K90).", &
                 default="MC76")

  call closeParameterBlock(param_file)

  call cvmix_init_ddiff(strat_param_max=CS%strat_param_max,          &
                        kappa_ddiff_s=US%Z2_T_to_m2_s*CS%kappa_ddiff_s, &
                        ddiff_exp1=CS%ddiff_exp1,                 &
                        ddiff_exp2=CS%ddiff_exp2,                 &
                        mol_diff=US%Z2_T_to_m2_s*CS%mol_diff,     &
                        kappa_ddiff_param1=CS%kappa_ddiff_param1, &
                        kappa_ddiff_param2=CS%kappa_ddiff_param2, &
                        kappa_ddiff_param3=CS%kappa_ddiff_param3, &
                        diff_conv_type=CS%diff_conv_type)

end function CVMix_ddiff_init

!> Subroutine for computing vertical diffusion coefficients for the
!! double diffusion mixing parameterization.
subroutine compute_ddiff_coeffs(h, tv, G, GV, US, j, Kd_T, Kd_S, CS, R_rho)

  type(ocean_grid_type),                      intent(in)    :: G    !< Grid structure.
  type(verticalGrid_type),                    intent(in)    :: GV   !< Vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h    !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                      intent(in)    :: tv   !< Thermodynamics structure.
  type(unit_scale_type),                      intent(in)    :: US   !< A dimensional unit scaling type
  integer,                                    intent(in)    :: j    !< Meridional grid index to work on.
  ! Kd_T and Kd_S are intent inout because only one j-row is set here, but they are essentially outputs.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kd_T !< Interface double diffusion diapycnal
                                                                    !! diffusivity for temperature
                                                                    !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kd_S !< Interface double diffusion diapycnal
                                                                    !! diffusivity for salinity
                                                                    !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  type(CVMix_ddiff_cs),                       pointer       :: CS   !< The control structure returned
                                                                    !! by a previous call to CVMix_ddiff_init.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                    optional, intent(inout) :: R_rho !< The density ratios at interfaces [nondim].

  ! Local variables
  real, dimension(SZK_(GV)) :: &
    cellHeight, &  !< Height of cell centers relative to the sea surface [H ~> m or kg m-2]
    dRho_dT,    &  !< partial derivatives of density with temperature [R C-1 ~> kg m-3 degC-1]
    dRho_dS,    &  !< partial derivatives of density with salinity [R S-1 ~> kg m-3 ppt-1]
    pres_int,   &  !< pressure at each interface [R L2 T-2 ~> Pa]
    temp_int,   &  !< temp and at interfaces [C ~> degC]
    salt_int,   &  !< salt at at interfaces [S ~> ppt]
    alpha_dT,   &  !< alpha*dT across interfaces [kg m-3]
    beta_dS,    &  !< beta*dS across interfaces [kg m-3]
    dT,         &  !< temperature difference between adjacent layers [C ~> degC]
    dS             !< salinity difference between adjacent layers [S ~> ppt]
  real, dimension(SZK_(GV)+1) :: &
    Kd1_T,      &  !< Diapycanal diffusivity of temperature [m2 s-1].
    Kd1_S          !< Diapycanal diffusivity of salinity [m2 s-1].

  real, dimension(SZK_(GV)+1) :: iFaceHeight !< Height of interfaces relative to the sea surface [H ~> m or kg m-2]
  real :: dh, hcorr ! Limited thicknesses and a cumulative correction [H ~> m or kg m-2]
  integer :: i, k

  ! initialize dummy variables
  pres_int(:) = 0.0; temp_int(:) = 0.0; salt_int(:) = 0.0
  alpha_dT(:) = 0.0; beta_dS(:) = 0.0; dRho_dT(:) = 0.0
  dRho_dS(:) = 0.0; dT(:) = 0.0; dS(:) = 0.0


  ! GMM, I am leaving some code commented below. We need to pass BLD to
  ! this subroutine to avoid adding diffusivity above that. This needs
  ! to be done once we re-structure the order of the calls.
  !if (.not. associated(hbl)) then
  !  allocate(hbl(SZI_(G), SZJ_(G)));
  !  hbl(:,:) = 0.0
  !endif

  do i = G%isc, G%iec

    ! skip calling at land points
    if (G%mask2dT(i,j) == 0.) cycle

    pres_int(1) = 0. ;  if (associated(tv%p_surf)) pres_int(1) = tv%p_surf(i,j)
    ! we don't have SST and SSS, so let's use values at top-most layer
    temp_int(1) = tv%T(i,j,1); salt_int(1) = tv%S(i,j,1)
    do K=2,GV%ke
      ! pressure at interface
      pres_int(K) = pres_int(K-1) + (GV%g_Earth * GV%H_to_RZ) * h(i,j,k-1)
      ! temp and salt at interface
      ! for temp: (t1*h1 + t2*h2)/(h1+h2)
      temp_int(K) = (tv%T(i,j,k-1)*h(i,j,k-1) + tv%T(i,j,k)*h(i,j,k)) / (h(i,j,k-1)+h(i,j,k))
      salt_int(K) = (tv%S(i,j,k-1)*h(i,j,k-1) + tv%S(i,j,k)*h(i,j,k)) / (h(i,j,k-1)+h(i,j,k))
      ! dT and dS
      dT(K) = (tv%T(i,j,k-1)-tv%T(i,j,k))
      dS(K) = (tv%S(i,j,k-1)-tv%S(i,j,k))
    enddo ! k-loop finishes

    call calculate_density_derivs(temp_int, salt_int, pres_int, drho_dT, drho_dS, tv%eqn_of_state)

    ! The "-1.0" below is needed so that the following criteria is satisfied:
    ! if ((alpha_dT > beta_dS) .and. (beta_dS > 0.0)) then "salt finger"
    ! if ((alpha_dT < 0.) .and. (beta_dS < 0.) .and. (alpha_dT > beta_dS)) then "diffusive convection"
    do k=1,GV%ke
      alpha_dT(k) = -1.0*US%R_to_kg_m3*drho_dT(k) * dT(k)
      beta_dS(k)  = US%R_to_kg_m3*drho_dS(k) * dS(k)
    enddo

    if (present(R_rho))  then
      do k=1,GV%ke
        ! Set R_rho using Adcroft's rule of reciprocals.
        R_rho(i,j,k) = 0.0 ; if (abs(beta_dS(k)) > 0.0) R_rho(i,j,k) = alpha_dT(k) / beta_dS(k)
        ! avoid NaN's again for safety, perhaps unnecessarily.
        if (R_rho(i,j,k) /= R_rho(i,j,k)) R_rho(i,j,k) = 0.0
      enddo
    endif

    iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
    hcorr = 0.0
    ! compute heights at cell center and interfaces
    do k=1,GV%ke
      dh = h(i,j,k) ! Nominal thickness to use for increment, in height units
      dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
      hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
      dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
      cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
      iFaceHeight(k+1) = iFaceHeight(k) - dh
    enddo

    ! gets index of the level and interface above hbl in [H ~> m or kg m-2]
    !kOBL = CVmix_kpp_compute_kOBL_depth(iFaceHeight, cellHeight, hbl(i,j))

    Kd1_T(:) = 0.0 ; Kd1_S(:) = 0.0
    call CVMix_coeffs_ddiff(Tdiff_out=Kd1_T(:), &
                            Sdiff_out=Kd1_S(:), &
                            strat_param_num=alpha_dT(:), &
                            strat_param_denom=beta_dS(:), &
                            nlev=GV%ke,    &
                            max_nlev=GV%ke)
    do K=1,GV%ke+1
      Kd_T(i,j,K) = GV%m2_s_to_HZ_T * Kd1_T(K)
      Kd_S(i,j,K) = GV%m2_s_to_HZ_T * Kd1_S(K)
    enddo

    ! Do not apply mixing due to convection within the boundary layer
    !do k=1,kOBL
    !  Kd_T(i,j,k) = 0.0
    !  Kd_S(i,j,k) = 0.0
    !enddo

  enddo ! i-loop

end subroutine compute_ddiff_coeffs

!> Reads the parameter "USE_CVMIX_DDIFF" and returns state.
!! This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function CVMix_ddiff_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  call get_param(param_file, mdl, "USE_CVMIX_DDIFF", CVMix_ddiff_is_used, &
                 default=.false., do_not_log=.true.)

end function CVMix_ddiff_is_used

!> Clear pointers and deallocate memory
! NOTE: Placeholder destructor
subroutine CVMix_ddiff_end(CS)
  type(CVMix_ddiff_cs), pointer :: CS !< Control structure for this module that
                                      !! will be deallocated in this subroutine
end subroutine CVMix_ddiff_end

end module MOM_CVMix_ddiff
