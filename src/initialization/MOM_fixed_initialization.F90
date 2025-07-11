!> Initializes fixed aspects of the model, such as horizontal grid metrics,
!! topography and Coriolis.
module MOM_fixed_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum, qchksum, uvchksum
use MOM_domains, only : pass_var
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_io, only : slasher
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_open_boundary, only : ocean_OBC_type
use MOM_open_boundary, only : open_boundary_config, open_boundary_query
use MOM_open_boundary, only : open_boundary_impose_normal_slope
use MOM_open_boundary, only : open_boundary_impose_land_mask
use MOM_shared_initialization, only : MOM_initialize_rotation, MOM_calculate_grad_Coriolis
use MOM_shared_initialization, only : initialize_topography_from_file, apply_topography_edits_from_file
use MOM_shared_initialization, only : initialize_topography_named, limit_topography, diagnoseMaximumDepth
use MOM_shared_initialization, only : set_rotation_planetary, set_rotation_beta_plane, initialize_grid_rotation_angle
use MOM_shared_initialization, only : reset_face_lengths_named, reset_face_lengths_file, reset_face_lengths_list
use MOM_shared_initialization, only : read_face_length_list, set_velocity_depth_max, set_velocity_depth_min
use MOM_shared_initialization, only : set_subgrid_topo_at_vel_from_file
use MOM_shared_initialization, only : compute_global_grid_integrals
use MOM_unit_scaling, only : unit_scale_type

use user_initialization, only : user_initialize_topography
use DOME_initialization, only : DOME_initialize_topography
use ISOMIP_initialization, only : ISOMIP_initialize_topography
use basin_builder, only : basin_builder_topography
use benchmark_initialization, only : benchmark_initialize_topography
use Neverworld_initialization, only : Neverworld_initialize_topography
use DOME2d_initialization, only : DOME2d_initialize_topography
use Kelvin_initialization, only : Kelvin_initialize_topography
use sloshing_initialization, only : sloshing_initialize_topography
use seamount_initialization, only : seamount_initialize_topography
use dumbbell_initialization, only : dumbbell_initialize_topography
use shelfwave_initialization, only : shelfwave_initialize_topography
use Phillips_initialization, only : Phillips_initialize_topography
use dense_water_initialization, only : dense_water_initialize_topography

implicit none ; private

public MOM_initialize_fixed, MOM_initialize_rotation, MOM_initialize_topography

contains

! -----------------------------------------------------------------------------
!> MOM_initialize_fixed sets up time-invariant quantities related to MOM6's
!!   horizontal grid, bathymetry, and the Coriolis parameter.
subroutine MOM_initialize_fixed(G, US, OBC, PF)
  type(dyn_horgrid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundary structure.
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.

  ! Local variables
  character(len=200) :: inputdir   ! The directory where NetCDF input files are.
  character(len=200) :: config
  logical            :: read_porous_file
  character(len=40)  :: mdl = "MOM_fixed_initialization" ! This module's name.
  integer :: I, J
  logical :: debug
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  call callTree_enter("MOM_initialize_fixed(), MOM_fixed_initialization.F90")
  call log_version(PF, mdl, version, "")
  call get_param(PF, mdl, "DEBUG", debug, default=.false.)

  call get_param(PF, mdl, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

  ! Set up the parameters of the physical domain (i.e. the grid), G
  call set_grid_metrics(G, PF, US)

  ! Set up the bottom depth, G%bathyT either analytically or from file
  ! This also sets G%max_depth based on the input parameter MAXIMUM_DEPTH,
  ! or, if absent, is diagnosed as G%max_depth = max( G%D(:,:) )
  call MOM_initialize_topography(G%bathyT, G%max_depth, G, PF, US)

  ! To initialize masks, the bathymetry in halo regions must be filled in
  call pass_var(G%bathyT, G%Domain)

  ! Determine the position of any open boundaries
  call open_boundary_config(G, US, PF, OBC)

  ! Make bathymetry consistent with open boundaries
  call open_boundary_impose_normal_slope(OBC, G, G%bathyT)

  ! This call sets masks that prohibit flow over any point interpreted as land
  call initialize_masks(G, PF, US)

  ! Make OBC mask consistent with land mask
  call open_boundary_impose_land_mask(OBC, G, G%areaCu, G%areaCv, US)

  if (debug) then
    call hchksum(G%bathyT, 'MOM_initialize_fixed: depth ', G%HI, haloshift=1, unscale=US%Z_to_m)
    call hchksum(G%mask2dT, 'MOM_initialize_fixed: mask2dT ', G%HI)
    call uvchksum('MOM_initialize_fixed: mask2dC[uv]', G%mask2dCu, &
                  G%mask2dCv, G%HI)
    call qchksum(G%mask2dBu, 'MOM_initialize_fixed: mask2dBu ', G%HI)
  endif

  ! Modulate geometric scales according to geography.
  call get_param(PF, mdl, "CHANNEL_CONFIG", config, &
                 "A parameter that determines which set of channels are \n"//&
                 "restricted to specific  widths.  Options are:\n"//&
                 " \t none - All channels have the grid width.\n"//&
                 " \t global_1deg - Sets 16 specific channels appropriate \n"//&
                 " \t\t for a 1-degree model, as used in CM2G.\n"//&
                 " \t list - Read the channel locations and widths from a \n"//&
                 " \t\t text file, like MOM_channel_list in the MOM_SIS \n"//&
                 " \t\t test case.\n"//&
                 " \t file - Read open face widths everywhere from a \n"//&
                 " \t\t NetCDF file on the model grid.", &
                 default="none")
  select case ( trim(config) )
    case ("none")
    case ("list") ; call reset_face_lengths_list(G, PF, US)
    case ("file") ; call reset_face_lengths_file(G, PF, US)
    case ("global_1deg") ; call reset_face_lengths_named(G, PF, trim(config), US)
    case default ; call MOM_error(FATAL, "MOM_initialize_fixed: "// &
      "Unrecognized channel configuration "//trim(config))
  end select

  !   This call sets the topography at velocity points.
  if (G%bathymetry_at_vel) then
    call get_param(PF, mdl, "VELOCITY_DEPTH_CONFIG", config, &
                   "A string that determines how the topography is set at "//&
                   "velocity points. This may be 'min' or 'max'.", &
                   default="max")
    select case ( trim(config) )
      case ("max") ; call set_velocity_depth_max(G)
      case ("min") ; call set_velocity_depth_min(G)
      case default ; call MOM_error(FATAL, "MOM_initialize_fixed: "// &
        "Unrecognized velocity depth configuration "//trim(config))
    end select
  endif

  ! Read sub-grid scale topography parameters at velocity points used for porous barrier calculation
  ! TODO: The following routine call may eventually be merged as one of the CHANNEL_CONFIG options
  call get_param(PF, mdl, "SUBGRID_TOPO_AT_VEL", read_porous_file, &
                 "If true, use variables from TOPO_AT_VEL_FILE as parameters for porous barrier.", &
                 default=.False.)
  if (read_porous_file) &
    call set_subgrid_topo_at_vel_from_file(G, PF, US)

!    Calculate the value of the Coriolis parameter at the latitude   !
!  of the q grid points [T-1 ~> s-1].
  call MOM_initialize_rotation(G%CoriolisBu, G, PF, US=US)
!   Calculate the components of grad f (beta)
  call MOM_calculate_grad_Coriolis(G%dF_dx, G%dF_dy, G, US=US)
!   Calculate the square of the Coriolis parameter
  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    G%Coriolis2Bu(I,J) = G%CoriolisBu(I,J)**2
  enddo ; enddo

  if (debug) then
    call qchksum(G%CoriolisBu, "MOM_initialize_fixed: f ", G%HI, unscale=US%s_to_T)
    call qchksum(G%Coriolis2Bu, "MOM_initialize_fixed: f2 ", G%HI, unscale=US%s_to_T**2)
    call hchksum(G%dF_dx, "MOM_initialize_fixed: dF_dx ", G%HI, unscale=US%m_to_L*US%s_to_T)
    call hchksum(G%dF_dy, "MOM_initialize_fixed: dF_dy ", G%HI, unscale=US%m_to_L*US%s_to_T)
  endif

  call initialize_grid_rotation_angle(G, PF)

! Compute global integrals of grid values for later use in scalar diagnostics !
  call compute_global_grid_integrals(G, US=US)

  call callTree_leave('MOM_initialize_fixed()')

end subroutine MOM_initialize_fixed

!> MOM_initialize_topography makes the appropriate call to set up the bathymetry in units of [Z ~> m].
subroutine MOM_initialize_topography(D, max_depth, G, PF, US)
  type(dyn_horgrid_type),           intent(in)  :: G  !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D  !< Ocean bottom depth [Z ~> m]
  type(param_file_type),            intent(in)  :: PF !< Parameter file structure
  real,                             intent(out) :: max_depth !< Maximum depth of model [Z ~> m]
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type

  ! This subroutine makes the appropriate call to set up the bottom depth.
  ! This is a separate subroutine so that it can be made public and shared with
  ! the ice-sheet code or other components.

  ! Local variables
  character(len=40)  :: mdl = "MOM_initialize_topography" ! This subroutine's name.
  character(len=200) :: config

  call get_param(PF, mdl, "TOPO_CONFIG", config, &
                 "This specifies how bathymetry is specified: \n"//&
                 " \t file - read bathymetric information from the file \n"//&
                 " \t\t specified by (TOPO_FILE).\n"//&
                 " \t flat - flat bottom set to MAXIMUM_DEPTH. \n"//&
                 " \t bowl - an analytically specified bowl-shaped basin \n"//&
                 " \t\t ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH. \n"//&
                 " \t spoon - a similar shape to 'bowl', but with an vertical \n"//&
                 " \t\t wall at the southern face. \n"//&
                 " \t halfpipe - a zonally uniform channel with a half-sine \n"//&
                 " \t\t profile in the meridional direction. \n"//&
                 " \t bbuilder - build topography from list of functions. \n"//&
                 " \t benchmark - use the benchmark test case topography. \n"//&
                 " \t Neverworld - use the Neverworld test case topography. \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t ISOMIP - use a slope and channel configuration for the \n"//&
                 " \t\t ISOMIP test case. \n"//&
                 " \t DOME2D - use a shelf and slope configuration for the \n"//&
                 " \t\t DOME2D gravity current/overflow test case. \n"//&
                 " \t Kelvin - flat but with rotated land mask.\n"//&
                 " \t seamount - Gaussian bump for spontaneous motion test case.\n"//&
                 " \t dumbbell - Sloshing channel with reservoirs on both ends.\n"//&
                 " \t shelfwave - exponential slope for shelfwave test case.\n"//&
                 " \t Phillips - ACC-like idealized topography used in the Phillips config.\n"//&
                 " \t dense - Denmark Strait-like dense water formation and overflow.\n"//&
                 " \t USER - call a user modified routine.", &
                 fail_if_missing=.true.)
  call get_param(PF, mdl, "MAXIMUM_DEPTH", max_depth, units="m", default=-1.e9, scale=US%m_to_Z, do_not_log=.true.)
  select case ( trim(config) )
    case ("file");      call initialize_topography_from_file(D, G, PF, US)
    case ("flat");      call initialize_topography_named(D, G, PF, config, max_depth, US)
    case ("spoon");     call initialize_topography_named(D, G, PF, config, max_depth, US)
    case ("bowl");      call initialize_topography_named(D, G, PF, config, max_depth, US)
    case ("halfpipe");  call initialize_topography_named(D, G, PF, config, max_depth, US)
    case ("DOME");      call DOME_initialize_topography(D, G, PF, max_depth, US)
    case ("ISOMIP");    call ISOMIP_initialize_topography(D, G, PF, max_depth, US)
    case ("bbuilder");  call basin_builder_topography(D, G, PF, max_depth)
    case ("benchmark"); call benchmark_initialize_topography(D, G, PF, max_depth, US)
    case ("Neverworld","Neverland"); call Neverworld_initialize_topography(D, G, PF, max_depth)
    case ("DOME2D");    call DOME2d_initialize_topography(D, G, PF, max_depth)
    case ("Kelvin");    call Kelvin_initialize_topography(D, G, PF, max_depth, US)
    case ("sloshing");  call sloshing_initialize_topography(D, G, PF, max_depth)
    case ("seamount");  call seamount_initialize_topography(D, G, PF, max_depth)
    case ("dumbbell");  call dumbbell_initialize_topography(D, G, PF, max_depth)
    case ("shelfwave"); call shelfwave_initialize_topography(D, G, PF, max_depth, US)
    case ("Phillips");  call Phillips_initialize_topography(D, G, PF, max_depth, US)
    case ("dense");     call dense_water_initialize_topography(D, G, PF, max_depth)
    case ("USER");      call user_initialize_topography(D, G, PF, max_depth, US)
    case default ;      call MOM_error(FATAL,"MOM_initialize_topography: "// &
      "Unrecognized topography setup '"//trim(config)//"'")
  end select
  if (max_depth>0.) then
    call log_param(PF, mdl, "MAXIMUM_DEPTH", max_depth, &
                   "The maximum depth of the ocean.", units="m", unscale=US%Z_to_m)
  else
    max_depth = diagnoseMaximumDepth(D,G)
    call log_param(PF, mdl, "!MAXIMUM_DEPTH", max_depth, &
                   "The (diagnosed) maximum depth of the ocean.", &
                   units="m", unscale=US%Z_to_m, like_default=.true.)
  endif
  if (trim(config) /= "DOME") then
    call limit_topography(D, G, PF, max_depth, US)
  endif

end subroutine MOM_initialize_topography

end module MOM_fixed_initialization
