!> Provides the ocean grid type
module MOM_grid

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_domains, only : MOM_domain_type, get_domain_extent, compute_block_extent
use MOM_domains, only : get_global_shape, deallocate_MOM_domain
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public MOM_grid_init, MOM_grid_end, set_derived_metrics, set_first_direction
public isPointInCell, hor_index_type, get_global_grid_size

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Ocean grid type. See mom_grid for details.
type, public :: ocean_grid_type
  type(MOM_domain_type), pointer :: Domain => NULL() !< Ocean model domain
  type(MOM_domain_type), pointer :: Domain_aux => NULL() !< A non-symmetric auxiliary domain type.
  type(hor_index_type) :: HI !< Horizontal index ranges
  type(hor_index_type) :: HId2 !< Horizontal index ranges for level-2-downsampling

  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: isd_global !< The value of isd in the global index space (decomposition invariant).
  integer :: jsd_global !< The value of isd in the global index space (decomposition invariant).
  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  integer :: ke         !< The number of layers in the vertical.
  logical :: symmetric  !< True if symmetric memory is used.
  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                  !! allowed.  The default is .false. (for now).
  integer :: first_direction !< An integer that indicates which direction is
                             !! to be updated first in directionally split
                             !! parts of the calculation.  This can be altered
                             !! during the course of the run via calls to
                             !! set_first_direction.

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    mask2dT, &   !< 0 for land points and 1 for ocean points on the h-grid [nondim].
    geoLatT, &   !< The geographic latitude at tracer (h) points [degrees_N] or [km] or [m]
    geoLonT, &   !< The geographic longitude at tracer (h) points [degrees_E] or [km] or [m]
    dxT, &       !< dxT is delta x at h points [L ~> m].
    IdxT, &      !< 1/dxT [L-1 ~> m-1].
    dyT, &       !< dyT is delta y at h points [L ~> m].
    IdyT, &      !< IdyT is 1/dyT [L-1 ~> m-1].
    areaT, &     !< The area of an h-cell [L2 ~> m2].
    IareaT, &    !< 1/areaT [L-2 ~> m-2].
    sin_rot, &   !< The sine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].
    cos_rot      !< The cosine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    mask2dCu, &  !< 0 for boundary points and 1 for ocean points on the u grid [nondim].
    OBCmaskCu, & !< 0 for boundary or OBC points and 1 for ocean points on the u grid [nondim].
    geoLatCu, &  !< The geographic latitude at u points [degrees_N] or [km] or [m]
    geoLonCu, &  !< The geographic longitude at u points [degrees_E] or [km] or [m].
    dxCu, &      !< dxCu is delta x at u points [L ~> m].
    IdxCu, &     !< 1/dxCu [L-1 ~> m-1].
    dyCu, &      !< dyCu is delta y at u points [L ~> m].
    IdyCu, &     !< 1/dyCu [L-1 ~> m-1].
    dy_Cu, &     !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
    IareaCu, &   !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    areaCu       !< The areas of the u-grid cells [L2 ~> m2].

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    mask2dCv, &  !< 0 for boundary points and 1 for ocean points on the v grid [nondim].
    OBCmaskCv, & !< 0 for boundary or OBC points and 1 for ocean points on the v grid [nondim].
    geoLatCv, &  !< The geographic latitude at v points [degrees_N] or [km] or [m]
    geoLonCv, &  !< The geographic longitude at v points [degrees_E] or [km] or [m].
    dxCv, &      !< dxCv is delta x at v points [L ~> m].
    IdxCv, &     !< 1/dxCv [L-1 ~> m-1].
    dyCv, &      !< dyCv is delta y at v points [L ~> m].
    IdyCv, &     !< 1/dyCv [L-1 ~> m-1].
    dx_Cv, &     !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
    IareaCv, &   !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
    areaCv       !< The areas of the v-grid cells [L2 ~> m2].

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    porous_DminU, & !< minimum topographic height (deepest) of U-face [Z ~> m]
    porous_DmaxU, & !< maximum topographic height (shallowest) of U-face [Z ~> m]
    porous_DavgU    !< average topographic height of U-face [Z ~> m]

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    porous_DminV, & !< minimum topographic height (deepest) of V-face [Z ~> m]
    porous_DmaxV, & !< maximum topographic height (shallowest) of V-face [Z ~> m]
    porous_DavgV    !< average topographic height of V-face [Z ~> m]

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    mask2dBu, &  !< 0 for boundary points and 1 for ocean points on the q grid [nondim].
    geoLatBu, &  !< The geographic latitude at q points [degrees_N] or [km] or [m]
    geoLonBu, &  !< The geographic longitude at q points [degrees_E] or [km] or [m].
    dxBu, &      !< dxBu is delta x at q points [L ~> m].
    IdxBu, &     !< 1/dxBu [L-1 ~> m-1].
    dyBu, &      !< dyBu is delta y at q points [L ~> m].
    IdyBu, &     !< 1/dyBu [L-1 ~> m-1].
    areaBu, &    !< areaBu is the area of a q-cell [L2 ~> m2]
    IareaBu      !< IareaBu = 1/areaBu [L-2 ~> m-2].

  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & !< The latitude of T points for the purpose of labeling the output axes,
                          !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLatT.
    gridLatB => NULL()    !< The latitude of B points for the purpose of labeling the output axes,
                          !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLatBu.
  real, pointer, dimension(:) :: &
    gridLonT => NULL(), & !< The longitude of T points for the purpose of labeling the output axes,
                          !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLonT.
    gridLonB => NULL()    !< The longitude of B points for the purpose of labeling the output axes,
                          !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
                          !! On many grids this is the same as geoLonBu.
  character(len=40) :: &
    ! Except on a Cartesian grid, these are usually some variant of "degrees".
    x_axis_units, &     !< The units that are used in labeling the x coordinate axes.
    y_axis_units, &     !< The units that are used in labeling the y coordinate axes.
    ! These are internally generated names, including "m", "km", "deg_E" and "deg_N".
    x_ax_unit_short, &  !< A short description of the x-axis units for documenting parameter units
    y_ax_unit_short     !< A short description of the y-axis units for documenting parameter units

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    bathyT           !< Ocean bottom depth at tracer points, in depth units [Z ~> m].
  real    :: Z_ref   !< A reference value for all geometric height fields, such as bathyT [Z ~> m].

  logical :: bathymetry_at_vel  !< If true, there are separate values for the
                  !! basin depths at velocity points.  Otherwise the effects of
                  !! of topography are entirely determined from thickness points.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Dblock_u, &   !< Topographic depths at u-points at which the flow is blocked [Z ~> m].
    Dopen_u       !< Topographic depths at u-points at which the flow is open at width dy_Cu [Z ~> m].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Dblock_v, &   !< Topographic depths at v-points at which the flow is blocked [Z ~> m].
    Dopen_v       !< Topographic depths at v-points at which the flow is open at width dx_Cv [Z ~> m].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    CoriolisBu, & !< The Coriolis parameter at corner points [T-1 ~> s-1].
    Coriolis2Bu   !< The square of the Coriolis parameter at corner points [T-2 ~> s-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    df_dx, &      !< Derivative d/dx f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].
    df_dy         !< Derivative d/dy f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].

  ! These variables are global sums that are useful for 1-d diagnostics.
  real :: areaT_global  !< Global sum of h-cell area [L2 ~> m2]
  real :: IareaT_global !< Global sum of inverse h-cell area (1/areaT_global) [L-2 ~> m-2].

  type(unit_scale_type), pointer :: US => NULL() !< A dimensional unit scaling type


  ! These variables are for block structures.
  integer :: nblocks  !< The number of sub-PE blocks on this PE
  type(hor_index_type), pointer :: Block(:) => NULL() !< Index ranges for each block

  ! These parameters are run-time parameters that are used during some
  ! initialization routines (but not all)
  real :: grid_unit_to_L !< A factor that converts a the geoLat and geoLon variables and related
                        !! variables like len_lat and len_lon into rescaled horizontal distance
                        !! units on a Cartesian grid, in [L km ~> 1000] or [L m-1 ~> 1] or
                        !! is 0 for a non-Cartesian grid.
  real :: south_lat     !< The latitude (or y-coordinate) of the first v-line [degrees_N] or [km] or [m]
  real :: west_lon      !< The longitude (or x-coordinate) of the first u-line [degrees_E] or [km] or [m]
  real :: len_lat       !< The latitudinal (or y-coord) extent of physical domain [degrees_N] or [km] or [m]
  real :: len_lon       !< The longitudinal (or x-coord) extent of physical domain [degrees_E] or [km] or [m]
  real :: Rad_Earth_L   !< The radius of the planet in rescaled units [L ~> m]
  real :: max_depth     !< The maximum depth of the ocean in depth units [Z ~> m]
end type ocean_grid_type

contains

!> MOM_grid_init initializes the ocean grid array sizes and grid memory.
subroutine MOM_grid_init(G, param_file, US, HI, global_indexing, bathymetry_at_vel)
  type(ocean_grid_type), intent(inout) :: G          !< The horizontal grid type
  type(param_file_type), intent(in)    :: param_file !< Parameter file handle
  type(unit_scale_type), optional, pointer :: US !< A dimensional unit scaling type
  type(hor_index_type), &
                  optional, intent(in) :: HI !< A hor_index_type for array extents
  logical,        optional, intent(in) :: global_indexing !< If true use global index
                             !! values instead of having the data domain on each
                             !! processor start at 1.
  logical,        optional, intent(in) :: bathymetry_at_vel !< If true, there are
                             !! separate values for the ocean bottom depths at
                             !! velocity points.  Otherwise the effects of topography
                             !! are entirely determined from thickness points.

  ! Local variables
  real :: mean_SeaLev_scale ! A scaling factor for the reference height variable [1] or [Z m-1 ~> 1]
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  integer :: ied_max, jed_max
  integer :: niblock, njblock, nihalo, njhalo, nblocks, n, i, j
  logical :: local_indexing  ! If false use global index values instead of having
                             ! the data domain on each processor start at 1.
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  integer, allocatable, dimension(:) :: ibegin, iend, jbegin, jend
  character(len=40)  :: mod_nm  = "MOM_grid" ! This module's name.

  mean_SeaLev_scale = 1.0 ;  if (associated(G%US)) mean_SeaLev_scale = G%US%m_to_Z

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mod_nm, "REFERENCE_HEIGHT", G%Z_ref, &
                 units="m", default=0.0, scale=mean_SeaLev_scale, do_not_log=.true.)
  call log_version(param_file, mod_nm, version, &
                   "Parameters providing information about the lateral grid.", &
                   log_to_all=.true., layout=.true., all_default=(G%Z_ref==0.0))

  call get_param(param_file, mod_nm, "NIBLOCK", niblock, "The number of blocks "// &
                 "in the x-direction on each processor (for openmp).", default=1, &
                 layoutParam=.true.)
  call get_param(param_file, mod_nm, "NJBLOCK", njblock, "The number of blocks "// &
                 "in the y-direction on each processor (for openmp).", default=1, &
                 layoutParam=.true.)
  if (present(US)) then ; if (associated(US)) G%US => US ; endif

  call get_param(param_file, mod_nm, "REFERENCE_HEIGHT", G%Z_ref, &
                 "A reference value for geometric height fields, such as bathyT.", &
                 units="m", default=0.0, scale=mean_SeaLev_scale)

  if (present(HI)) then
    G%HI = HI

    G%isc = HI%isc ; G%iec = HI%iec ; G%jsc = HI%jsc ; G%jec = HI%jec
    G%isd = HI%isd ; G%ied = HI%ied ; G%jsd = HI%jsd ; G%jed = HI%jed
    G%isg = HI%isg ; G%ieg = HI%ieg ; G%jsg = HI%jsg ; G%jeg = HI%jeg

    G%IscB = HI%IscB ; G%IecB = HI%IecB ; G%JscB = HI%JscB ; G%JecB = HI%JecB
    G%IsdB = HI%IsdB ; G%IedB = HI%IedB ; G%JsdB = HI%JsdB ; G%JedB = HI%JedB
    G%IsgB = HI%IsgB ; G%IegB = HI%IegB ; G%JsgB = HI%JsgB ; G%JegB = HI%JegB

    G%idg_offset = HI%idg_offset ; G%jdg_offset = HI%jdg_offset
    G%isd_global = G%isd + HI%idg_offset ; G%jsd_global = G%jsd + HI%jdg_offset
    G%symmetric = HI%symmetric
  else
    local_indexing = .true.
    if (present(global_indexing)) local_indexing = .not.global_indexing
    call hor_index_init(G%Domain, G%HI, param_file, &
                        local_indexing=local_indexing)

    ! get_domain_extent ensures that domains start at 1 for compatibility between
    ! static and dynamically allocated arrays, unless global_indexing is true.
    call get_domain_extent(G%Domain, G%isc, G%iec, G%jsc, G%jec, &
                           G%isd, G%ied, G%jsd, G%jed, &
                           G%isg, G%ieg, G%jsg, G%jeg, &
                           G%idg_offset, G%jdg_offset, G%symmetric, &
                           local_indexing=local_indexing)
    G%isd_global = G%isd+G%idg_offset ; G%jsd_global = G%jsd+G%jdg_offset
  endif

  G%nonblocking_updates = G%Domain%nonblocking_updates

  ! Set array sizes for fields that are discretized at tracer cell boundaries.
  G%IscB = G%isc ; G%JscB = G%jsc
  G%IsdB = G%isd ; G%JsdB = G%jsd
  G%IsgB = G%isg ; G%JsgB = G%jsg
  if (G%symmetric) then
    G%IscB = G%isc-1 ; G%JscB = G%jsc-1
    G%IsdB = G%isd-1 ; G%JsdB = G%jsd-1
    G%IsgB = G%isg-1 ; G%JsgB = G%jsg-1
  endif
  G%IecB = G%iec ; G%JecB = G%jec
  G%IedB = G%ied ; G%JedB = G%jed
  G%IegB = G%ieg ; G%JegB = G%jeg

  call MOM_mesg("  MOM_grid.F90, MOM_grid_init: allocating metrics", 5)

  call allocate_metrics(G)

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  G%bathymetry_at_vel = .false.
  if (present(bathymetry_at_vel)) G%bathymetry_at_vel = bathymetry_at_vel
  if (G%bathymetry_at_vel) then
    ALLOC_(G%Dblock_u(IsdB:IedB, jsd:jed)) ; G%Dblock_u(:,:) = -G%Z_ref
    ALLOC_(G%Dopen_u(IsdB:IedB, jsd:jed))  ; G%Dopen_u(:,:) = -G%Z_ref
    ALLOC_(G%Dblock_v(isd:ied, JsdB:JedB)) ; G%Dblock_v(:,:) = -G%Z_ref
    ALLOC_(G%Dopen_v(isd:ied, JsdB:JedB))  ; G%Dopen_v(:,:) = -G%Z_ref
  endif

! setup block indices.
  nihalo = G%Domain%nihalo
  njhalo = G%Domain%njhalo
  nblocks = niblock * njblock
  if (nblocks < 1) call MOM_error(FATAL, "MOM_grid_init: " // &
       "nblocks(=NI_BLOCK*NJ_BLOCK) must be no less than 1")

  allocate(ibegin(niblock), iend(niblock), jbegin(njblock), jend(njblock))
  call compute_block_extent(G%HI%isc,G%HI%iec,niblock,ibegin,iend)
  call compute_block_extent(G%HI%jsc,G%HI%jec,njblock,jbegin,jend)
  !-- make sure the last block is the largest.
  do i = 1, niblock-1
    if (iend(i)-ibegin(i) > iend(niblock)-ibegin(niblock) ) call MOM_error(FATAL, &
       "MOM_grid_init: the last block size in x-direction is not the largest")
  enddo
  do j = 1, njblock-1
    if (jend(j)-jbegin(j) > jend(njblock)-jbegin(njblock) ) call MOM_error(FATAL, &
       "MOM_grid_init: the last block size in y-direction is not the largest")
  enddo

  G%nblocks = nblocks
  allocate(G%Block(nblocks))
  ied_max = 1 ; jed_max = 1
  do n = 1,nblocks
    ! Copy all information from the array index type describing the local grid.
    G%Block(n) = G%HI

    i = mod((n-1), niblock) + 1
    j = (n-1)/niblock + 1
    !--- isd and jsd are always 1 for each block to permit array reuse.
    G%Block(n)%isd = 1 ; G%Block(n)%jsd = 1
    G%Block(n)%isc = G%Block(n)%isd+nihalo
    G%Block(n)%jsc = G%Block(n)%jsd+njhalo
    G%Block(n)%iec = G%Block(n)%isc + iend(i) - ibegin(i)
    G%Block(n)%jec = G%Block(n)%jsc + jend(j) - jbegin(j)
    G%Block(n)%ied = G%Block(n)%iec + nihalo
    G%Block(n)%jed = G%Block(n)%jec + njhalo
    G%Block(n)%IscB = G%Block(n)%isc; G%Block(n)%IecB = G%Block(n)%iec
    G%Block(n)%JscB = G%Block(n)%jsc; G%Block(n)%JecB = G%Block(n)%jec
    !   For symmetric memory domains, the first block will have the extra point
    ! at the lower boundary of its computational domain.
    if (G%symmetric) then
      if (i==1) G%Block(n)%IscB = G%Block(n)%IscB-1
      if (j==1) G%Block(n)%JscB = G%Block(n)%JscB-1
    endif
    G%Block(n)%IsdB = G%Block(n)%isd; G%Block(n)%IedB = G%Block(n)%ied
    G%Block(n)%JsdB = G%Block(n)%jsd; G%Block(n)%JedB = G%Block(n)%jed
    !--- For symmetric memory domain, every block will have an extra point
    !--- at the lower boundary of its data domain.
    if (G%symmetric) then
      G%Block(n)%IsdB = G%Block(n)%IsdB-1
      G%Block(n)%JsdB = G%Block(n)%JsdB-1
    endif
    G%Block(n)%idg_offset = (ibegin(i) - G%Block(n)%isc) + G%HI%idg_offset
    G%Block(n)%jdg_offset = (jbegin(j) - G%Block(n)%jsc) + G%HI%jdg_offset
    ! Find the largest values of ied and jed so that all blocks will have the
    ! same size in memory.
    ied_max = max(ied_max, G%Block(n)%ied)
    jed_max = max(jed_max, G%Block(n)%jed)
  enddo

  ! Reset all of the data domain sizes to match the largest for array reuse,
  ! recalling that all block have isd=jed=1 for array reuse.
  do n = 1,nblocks
    G%Block(n)%ied = ied_max ; G%Block(n)%IedB = ied_max
    G%Block(n)%jed = jed_max ; G%Block(n)%JedB = jed_max
  enddo

  !-- do some bounds error checking
  if ( G%block(nblocks)%ied+G%block(nblocks)%idg_offset > G%HI%ied + G%HI%idg_offset ) &
        call MOM_error(FATAL, "MOM_grid_init: G%ied_bk > G%ied")
  if ( G%block(nblocks)%jed+G%block(nblocks)%jdg_offset > G%HI%jed + G%HI%jdg_offset ) &
        call MOM_error(FATAL, "MOM_grid_init: G%jed_bk > G%jed")

  call get_domain_extent(G%Domain, G%HId2%isc, G%HId2%iec, G%HId2%jsc, G%HId2%jec, &
                         G%HId2%isd, G%HId2%ied, G%HId2%jsd, G%HId2%jed, &
                         G%HId2%isg, G%HId2%ieg, G%HId2%jsg, G%HId2%jeg, coarsen=2)

  ! Set array sizes for fields that are discretized at tracer cell boundaries.
  G%HId2%IscB = G%HId2%isc ; G%HId2%JscB = G%HId2%jsc
  G%HId2%IsdB = G%HId2%isd ; G%HId2%JsdB = G%HId2%jsd
  G%HId2%IsgB = G%HId2%isg ; G%HId2%JsgB = G%HId2%jsg
  if (G%symmetric) then
    G%HId2%IscB = G%HId2%isc-1 ; G%HId2%JscB = G%HId2%jsc-1
    G%HId2%IsdB = G%HId2%isd-1 ; G%HId2%JsdB = G%HId2%jsd-1
    G%HId2%IsgB = G%HId2%isg-1 ; G%HId2%JsgB = G%HId2%jsg-1
  endif
  G%HId2%IecB = G%HId2%iec ; G%HId2%JecB = G%HId2%jec
  G%HId2%IedB = G%HId2%ied ; G%HId2%JedB = G%HId2%jed
  G%HId2%IegB = G%HId2%ieg ; G%HId2%JegB = G%HId2%jeg

end subroutine MOM_grid_init

!> set_derived_metrics calculates metric terms that are derived from other metrics.
subroutine set_derived_metrics(G, US)
  type(ocean_grid_type), intent(inout) :: G  !< The horizontal grid structure
  type(unit_scale_type), intent(in)    :: US !< A dimensional unit scaling type
!    Various inverse grid spacings and derived areas are calculated within this
!  subroutine.
  integer :: i, j, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do j=jsd,jed ; do i=isd,ied
    if (G%dxT(i,j) < 0.0) G%dxT(i,j) = 0.0
    if (G%dyT(i,j) < 0.0) G%dyT(i,j) = 0.0
    G%IdxT(i,j) = Adcroft_reciprocal(G%dxT(i,j))
    G%IdyT(i,j) = Adcroft_reciprocal(G%dyT(i,j))
    G%IareaT(i,j) = Adcroft_reciprocal(G%areaT(i,j))
  enddo ; enddo

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dxCu(I,j) < 0.0) G%dxCu(I,j) = 0.0
    if (G%dyCu(I,j) < 0.0) G%dyCu(I,j) = 0.0
    G%IdxCu(I,j) = Adcroft_reciprocal(G%dxCu(I,j))
    G%IdyCu(I,j) = Adcroft_reciprocal(G%dyCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dxCv(i,J) < 0.0) G%dxCv(i,J) = 0.0
    if (G%dyCv(i,J) < 0.0) G%dyCv(i,J) = 0.0
    G%IdxCv(i,J) = Adcroft_reciprocal(G%dxCv(i,J))
    G%IdyCv(i,J) = Adcroft_reciprocal(G%dyCv(i,J))
  enddo ; enddo

  do J=JsdB,JedB ; do I=IsdB,IedB
    if (G%dxBu(I,J) < 0.0) G%dxBu(I,J) = 0.0
    if (G%dyBu(I,J) < 0.0) G%dyBu(I,J) = 0.0

    G%IdxBu(I,J) = Adcroft_reciprocal(G%dxBu(I,J))
    G%IdyBu(I,J) = Adcroft_reciprocal(G%dyBu(I,J))
    ! areaBu has usually been set to a positive area elsewhere.
    if (G%areaBu(I,J) <= 0.0) G%areaBu(I,J) = G%dxBu(I,J) * G%dyBu(I,J)
    G%IareaBu(I,J) =  Adcroft_reciprocal(G%areaBu(I,J))
  enddo ; enddo
end subroutine set_derived_metrics

!> Adcroft_reciprocal(x) = 1/x for |x|>0 or 0 for x=0.
function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted [A].
  real :: I_val            !< The Adcroft reciprocal of val [A-1].

  I_val = 0.0 ; if (val /= 0.0) I_val = 1.0/val
end function Adcroft_reciprocal

!> Returns true if the coordinates (x,y) are within the h-cell (i,j)
logical function isPointInCell(G, i, j, x, y)
  type(ocean_grid_type), intent(in) :: G !< Grid type
  integer,               intent(in) :: i !< i index of cell to test
  integer,               intent(in) :: j !< j index of cell to test
  real,                  intent(in) :: x !< x coordinate of point [degrees_E]
  real,                  intent(in) :: y !< y coordinate of point [degrees_N]
  ! Local variables
  real :: xNE, xNW, xSE, xSW ! Longitudes of cell corners [degrees_E]
  real :: yNE, yNW, ySE, ySW ! Latitudes of cell corners [degrees_N]
  real :: l0, l1, l2, l3 ! Crossed products of differences in position [degrees_E degrees_N]
  real :: p0, p1, p2, p3 ! Trinary unitary values reflecting the signs of the crossed products [nondim]
  isPointInCell = .false.
  xNE = G%geoLonBu(i  ,j  ) ; yNE = G%geoLatBu(i  ,j  )
  xNW = G%geoLonBu(i-1,j  ) ; yNW = G%geoLatBu(i-1,j  )
  xSE = G%geoLonBu(i  ,j-1) ; ySE = G%geoLatBu(i  ,j-1)
  xSW = G%geoLonBu(i-1,j-1) ; ySW = G%geoLatBu(i-1,j-1)
  ! This is a crude calculation that assumes a geographic coordinate system
  if (x<min(xNE,xNW,xSE,xSW) .or. x>max(xNE,xNW,xSE,xSW) .or. &
      y<min(yNE,yNW,ySE,ySW) .or. y>max(yNE,yNW,ySE,ySW) ) then
    return ! Avoid the more complicated calculation
  endif
  l0 = (x-xSW)*(ySE-ySW) - (y-ySW)*(xSE-xSW)
  l1 = (x-xSE)*(yNE-ySE) - (y-ySE)*(xNE-xSE)
  l2 = (x-xNE)*(yNW-yNE) - (y-yNE)*(xNW-xNE)
  l3 = (x-xNW)*(ySW-yNW) - (y-yNW)*(xSW-xNW)

  p0 = sign(1., l0) ; if (l0 == 0.) p0=0.
  p1 = sign(1., l1) ; if (l1 == 0.) p1=0.
  p2 = sign(1., l2) ; if (l2 == 0.) p2=0.
  p3 = sign(1., l3) ; if (l3 == 0.) p3=0.

  if ( (abs(p0)+abs(p2)) + (abs(p1)+abs(p3)) == abs((p0+p2) + (p1+p3)) ) then
    isPointInCell=.true.
  endif
end function isPointInCell

!> Store an integer indicating which direction to work on first.
subroutine set_first_direction(G, y_first)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  integer,               intent(in) :: y_first !< The first direction to store

  G%first_direction = y_first
end subroutine set_first_direction

!> Return global shape of horizontal grid
subroutine get_global_grid_size(G, niglobal, njglobal)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type
  integer,               intent(out)   :: niglobal !< i-index global size of grid
  integer,               intent(out)   :: njglobal !< j-index global size of grid

  call get_global_shape(G%domain, niglobal, njglobal)

end subroutine get_global_grid_size

!> Allocate memory used by the ocean_grid_type and related structures.
subroutine allocate_metrics(G)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg

  ! This subroutine allocates the lateral elements of the ocean_grid_type that
  ! are always used and zeros them out.

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  ALLOC_(G%dxT(isd:ied,jsd:jed))       ; G%dxT(:,:) = 0.0
  ALLOC_(G%dxCu(IsdB:IedB,jsd:jed))    ; G%dxCu(:,:) = 0.0
  ALLOC_(G%dxCv(isd:ied,JsdB:JedB))    ; G%dxCv(:,:) = 0.0
  ALLOC_(G%dxBu(IsdB:IedB,JsdB:JedB))  ; G%dxBu(:,:) = 0.0
  ALLOC_(G%IdxT(isd:ied,jsd:jed))      ; G%IdxT(:,:) = 0.0
  ALLOC_(G%IdxCu(IsdB:IedB,jsd:jed))   ; G%IdxCu(:,:) = 0.0
  ALLOC_(G%IdxCv(isd:ied,JsdB:JedB))   ; G%IdxCv(:,:) = 0.0
  ALLOC_(G%IdxBu(IsdB:IedB,JsdB:JedB)) ; G%IdxBu(:,:) = 0.0

  ALLOC_(G%dyT(isd:ied,jsd:jed))       ; G%dyT(:,:) = 0.0
  ALLOC_(G%dyCu(IsdB:IedB,jsd:jed))    ; G%dyCu(:,:) = 0.0
  ALLOC_(G%dyCv(isd:ied,JsdB:JedB))    ; G%dyCv(:,:) = 0.0
  ALLOC_(G%dyBu(IsdB:IedB,JsdB:JedB))  ; G%dyBu(:,:) = 0.0
  ALLOC_(G%IdyT(isd:ied,jsd:jed))      ; G%IdyT(:,:) = 0.0
  ALLOC_(G%IdyCu(IsdB:IedB,jsd:jed))   ; G%IdyCu(:,:) = 0.0
  ALLOC_(G%IdyCv(isd:ied,JsdB:JedB))   ; G%IdyCv(:,:) = 0.0
  ALLOC_(G%IdyBu(IsdB:IedB,JsdB:JedB)) ; G%IdyBu(:,:) = 0.0

  ALLOC_(G%areaT(isd:ied,jsd:jed))       ; G%areaT(:,:) = 0.0
  ALLOC_(G%IareaT(isd:ied,jsd:jed))      ; G%IareaT(:,:) = 0.0
  ALLOC_(G%areaBu(IsdB:IedB,JsdB:JedB))  ; G%areaBu(:,:) = 0.0
  ALLOC_(G%IareaBu(IsdB:IedB,JsdB:JedB)) ; G%IareaBu(:,:) = 0.0

  ALLOC_(G%mask2dT(isd:ied,jsd:jed))      ; G%mask2dT(:,:) = 0.0
  ALLOC_(G%mask2dCu(IsdB:IedB,jsd:jed))   ; G%mask2dCu(:,:) = 0.0
  ALLOC_(G%OBCmaskCu(IsdB:IedB,jsd:jed))  ; G%OBCmaskCu(:,:) = 0.0
  ALLOC_(G%mask2dCv(isd:ied,JsdB:JedB))   ; G%mask2dCv(:,:) = 0.0
  ALLOC_(G%OBCmaskCv(isd:ied,JsdB:JedB))  ; G%OBCmaskCv(:,:) = 0.0
  ALLOC_(G%mask2dBu(IsdB:IedB,JsdB:JedB)) ; G%mask2dBu(:,:) = 0.0
  ALLOC_(G%geoLatT(isd:ied,jsd:jed))      ; G%geoLatT(:,:) = 0.0
  ALLOC_(G%geoLatCu(IsdB:IedB,jsd:jed))   ; G%geoLatCu(:,:) = 0.0
  ALLOC_(G%geoLatCv(isd:ied,JsdB:JedB))   ; G%geoLatCv(:,:) = 0.0
  ALLOC_(G%geoLatBu(IsdB:IedB,JsdB:JedB)) ; G%geoLatBu(:,:) = 0.0
  ALLOC_(G%geoLonT(isd:ied,jsd:jed))      ; G%geoLonT(:,:) = 0.0
  ALLOC_(G%geoLonCu(IsdB:IedB,jsd:jed))   ; G%geoLonCu(:,:) = 0.0
  ALLOC_(G%geoLonCv(isd:ied,JsdB:JedB))   ; G%geoLonCv(:,:) = 0.0
  ALLOC_(G%geoLonBu(IsdB:IedB,JsdB:JedB)) ; G%geoLonBu(:,:) = 0.0

  ALLOC_(G%dx_Cv(isd:ied,JsdB:JedB))     ; G%dx_Cv(:,:) = 0.0
  ALLOC_(G%dy_Cu(IsdB:IedB,jsd:jed))     ; G%dy_Cu(:,:) = 0.0

  ALLOC_(G%porous_DminU(IsdB:IedB,jsd:jed)); G%porous_DminU(:,:) = 0.0
  ALLOC_(G%porous_DmaxU(IsdB:IedB,jsd:jed)); G%porous_DmaxU(:,:) = 0.0
  ALLOC_(G%porous_DavgU(IsdB:IedB,jsd:jed)); G%porous_DavgU(:,:) = 0.0

  ALLOC_(G%porous_DminV(isd:ied,JsdB:JedB)); G%porous_DminV(:,:) = 0.0
  ALLOC_(G%porous_DmaxV(isd:ied,JsdB:JedB)); G%porous_DmaxV(:,:) = 0.0
  ALLOC_(G%porous_DavgV(isd:ied,JsdB:JedB)); G%porous_DavgV(:,:) = 0.0

  ALLOC_(G%areaCu(IsdB:IedB,jsd:jed))  ; G%areaCu(:,:) = 0.0
  ALLOC_(G%areaCv(isd:ied,JsdB:JedB))  ; G%areaCv(:,:) = 0.0
  ALLOC_(G%IareaCu(IsdB:IedB,jsd:jed)) ; G%IareaCu(:,:) = 0.0
  ALLOC_(G%IareaCv(isd:ied,JsdB:JedB)) ; G%IareaCv(:,:) = 0.0

  ALLOC_(G%bathyT(isd:ied, jsd:jed)) ; G%bathyT(:,:) = -G%Z_ref
  ALLOC_(G%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; G%CoriolisBu(:,:) = 0.0
  ALLOC_(G%Coriolis2Bu(IsdB:IedB, JsdB:JedB)) ; G%Coriolis2Bu(:,:) = 0.0
  ALLOC_(G%dF_dx(isd:ied, jsd:jed)) ; G%dF_dx(:,:) = 0.0
  ALLOC_(G%dF_dy(isd:ied, jsd:jed)) ; G%dF_dy(:,:) = 0.0

  ALLOC_(G%sin_rot(isd:ied,jsd:jed)) ; G%sin_rot(:,:) = 0.0
  ALLOC_(G%cos_rot(isd:ied,jsd:jed)) ; G%cos_rot(:,:) = 1.0

  allocate(G%gridLonT(isg:ieg), source=0.0)
  allocate(G%gridLonB(G%IsgB:G%IegB), source=0.0)
  allocate(G%gridLatT(jsg:jeg), source=0.0)
  allocate(G%gridLatB(G%JsgB:G%JegB), source=0.0)

end subroutine allocate_metrics

!> Release memory used by the ocean_grid_type and related structures.
subroutine MOM_grid_end(G)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type

  deallocate(G%Block)

  if (G%bathymetry_at_vel) then
    DEALLOC_(G%Dblock_u) ; DEALLOC_(G%Dopen_u)
    DEALLOC_(G%Dblock_v) ; DEALLOC_(G%Dopen_v)
  endif

  DEALLOC_(G%dxT)  ; DEALLOC_(G%dxCu)  ; DEALLOC_(G%dxCv)  ; DEALLOC_(G%dxBu)
  DEALLOC_(G%IdxT) ; DEALLOC_(G%IdxCu) ; DEALLOC_(G%IdxCv) ; DEALLOC_(G%IdxBu)

  DEALLOC_(G%dyT)  ; DEALLOC_(G%dyCu)  ; DEALLOC_(G%dyCv)  ; DEALLOC_(G%dyBu)
  DEALLOC_(G%IdyT) ; DEALLOC_(G%IdyCu) ; DEALLOC_(G%IdyCv) ; DEALLOC_(G%IdyBu)

  DEALLOC_(G%areaT)  ; DEALLOC_(G%IareaT)
  DEALLOC_(G%areaBu) ; DEALLOC_(G%IareaBu)
  DEALLOC_(G%areaCu) ; DEALLOC_(G%IareaCu)
  DEALLOC_(G%areaCv)  ; DEALLOC_(G%IareaCv)

  DEALLOC_(G%mask2dT)  ; DEALLOC_(G%mask2dCu) ; DEALLOC_(G%OBCmaskCu)
  DEALLOC_(G%mask2dCv) ; DEALLOC_(G%OBCmaskCv) ; DEALLOC_(G%mask2dBu)

  DEALLOC_(G%geoLatT)  ; DEALLOC_(G%geoLatCu)
  DEALLOC_(G%geoLatCv) ; DEALLOC_(G%geoLatBu)
  DEALLOC_(G%geoLonT)  ; DEALLOC_(G%geoLonCu)
  DEALLOC_(G%geoLonCv) ; DEALLOC_(G%geoLonBu)

  DEALLOC_(G%dx_Cv) ; DEALLOC_(G%dy_Cu)

  DEALLOC_(G%bathyT)  ; DEALLOC_(G%CoriolisBu) ; DEALLOC_(G%Coriolis2Bu)
  DEALLOC_(G%dF_dx)   ; DEALLOC_(G%dF_dy)
  DEALLOC_(G%sin_rot) ; DEALLOC_(G%cos_rot)

  DEALLOC_(G%porous_DminU) ; DEALLOC_(G%porous_DmaxU) ; DEALLOC_(G%porous_DavgU)
  DEALLOC_(G%porous_DminV) ; DEALLOC_(G%porous_DmaxV) ; DEALLOC_(G%porous_DavgV)

  deallocate(G%gridLonT) ; deallocate(G%gridLatT)
  deallocate(G%gridLonB) ; deallocate(G%gridLatB)

  ! The cursory flag avoids doing any deallocation of memory in the underlying
  ! infrastructure to avoid problems due to shared pointers.
  call deallocate_MOM_domain(G%Domain, cursory=.true.)

end subroutine MOM_grid_end

!> \namespace mom_grid
!!
!! Grid metrics and their inverses are labelled according to their staggered location on a Arakawa C (or B) grid.
!! - Metrics centered on h- or T-points are labelled T, e.g. dxT is the distance across the cell in the x-direction.
!! - Metrics centered on u-points are labelled Cu (C-grid u location). e.g. dyCu is the y-distance between
!!   two corners of a T-cell.
!! - Metrics centered on v-points are labelled Cv (C-grid v location). e.g. dyCv is the y-distance between two -points.
!! - Metrics centered on q-points are labelled Bu (B-grid u,v location). e.g. areaBu is the area centered on a q-point.
!!
!! \image html Grid_metrics.png "The labelling of distances (grid metrics) at various staggered
!! location on an T-cell and around a q-point."
!!
!! Areas centered at T-, u-, v- and q- points are `areaT`, `areaCu`, `areaCv` and `areaBu` respectively.
!!
!! The reciprocal of metrics are pre-calculated and also stored in the ocean_grid_type with a I prepended to the name.
!! For example, `1./areaT` is called `IareaT`, and `1./dyCv` is `IdyCv`.
!!
!! Geographic latitude and longitude (or model coordinates if not on a sphere) are stored in
!! `geoLatT`, `geoLonT` for T-points.
!! u-, v- and q- point coordinates are follow same pattern of replacing T with Cu, Cv and Bu respectively.
!!
!! Each location also has a 2D mask indicating whether the entire column is land or ocean.
!! `mask2dT` is 1 if the column is wet or 0 if the T-cell is land.
!! `mask2dCu` is 1 if both neighboring columns are ocean, and 0 if either is land.
!! `OBCmasku` is 1 if both neighboring columns are ocean, and 0 if either is land of if this is OBC point.

end module MOM_grid
