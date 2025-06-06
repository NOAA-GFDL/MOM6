!> Tidal contributions to geopotential
module MOM_tidal_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_harmonic_analysis, &
                       only : HA_init, HA_register, harmonic_analysis_CS
use MOM_io,            only : field_exists, file_exists, MOM_read_data
use MOM_time_manager,  only : set_date, time_type, time_type_to_real, operator(-)
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public calc_tidal_forcing, tidal_forcing_init, tidal_forcing_end
public calc_tidal_forcing_legacy
! MOM_open_boundary uses the following to set tides on the boundary.
public astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency

#include <MOM_memory.h>

integer, parameter :: MAX_CONSTITUENTS = 10 !< The maximum number of tidal
                                            !! constituents that could be used.
!> Simple type to store astronomical longitudes used to calculate tidal phases.
type, public :: astro_longitudes
  real :: s  !< Mean longitude of moon [rad]
  real :: h  !< Mean longitude of sun [rad]
  real :: p  !< Mean longitude of lunar perigee [rad]
  real :: N  !< Longitude of ascending node [rad]
end type astro_longitudes

!> The control structure for the MOM_tidal_forcing module
type, public :: tidal_forcing_CS ; private
  logical :: use_tidal_sal_file !< If true, Read the tidal self-attraction
                      !! and loading from input files, specified
                      !! by TIDAL_INPUT_FILE.
  logical :: use_tidal_sal_prev !< If true, use the SAL from the previous
                      !! iteration of the tides to facilitate convergence.
  logical :: use_eq_phase !< If true, tidal forcing is phase-shifted to match
                      !! equilibrium tide. Set to false if providing tidal phases
                      !! that have already been shifted by the
                      !! astronomical/equilibrium argument.
  real    :: sal_scalar = 0.0 !< The constant of proportionality between self-attraction and
                      !! loading (SAL) geopotential anomaly and total geopotential geopotential
                      !! anomalies. This is only used if USE_PREVIOUS_TIDES is true. [nondim].
  integer :: nc       !< The number of tidal constituents in use.
  real, dimension(MAX_CONSTITUENTS) :: &
    freq, &           !< The frequency of a tidal constituent [rad T-1 ~> rad s-1].
    phase0, &         !< The phase of a tidal constituent at time 0 [rad].
    amp, &            !< The amplitude of a tidal constituent at time 0 [Z ~> m].
    love_no           !< The Love number of a tidal constituent at time 0 [nondim].
  integer :: struct(MAX_CONSTITUENTS) !< An encoded spatial structure for each constituent
  character (len=16) :: const_name(MAX_CONSTITUENTS) !< The name of each constituent

  type(time_type) :: time_ref !< Reference time (t = 0) used to calculate tidal forcing.
  type(astro_longitudes) :: tidal_longitudes !< Astronomical longitudes used to calculate
                                   !! tidal phases at t = 0.
  real, allocatable :: &
    sin_struct(:,:,:), &    !< The sine based structures that can be associated with
                            !! the astronomical forcing [nondim].
    cos_struct(:,:,:), &    !< The cosine based structures that can be associated with
                            !! the astronomical forcing [nondim].
    cosphasesal(:,:,:), &   !< The cosine of the phase of the self-attraction and loading amphidromes [nondim].
    sinphasesal(:,:,:), &   !< The sine of the phase of the self-attraction and loading amphidromes [nondim].
    ampsal(:,:,:), &        !< The amplitude of the SAL [Z ~> m].
    cosphase_prev(:,:,:), & !< The cosine of the phase of the amphidromes in the previous tidal solutions [nondim].
    sinphase_prev(:,:,:), & !< The sine of the phase of the amphidromes in the previous tidal solutions [nondim].
    amp_prev(:,:,:), &      !< The amplitude of the previous tidal solution [Z ~> m].
    tide_fn(:), &           !< Amplitude modulation of tides by nodal cycle [nondim].
    tide_un(:)              !< Phase modulation of tides by nodal cycle [rad].
end type tidal_forcing_CS

integer :: id_clock_tides !< CPU clock for tides

contains

!> Finds astronomical longitudes s, h, p, and N,
!! the mean longitude of the moon, sun, lunar perigee, and ascending node, respectively,
!! at the specified reference time time_ref.
!! These formulas were obtained from
!! Kowalik and Luick, "Modern Theory and Practice of Tide Analysis and Tidal Power", 2019
!! (their Equation I.71), which are based on Schureman, 1958.
!! For simplicity, the time associated with time_ref should
!! be at midnight. These formulas also only make sense if
!! the calendar is Gregorian.
subroutine astro_longitudes_init(time_ref, longitudes)
  type(time_type), intent(in) :: time_ref            !> Time to calculate longitudes for.
  type(astro_longitudes), intent(out) :: longitudes  !> Lunar and solar longitudes at time_ref.

  ! Local variables
  real :: D                                          !> Time since the reference date [days]
  real :: T                                          !> Time in Julian centuries [centuries]
  real, parameter :: PI = 4.0 * atan(1.0)            !> 3.14159... [nondim]

  ! Find date at time_ref in days since midnight at the start of 1900-01-01
  D = time_type_to_real(time_ref - set_date(1900, 1, 1, 0, 0, 0)) / (24.0 * 3600.0)
  ! Time since 1900-01-01 in Julian centuries
  ! Kowalik and Luick use 36526, but Schureman uses 36525 which I think is correct.
  T = D / 36525.0
  ! Calculate longitudes, including converting to radians on [0, 2pi)
  ! s: Mean longitude of moon
  longitudes%s = mod((277.0248 + 481267.8906 * T) + 0.0011 * (T**2), 360.0) * PI / 180.0
  ! h: Mean longitude of sun
  longitudes%h = mod((280.1895 + 36000.7689 * T) + 3.0310e-4 * (T**2), 360.0) * PI / 180.0
  ! p: Mean longitude of lunar perigee
  longitudes%p = mod((334.3853 + 4069.0340 * T) - 0.0103 * (T**2), 360.0) * PI / 180.0
  ! n: Longitude of ascending node
  longitudes%N = mod((259.1568 - 1934.142 * T) + 0.0021 * (T**2), 360.0) * PI / 180.0
end subroutine astro_longitudes_init

!> Calculates the equilibrium phase argument for the given tidal
!! constituent constit and the astronomical longitudes and the reference time.
!! These formulas follow Table I.4 of Kowalik and Luick,
!! "Modern Theory and Practice of Tide Analysis and Tidal Power", 2019.
function eq_phase(constit, longitudes)
  character (len=2), intent(in) :: constit !> Name of constituent (e.g., M2).
  type(astro_longitudes), intent(in) :: longitudes   !> Mean longitudes calculated using astro_longitudes_init
  real, parameter :: PI = 4.0 * atan(1.0)  !> 3.14159... [nondim]
  real :: eq_phase                         !> The equilibrium phase argument for the constituent [rad].

  select case (constit)
    case ("M2")
      eq_phase = 2 * (longitudes%h - longitudes%s)
    case ("S2")
      eq_phase = 0.0
    case ("N2")
      eq_phase = (- 3 * longitudes%s + 2 * longitudes%h) + longitudes%p
    case ("K2")
      eq_phase = 2 * longitudes%h
    case ("K1")
      eq_phase = longitudes%h + PI / 2.0
    case ("O1")
      eq_phase = (- 2 * longitudes%s + longitudes%h) - PI / 2.0
    case ("P1")
      eq_phase = - longitudes%h - PI / 2.0
    case ("Q1")
      eq_phase = ((- 3 * longitudes%s + longitudes%h) + longitudes%p) - PI / 2.0
    case ("MF")
      eq_phase = 2 * longitudes%s
    case ("MM")
      eq_phase = longitudes%s - longitudes%p
    case default
      call MOM_error(FATAL, "eq_phase: unrecognized constituent")
  end select
end function eq_phase

!> Looks up angular frequencies for the main tidal constituents.
!! Values used here are from previous versions of MOM.
function tidal_frequency(constit)
  character (len=2), intent(in) :: constit !> Constituent to look up
  real :: tidal_frequency                  !> Angular frequency [rad s-1]

  select case (constit)
    case ("M2")
      tidal_frequency = 1.4051890e-4
    case ("S2")
      tidal_frequency = 1.4544410e-4
    case ("N2")
      tidal_frequency = 1.3787970e-4
    case ("K2")
      tidal_frequency = 1.4584234e-4
    case ("K1")
      tidal_frequency = 0.7292117e-4
    case ("O1")
      tidal_frequency = 0.6759774e-4
    case ("P1")
      tidal_frequency = 0.7252295e-4
    case ("Q1")
      tidal_frequency = 0.6495854e-4
    case ("MF")
      tidal_frequency = 0.053234e-4
    case ("MM")
      tidal_frequency = 0.026392e-4
    case default
      call MOM_error(FATAL, "tidal_frequency: unrecognized constituent")
  end select
end function tidal_frequency

!> Find amplitude (f) and phase (u) modulation of tidal constituents by the 18.6
!! year nodal cycle. Values here follow Table I.6 in Kowalik and Luick,
!! "Modern Theory and Practice of Tide Analysis and Tidal Power", 2019.
subroutine nodal_fu(constit, nodelon, fn, un)
  character (len=2), intent(in)  :: constit !> Tidal constituent to find modulation for.
  real,              intent(in)  :: nodelon !> Longitude of ascending node [rad], which
                                            !! can be calculated using astro_longitudes_init.
  real,              intent(out) :: fn      !> Amplitude modulation [nondim]
  real,              intent(out) :: un      !> Phase modulation [rad]

  real, parameter :: RADIANS = 4.0 * atan(1.0) / 180.0  !> Converts degrees to radians [nondim]

  select case (constit)
    case ("M2")
      fn = 1.0 - 0.037 * cos(nodelon)
      un = -2.1 * RADIANS * sin(nodelon)
    case ("S2")
      fn = 1.0  ! Solar S2 has no amplitude modulation.
      un = 0.0  ! S2 has no phase modulation.
    case ("N2")
      fn = 1.0 - 0.037 * cos(nodelon)
      un = -2.1 * RADIANS * sin(nodelon)
    case ("K2")
      fn = 1.024 + 0.286 * cos(nodelon)
      un = -17.7 * RADIANS * sin(nodelon)
    case ("K1")
      fn = 1.006 + 0.115 * cos(nodelon)
      un = -8.9 * RADIANS * sin(nodelon)
    case ("O1")
      fn = 1.009 + 0.187 * cos(nodelon)
      un = 10.8 * RADIANS * sin(nodelon)
    case ("P1")
      fn = 1.0  ! P1 has no amplitude modulation.
      un = 0.0  ! P1 has no phase modulation.
    case ("Q1")
      fn = 1.009 + 0.187 * cos(nodelon)
      un = 10.8 * RADIANS * sin(nodelon)
    case ("MF")
      fn = 1.043 + 0.414 * cos(nodelon)
      un = -23.7 * RADIANS * sin(nodelon)
    case ("MM")
      fn = 1.0 - 0.130 * cos(nodelon)
      un = 0.0  ! MM has no phase modulation.
    case default
      call MOM_error(FATAL, "nodal_fu: unrecognized constituent")
  end select

end subroutine nodal_fu

!> This subroutine allocates space for the static variables used
!! by this module.  The metrics may be effectively 0, 1, or 2-D arrays,
!! while fields like the background viscosities are 2-D arrays.
!! ALLOC is a macro defined in MOM_memory.h for allocate or nothing with
!! static memory.
subroutine tidal_forcing_init(Time, G, US, param_file, CS, HA_CS)
  type(time_type),        intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(tidal_forcing_CS), intent(inout) :: CS   !< Tidal forcing control structure
  type(harmonic_analysis_CS), optional, intent(out) :: HA_CS !< Control structure for harmonic analysis

  ! Local variables
  real, dimension(SZI_(G), SZJ_(G)) :: &
    phase, &          ! The phase of some tidal constituent [radians].
    lat_rad, lon_rad  ! Latitudes and longitudes of h-points [radians].
  real :: deg_to_rad  ! A conversion factor from degrees to radians [radian degree-1]
  real, dimension(MAX_CONSTITUENTS) :: freq_def ! Default frequency for each tidal constituent [rad s-1]
  real, dimension(MAX_CONSTITUENTS) :: phase0_def ! Default reference phase for each tidal constituent [rad]
  real, dimension(MAX_CONSTITUENTS) :: amp_def  ! Default amplitude for each tidal constituent [m]
  real, dimension(MAX_CONSTITUENTS) :: love_def ! Default love number for each constituent [nondim]
  integer, dimension(3) :: tide_ref_date !< Reference date (t = 0) for tidal forcing.
  integer, dimension(3) :: nodal_ref_date !< Reference date for calculating nodal modulation for tidal forcing.
  logical :: use_M2, use_S2, use_N2, use_K2, use_K1, use_O1, use_P1, use_Q1
  logical :: use_MF, use_MM
  logical :: tides      ! True if a tidal forcing is to be used.
  logical :: add_nodal_terms = .false.        !< If true, insert terms for the 18.6 year modulation when
                                              !! calculating tidal forcing.
  type(time_type) :: nodal_time               !< Model time to calculate nodal modulation for.
  type(astro_longitudes) :: nodal_longitudes  !< Solar and lunar longitudes for tidal forcing
  logical :: HA_ssh, HA_ubt, HA_vbt
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tidal_forcing" ! This module's name.
  character(len=128) :: mesg
  character(len=200) :: tidal_input_files(4*MAX_CONSTITUENTS)
  integer :: i, j, c, is, ie, js, je, isd, ied, jsd, jed, nc

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd; jed = G%jed

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TIDES", tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)

  if (.not.tides) return

  ! Set up the spatial structure functions for the diurnal, semidiurnal, and
  ! low-frequency tidal components.
  allocate(CS%sin_struct(isd:ied,jsd:jed,3), source=0.0)
  allocate(CS%cos_struct(isd:ied,jsd:jed,3), source=0.0)
  deg_to_rad = 4.0*ATAN(1.0)/180.0
  do j=js-1,je+1 ; do i=is-1,ie+1
    lat_rad(i,j) = G%geoLatT(i,j)*deg_to_rad
    lon_rad(i,j) = G%geoLonT(i,j)*deg_to_rad
  enddo ; enddo
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%sin_struct(i,j,1) = -sin(2.0*lat_rad(i,j)) * sin(lon_rad(i,j))
    CS%cos_struct(i,j,1) =  sin(2.0*lat_rad(i,j)) * cos(lon_rad(i,j))
    CS%sin_struct(i,j,2) = -cos(lat_rad(i,j))**2 * sin(2.0*lon_rad(i,j))
    CS%cos_struct(i,j,2) =  cos(lat_rad(i,j))**2 * cos(2.0*lon_rad(i,j))
    CS%sin_struct(i,j,3) =  0.0
    CS%cos_struct(i,j,3) = (0.5-1.5*sin(lat_rad(i,j))**2)
  enddo ; enddo

  call get_param(param_file, mdl, "TIDE_M2", use_M2, &
                 "If true, apply tidal momentum forcing at the M2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_S2", use_S2, &
                 "If true, apply tidal momentum forcing at the S2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_N2", use_N2, &
                 "If true, apply tidal momentum forcing at the N2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_K2", use_K2, &
                 "If true, apply tidal momentum forcing at the K2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_K1", use_K1, &
                 "If true, apply tidal momentum forcing at the K1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_O1", use_O1, &
                 "If true, apply tidal momentum forcing at the O1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_P1", use_P1, &
                 "If true, apply tidal momentum forcing at the P1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_Q1", use_Q1, &
                 "If true, apply tidal momentum forcing at the Q1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_MF", use_MF, &
                 "If true, apply tidal momentum forcing at the MF "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_MM", use_MM, &
                 "If true, apply tidal momentum forcing at the MM "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)

  ! Determine how many tidal components are to be used.
  nc = 0
  if (use_M2) nc=nc+1 ; if (use_S2) nc=nc+1
  if (use_N2) nc=nc+1 ; if (use_K2) nc=nc+1
  if (use_K1) nc=nc+1 ; if (use_O1) nc=nc+1
  if (use_P1) nc=nc+1 ; if (use_Q1) nc=nc+1
  if (use_MF) nc=nc+1 ; if (use_MM) nc=nc+1
  CS%nc = nc

  if (nc == 0) then
    call MOM_error(FATAL, "tidal_forcing_init: "// &
        "TIDES are defined, but no tidal constituents are used.")
    return
  endif

  call get_param(param_file, mdl, "TIDAL_SAL_FROM_FILE", CS%use_tidal_sal_file, &
                 "If true, read the tidal self-attraction and loading "//&
                 "from input files, specified by TIDAL_INPUT_FILE. "//&
                 "This is only used if TIDES is true.", default=.false.)
  call get_param(param_file, mdl, "USE_PREVIOUS_TIDES", CS%use_tidal_sal_prev, &
                 "If true, use the SAL from the previous iteration of the "//&
                 "tides to facilitate convergent iteration. "//&
                 "This is only used if TIDES is true.", default=.false.)
  if (CS%use_tidal_sal_prev) &
    call get_param(param_file, mdl, "SAL_SCALAR_VALUE", CS%sal_scalar, "The constant of "//&
                   "proportionality between self-attraction and loading (SAL) geopotential "//&
                   "anomaly and barotropic geopotential anomalies. This is only used if "//&
                   "SAL_SCALAR_APPROX is true or USE_PREVIOUS_TIDES is true.", default=0.0, &
                   units="m m-1", do_not_log=(.not.CS%use_tidal_sal_prev), &
                   old_name='TIDE_SAL_SCALAR_VALUE')

  if (nc > MAX_CONSTITUENTS) then
    write(mesg,'("Increase MAX_CONSTITUENTS in MOM_tidal_forcing.F90 to at least",I3, &
                &"to accommodate all the registered tidal constituents.")') nc
    call MOM_error(FATAL, "MOM_tidal_forcing"//mesg)
  endif

  do c=1,4*MAX_CONSTITUENTS ; tidal_input_files(c) = "" ; enddo

  if (CS%use_tidal_sal_file .or. CS%use_tidal_sal_prev) then
    call get_param(param_file, mdl, "TIDAL_INPUT_FILE", tidal_input_files, &
                   "A list of input files for tidal information.",         &
                   default="", fail_if_missing=.true.)
  endif

  call get_param(param_file, mdl, "TIDE_REF_DATE", tide_ref_date, &
                 "Year,month,day to use as reference date for tidal forcing. "//&
                 "If not specified, defaults to 0.", &
                 old_name="OBC_TIDE_REF_DATE", defaults=(/0, 0, 0/))

  call get_param(param_file, mdl, "TIDE_USE_EQ_PHASE", CS%use_eq_phase, &
                 "Correct phases by calculating equilibrium phase arguments for TIDE_REF_DATE. ", &
                 old_name="OBC_TIDE_ADD_EQ_PHASE", default=.false., fail_if_missing=.false.)

  if (sum(tide_ref_date) == 0) then  ! tide_ref_date defaults to 0.
    CS%time_ref = set_date(1, 1, 1, 0, 0, 0)
  else
    if (.not. CS%use_eq_phase) then
      ! Using a reference date but not using phase relative to equilibrium.
      ! This makes sense as long as either phases are overridden, or
      ! correctly simulating tidal phases is not desired.
      call MOM_mesg('Tidal phases will *not* be corrected with equilibrium arguments.')
    endif
    CS%time_ref = set_date(tide_ref_date(1), tide_ref_date(2), tide_ref_date(3), 0, 0, 0)
  endif

  ! Initialize reference time for tides and find relevant lunar and solar
  ! longitudes at the reference time.
  if (CS%use_eq_phase) call astro_longitudes_init(CS%time_ref, CS%tidal_longitudes)

  ! Set the parameters for all components that are in use.
  c=0
  if (use_M2) then
    c=c+1 ; CS%const_name(c) = "M2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.242334 ! Default amplitude in m.
  endif

  if (use_S2) then
    c=c+1 ; CS%const_name(c) = "S2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.112743 ! Default amplitude in m.
  endif

  if (use_N2) then
    c=c+1 ; CS%const_name(c) = "N2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.046397 ! Default amplitude in m.
  endif

  if (use_K2) then
    c=c+1 ; CS%const_name(c) = "K2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.030684 ! Default amplitude in m.
  endif

  if (use_K1) then
    c=c+1 ; CS%const_name(c) = "K1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.736 ; amp_def(c) = 0.141565 ! Default amplitude in m.
  endif

  if (use_O1) then
    c=c+1 ; CS%const_name(c) = "O1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.695 ; amp_def(c) = 0.100661 ! Default amplitude in m.
  endif

  if (use_P1) then
    c=c+1 ; CS%const_name(c) = "P1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.706 ; amp_def(c) = 0.046848 ! Default amplitude in m.
  endif

  if (use_Q1) then
    c=c+1 ; CS%const_name(c) = "Q1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.695 ; amp_def(c) = 0.019273 ! Default amplitude in m.
  endif

  if (use_MF) then
    c=c+1 ; CS%const_name(c) = "MF" ; CS%struct(c) = 3
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.042041 ! Default amplitude in m.
  endif

  if (use_MM) then
    c=c+1 ; CS%const_name(c) = "MM" ; CS%struct(c) = 3
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.022191 ! Default amplitude in m.
  endif

  ! Set defaults for all included constituents
  ! and things that can be set by functions
  do c=1,nc
    freq_def(c) = tidal_frequency(CS%const_name(c))
    love_def(c) = CS%love_no(c)
    CS%phase0(c) = 0.0
    if (CS%use_eq_phase) then
      phase0_def(c) = eq_phase(CS%const_name(c), CS%tidal_longitudes)
    else
      phase0_def(c) = 0.0
    endif
  enddo

  !  Parse the input file to potentially override the default values for the
  ! frequency, amplitude and initial phase of each constituent, and log the
  ! values that are actually used.
  do c=1,nc
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_FREQ", CS%freq(c), &
                   "Frequency of the "//trim(CS%const_name(c))//" tidal constituent. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true, or if OBC_TIDE_N_CONSTITUENTS > 0 and "//trim(CS%const_name(c))// &
                   " is in OBC_TIDE_CONSTITUENTS.", units="rad s-1", default=freq_def(c), &
                   scale=US%T_to_s)
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_AMP", CS%amp(c), &
                   "Amplitude of the "//trim(CS%const_name(c))//" tidal constituent. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true.", units="m", default=amp_def(c), scale=US%m_to_Z)
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_PHASE_T0", CS%phase0(c), &
                   "Phase of the "//trim(CS%const_name(c))//" tidal constituent at time 0. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true.", units="radians", default=phase0_def(c))
  enddo

  if (CS%use_tidal_sal_file) then
    allocate(CS%cosphasesal(isd:ied,jsd:jed,nc))
    allocate(CS%sinphasesal(isd:ied,jsd:jed,nc))
    allocate(CS%ampsal(isd:ied,jsd:jed,nc))
    do c=1,nc
      ! Read variables with names like PHASE_SAL_M2 and AMP_SAL_M2.
      call find_in_files(tidal_input_files, "PHASE_SAL_"//trim(CS%const_name(c)), phase, G)
      call find_in_files(tidal_input_files, "AMP_SAL_"//trim(CS%const_name(c)), CS%ampsal(:,:,c), &
                         G, scale=US%m_to_Z)
      call pass_var(phase,           G%domain,complete=.false.)
      call pass_var(CS%ampsal(:,:,c),G%domain,complete=.true.)
      do j=js-1,je+1 ; do i=is-1,ie+1
        CS%cosphasesal(i,j,c) = cos(phase(i,j)*deg_to_rad)
        CS%sinphasesal(i,j,c) = sin(phase(i,j)*deg_to_rad)
      enddo ; enddo
    enddo
  endif

  if (CS%use_tidal_sal_prev) then
    allocate(CS%cosphase_prev(isd:ied,jsd:jed,nc))
    allocate(CS%sinphase_prev(isd:ied,jsd:jed,nc))
    allocate(CS%amp_prev(isd:ied,jsd:jed,nc))
    do c=1,nc
      ! Read variables with names like PHASE_PREV_M2 and AMP_PREV_M2.
      call find_in_files(tidal_input_files, "PHASE_PREV_"//trim(CS%const_name(c)), phase, G)
      call find_in_files(tidal_input_files, "AMP_PREV_"//trim(CS%const_name(c)), CS%amp_prev(:,:,c), &
                         G, scale=US%m_to_Z)
      call pass_var(phase,             G%domain,complete=.false.)
      call pass_var(CS%amp_prev(:,:,c),G%domain,complete=.true.)
      do j=js-1,je+1 ; do i=is-1,ie+1
        CS%cosphase_prev(i,j,c) = cos(phase(i,j)*deg_to_rad)
        CS%sinphase_prev(i,j,c) = sin(phase(i,j)*deg_to_rad)
      enddo ; enddo
    enddo
  endif

  call get_param(param_file, mdl, "TIDE_ADD_NODAL", add_nodal_terms, &
                 "If true, include 18.6 year nodal modulation in the astronomical tidal forcing.", &
                 old_name="OBC_TIDE_ADD_NODAL", default=.false.)
  call get_param(param_file, mdl, "TIDE_NODAL_REF_DATE", nodal_ref_date, &
                 "Fixed reference date to use for nodal modulation of astronomical tidal forcing.", &
                 old_name="OBC_TIDE_REF_DATE", fail_if_missing=.false., defaults=(/0, 0, 0/))

  ! If the nodal correction is based on a different time, initialize that.
  ! Otherwise, it can use N from the time reference.
  if (add_nodal_terms) then
    if (sum(nodal_ref_date) /= 0) then
      ! A reference date was provided for the nodal correction
      nodal_time = set_date(nodal_ref_date(1), nodal_ref_date(2), nodal_ref_date(3))
      call astro_longitudes_init(nodal_time, nodal_longitudes)
    elseif (CS%use_eq_phase) then
      ! Astronomical longitudes were already calculated for use in equilibrium phases,
      ! so use nodal longitude from that.
      nodal_longitudes = CS%tidal_longitudes
    else
      ! Tidal reference time is a required parameter, so calculate the longitudes from that.
      call astro_longitudes_init(CS%time_ref, nodal_longitudes)
    endif
  endif

  allocate(CS%tide_fn(nc))
  allocate(CS%tide_un(nc))

  do c=1,nc
    ! Find nodal corrections if needed
    if (add_nodal_terms) then
      call nodal_fu(trim(CS%const_name(c)), nodal_longitudes%N, CS%tide_fn(c), CS%tide_un(c))
    else
      CS%tide_fn(c) = 1.0
      CS%tide_un(c) = 0.0
    endif
  enddo

  if (present(HA_CS)) then
    call HA_init(Time, US, param_file, CS%time_ref, CS%nc, CS%freq, CS%phase0, CS%const_name, &
                 CS%tide_fn, CS%tide_un, HA_CS)
    call get_param(param_file, mdl, "HA_SSH", HA_ssh, &
                   "If true, perform harmonic analysis of sea serface height.", default=.false.)
    if (HA_ssh) call HA_register('ssh', 'h', HA_CS)
    call get_param(param_file, mdl, "HA_UBT", HA_ubt, &
                   "If true, perform harmonic analysis of zonal barotropic velocity.", default=.false.)
    if (HA_ubt) call HA_register('ubt', 'u', HA_CS)
    call get_param(param_file, mdl, "HA_VBT", HA_vbt, &
                   "If true, perform harmonic analysis of meridional barotropic velocity.", default=.false.)
    if (HA_vbt) call HA_register('vbt', 'v', HA_CS)
  endif

  id_clock_tides = cpu_clock_id('(Ocean tides)', grain=CLOCK_MODULE)

end subroutine tidal_forcing_init

!> This subroutine finds a named variable in a list of files and reads its
!! values into a domain-decomposed 2-d array
subroutine find_in_files(filenames, varname, array, G, scale)
  character(len=*), dimension(:),   intent(in)  :: filenames !< The names of the files to search for the named variable
  character(len=*),                 intent(in)  :: varname   !< The name of the variable to read
  type(ocean_grid_type),            intent(in)  :: G         !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: array     !< The array to fill with the data [arbitrary]
  real,                   optional, intent(in)  :: scale     !< A factor by which to rescale the array to translate it
                                                             !! into its desired units [arbitrary]
  ! Local variables
  integer :: nf

  do nf=1,size(filenames)
    if (LEN_TRIM(filenames(nf)) == 0) cycle
    if (field_exists(filenames(nf), varname, MOM_domain=G%Domain)) then
      call MOM_read_data(filenames(nf), varname, array, G%Domain, scale=scale)
      return
    endif
  enddo

  do nf=size(filenames),1,-1
    if (file_exists(filenames(nf), G%Domain)) then
      call MOM_error(FATAL, "MOM_tidal_forcing.F90: Unable to find "// &
         trim(varname)//" in any of the tidal input files, last tried "// &
         trim(filenames(nf)))
    endif
  enddo

  call MOM_error(FATAL, "MOM_tidal_forcing.F90: Unable to find any of the "// &
                  "tidal input files, including "//trim(filenames(1)))

end subroutine find_in_files

!>   This subroutine calculates the geopotential anomalies that drive the tides,
!! including tidal self-attraction and loading from previous solutions.
subroutine calc_tidal_forcing(Time, e_tide_eq, e_tide_sal, G, US, CS)
  type(ocean_grid_type),            intent(in)  :: G          !< The ocean's grid structure.
  type(time_type),                  intent(in)  :: Time       !< The time for the caluculation.
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_eq  !< The geopotential height anomalies
                                                              !! due to the equilibrium tides [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_sal !< The geopotential height anomalies
                                                              !! due to the tidal SAL [Z ~> m].
  type(unit_scale_type),            intent(in)  :: US         !< A dimensional unit scaling type
  type(tidal_forcing_CS),           intent(in)  :: CS         !< The control structure returned by a
                                                              !! previous call to tidal_forcing_init.

  ! Local variables
  real :: now       ! The relative time compared with the tidal reference [T ~> s]
  real :: amp_cosomegat, amp_sinomegat ! The tidal amplitudes times the components of phase [Z ~> m]
  real :: cosomegat, sinomegat ! The components of the phase [nondim]
  integer :: i, j, c, m, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call cpu_clock_begin(id_clock_tides)

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    e_tide_eq(i,j) = 0.0
    e_tide_sal(i,j) = 0.0
  enddo ; enddo

  if (CS%nc == 0) then
    return
  endif

  now = US%s_to_T * time_type_to_real(Time - cs%time_ref)

  do c=1,CS%nc
    m = CS%struct(c)
    amp_cosomegat = CS%amp(c)*CS%love_no(c)*CS%tide_fn(c) * cos(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    amp_sinomegat = CS%amp(c)*CS%love_no(c)*CS%tide_fn(c) * sin(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_tide_eq(i,j) = e_tide_eq(i,j) + (amp_cosomegat*CS%cos_struct(i,j,m) + &
                                         amp_sinomegat*CS%sin_struct(i,j,m))
    enddo ; enddo
  enddo

  if (CS%use_tidal_sal_file) then ; do c=1,CS%nc
    cosomegat = CS%tide_fn(c) * cos(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    sinomegat = CS%tide_fn(c) * sin(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_tide_sal(i,j) = e_tide_sal(i,j) + CS%ampsal(i,j,c) * &
          (cosomegat*CS%cosphasesal(i,j,c) + sinomegat*CS%sinphasesal(i,j,c))
    enddo ; enddo
  enddo ; endif

  if (CS%use_tidal_sal_prev) then ; do c=1,CS%nc
    cosomegat = CS%tide_fn(c) * cos(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    sinomegat = CS%tide_fn(c) * sin(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_tide_sal(i,j) = e_tide_sal(i,j) - CS%sal_scalar * CS%amp_prev(i,j,c) * &
          (cosomegat*CS%cosphase_prev(i,j,c) + sinomegat*CS%sinphase_prev(i,j,c))
    enddo ; enddo
  enddo ; endif

  call cpu_clock_end(id_clock_tides)

end subroutine calc_tidal_forcing

!>   This subroutine functions the same as calc_tidal_forcing but outputs a field that combines
!! previously calculated self-attraction and loading (SAL) and tidal forcings, so that old answers
!! can be preserved bitwise before SAL is separated out as an individual module.
subroutine calc_tidal_forcing_legacy(Time, e_sal, e_sal_tide, e_tide_eq, e_tide_sal, G, US, CS)
  type(ocean_grid_type),            intent(in)  :: G          !< The ocean's grid structure.
  type(time_type),                  intent(in)  :: Time       !< The time for the caluculation.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: e_sal      !< The self-attraction and loading fields
                                                              !! calculated previously used to
                                                              !! initialized e_sal_tide [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_sal_tide !< The total geopotential height anomalies
                                                              !! due to both SAL and tidal forcings [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_eq  !< The geopotential height anomalies
                                                              !! due to the equilibrium tides [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: e_tide_sal !< The geopotential height anomalies
                                                              !! due to the tidal SAL [Z ~> m].
  type(unit_scale_type),            intent(in)  :: US         !< A dimensional unit scaling type
  type(tidal_forcing_CS),           intent(in)  :: CS         !< The control structure returned by a
                                                              !! previous call to tidal_forcing_init.

  ! Local variables
  real :: now       ! The relative time compared with the tidal reference [T ~> s]
  real :: amp_cosomegat, amp_sinomegat ! The tidal amplitudes times the components of phase [Z ~> m]
  real :: cosomegat, sinomegat ! The components of the phase [nondim]
  real :: amp_cossin ! A temporary field that adds cosines and sines [nondim]
  integer :: i, j, c, m, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call cpu_clock_begin(id_clock_tides)

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    e_sal_tide(i,j) = 0.0
    e_tide_eq(i,j) = 0.0
    e_tide_sal(i,j) = 0.0
  enddo ; enddo

  if (CS%nc == 0) then
    return
  endif

  now = US%s_to_T * time_type_to_real(Time - cs%time_ref)

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    e_sal_tide(i,j) = e_sal(i,j)
  enddo ; enddo

  do c=1,CS%nc
    m = CS%struct(c)
    amp_cosomegat = CS%amp(c)*CS%love_no(c)*CS%tide_fn(c) * cos(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    amp_sinomegat = CS%amp(c)*CS%love_no(c)*CS%tide_fn(c) * sin(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      amp_cossin = (amp_cosomegat*CS%cos_struct(i,j,m) + amp_sinomegat*CS%sin_struct(i,j,m))
      e_sal_tide(i,j) = e_sal_tide(i,j) + amp_cossin
      e_tide_eq(i,j) = e_tide_eq(i,j) + amp_cossin
    enddo ; enddo
  enddo

  if (CS%use_tidal_sal_file) then ; do c=1,CS%nc
    cosomegat = CS%tide_fn(c) * cos(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    sinomegat = CS%tide_fn(c) * sin(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      amp_cossin = CS%ampsal(i,j,c) &
        * (cosomegat*CS%cosphasesal(i,j,c) + sinomegat*CS%sinphasesal(i,j,c))
      e_sal_tide(i,j) = e_sal_tide(i,j) + amp_cossin
      e_tide_sal(i,j) = e_tide_sal(i,j) + amp_cossin
    enddo ; enddo
  enddo ; endif

  if (CS%use_tidal_sal_prev) then ; do c=1,CS%nc
    cosomegat = CS%tide_fn(c) * cos(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    sinomegat = CS%tide_fn(c) * sin(CS%freq(c)*now + (CS%phase0(c) + CS%tide_un(c)))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      amp_cossin = -CS%sal_scalar * CS%amp_prev(i,j,c) &
        * (cosomegat*CS%cosphase_prev(i,j,c) + sinomegat*CS%sinphase_prev(i,j,c))
      e_sal_tide(i,j) = e_sal_tide(i,j) + amp_cossin
      e_tide_sal(i,j) = e_tide_sal(i,j) + amp_cossin
    enddo ; enddo
  enddo ; endif
  call cpu_clock_end(id_clock_tides)

end subroutine calc_tidal_forcing_legacy

!> This subroutine deallocates memory associated with the tidal forcing module.
subroutine tidal_forcing_end(CS)
  type(tidal_forcing_CS), intent(inout) :: CS !< The control structure returned by a previous call
                                              !! to tidal_forcing_init; it is deallocated here.

  if (allocated(CS%sin_struct)) deallocate(CS%sin_struct)
  if (allocated(CS%cos_struct)) deallocate(CS%cos_struct)

  if (allocated(CS%cosphasesal)) deallocate(CS%cosphasesal)
  if (allocated(CS%sinphasesal)) deallocate(CS%sinphasesal)
  if (allocated(CS%ampsal))      deallocate(CS%ampsal)

  if (allocated(CS%cosphase_prev)) deallocate(CS%cosphase_prev)
  if (allocated(CS%sinphase_prev)) deallocate(CS%sinphase_prev)
  if (allocated(CS%amp_prev))      deallocate(CS%amp_prev)
end subroutine tidal_forcing_end

!> \namespace tidal_forcing
!!
!! \section section_tides Tidal forcing
!!
!! Code by Robert Hallberg, August 2005, based on C-code by Harper
!! Simmons, February, 2003, in turn based on code by Brian Arbic.
!!
!!   The main subroutine in this file calculates the total tidal
!! contribution to the geopotential, including self-attraction and
!! loading terms and the astronomical contributions.  All options
!! are selected with entries in a file that is parsed at run-time.
!! Overall tides are enabled with the run-time parameter 'TIDES=True'.
!! Tidal constituents must be individually enabled with lines like
!! 'TIDE_M2=True'.  This file has default values of amplitude,
!! frequency, Love number, and phase at time 0 for the Earth's M2,
!! S2, N2, K2, K1, O1, P1, Q1,  MF, and MM tidal constituents, but
!! the frequency, amplitude and phase ant time 0 for each constituent
!! can be changed at run time by setting variables like TIDE_M2_FREQ,
!! TIDE_M2_AMP and TIDE_M2_PHASE_T0 (for M2).
!!
!!   In addition, approaches to calculate self-attraction and loading
!! due to tides (harmonics of astronomical forcing frequencies)
!! are provided. <code>TIDAL_SAL_FROM_FILE</code> can be set to read the phase and
!! amplitude of the tidal SAL. <code>USE_PREVIOUS_TIDES</code> may be useful in
!! combination with the scalar approximation to iterate the SAL to
!! convergence (for details, see \cite Arbic2004). With
!! <code>TIDAL_SAL_FROM_FILE</code> or <code>USE_PREVIOUS_TIDES</code>, a list of input
!! files must be provided to describe each constituent's properties from
!! a previous solution. The online SAL calculations that are functions
!! of SSH (rather should be bottom pressure anmoaly), either a scalar
!! approximation or with spherical harmonic transforms, are located in
!! <code>MOM_self_attr_load</code>.
end module MOM_tidal_forcing
