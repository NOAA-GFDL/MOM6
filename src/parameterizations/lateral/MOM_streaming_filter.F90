!> Streaming band-pass filter for detecting the instantaneous tidal signals in the simulation
!!
!! Dec 5, 2024: Major revision.
!!
!! The filters and their target frequencies are no longer hard-coded. Instead, up to 10 filters
!! with tidal frequencies as their target frequencies and an unspecified number of filters with
!! arbitrary target frequencies can be turned on. The filter names are specified in MOM_input
!! and must consist of two letters/numbers. If a filter name is the same as the name of a tidal
!! constituent, then the corresponding tidal frequency will be used as its target frequency.
!! Otherwise, the user must provide the target frequency. In either case, the target frequency
!! is specified by "TIDE_${FILTER_NAME}_FREQ" in MOM_input.

module MOM_streaming_filter

use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE, FATAL
use MOM_file_parser,   only : get_param, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_tidal_forcing, only : tidal_frequency
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public Filt_register, Filt_accum

#include <MOM_memory.h>

!> Control structure for the MOM_streaming_filter module
type, public :: Filter_CS ; private
  integer :: nf                       !< Number of filters to be used in the simulation
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
  character(len=2), allocatable, dimension(:) :: filter_names !< Names of filters
  real, allocatable, dimension(:)     :: filter_omega !< Target frequencies of filters [T-1 ~> s-1]
  real, allocatable, dimension(:)     :: filter_alpha !< Bandwidth parameters of filters [nondim]
  real, allocatable, dimension(:,:,:) :: s1, &        !< Dummy variable [A]
                                         u1           !< Filtered data [A]
  real :: old_time = -1.0             !< The time of the previous accumulating step [T ~> s]
end type Filter_CS

contains

!> This subroutine initializes the filters given the number of filters and the grid
subroutine Filt_register(nf, grid, HI, US, param_file, CS)
  integer,               intent(in)   :: nf           !< Number of filters to be used in the simulation
  character(len=*),      intent(in)   :: grid         !< Horizontal grid location: h, u, or v
  type(hor_index_type),  intent(in)   :: HI           !< Horizontal index type structure
  type(unit_scale_type), intent(in)   :: US           !< A dimensional unit scaling type
  type(param_file_type), intent(in)   :: param_file   !< A structure to parse for run-time parameters
  type(Filter_CS),       intent(out)  :: CS           !< Control structure for MOM_streaming_filter

  ! Local variables
  character(len=40)  :: mdl = "MOM_streaming_filter"  !< This module's name
  character(len=50)  :: filter_name_str               !< List of filters to be registered
  character(len=200) :: mesg
  integer :: c

  CS%nf = nf

  select case (trim(grid))
    case ('h')
      CS%is = HI%isd  ; CS%ie = HI%ied  ; CS%js = HI%jsd  ; CS%je = HI%jed
    case ('u')
      CS%is = HI%IsdB ; CS%ie = HI%IedB ; CS%js = HI%jsd  ; CS%je = HI%jed
    case ('v')
      CS%is = HI%isd  ; CS%ie = HI%ied  ; CS%js = HI%JsdB ; CS%je = HI%JedB
    case default
      call MOM_error(FATAL, "MOM_streaming_filter: horizontal grid not supported")
  end select

  allocate(CS%s1(CS%is:CS%ie,CS%js:CS%je,nf)) ; CS%s1(:,:,:) = 0.0
  allocate(CS%u1(CS%is:CS%ie,CS%js:CS%je,nf)) ; CS%u1(:,:,:) = 0.0

  call get_param(param_file, mdl, "FILTER_NAMES", filter_name_str, &
                 "Names of streaming band-pass filters to be used in the simulation.", &
                 fail_if_missing=.true.)
  allocate(CS%filter_names(nf))
  allocate(CS%filter_omega(nf))
  allocate(CS%filter_alpha(nf))
  read(filter_name_str, *) CS%filter_names

  do c=1,nf
    ! If filter_name_str consists of tidal constituents, use tidal frequencies.
    call get_param(param_file, mdl, "TIDE_"//trim(CS%filter_names(c))//"_FREQ", CS%filter_omega(c), &
                   "Target frequency of the "//trim(CS%filter_names(c))//" filter. "//&
                   "This is used if USE_FILTER is true and "//trim(CS%filter_names(c))//&
                   " is in FILTER_NAMES, even if TIDES and TIDE_"//trim(CS%filter_names(c))//&
                   " are false.", units="s-1", &
                   default=tidal_frequency(trim(CS%filter_names(c))), scale=US%T_to_s)
    call get_param(param_file, mdl, "FILTER_"//trim(CS%filter_names(c))//"_ALPHA", CS%filter_alpha(c), &
                   "Bandwidth parameter of the "//trim(CS%filter_names(c))//" filter. "//&
                   "Must be positive.", units="nondim", fail_if_missing=.true.)
    if (CS%filter_omega(c)<=0.0) call MOM_error(FATAL, "MOM_streaming_filter: target frequency <= 0")
    if (CS%filter_alpha(c)<=0.0) call MOM_error(FATAL, "MOM_streaming_filter: bandwidth <= 0")

    write(mesg,*) "MOM_streaming_filter: ", trim(CS%filter_names(c)), &
                  " filter registered, target frequency = ", CS%filter_omega(c), &
                  ", bandwidth = ", CS%filter_alpha(c)
    call MOM_error(NOTE, trim(mesg))
  enddo

end subroutine Filt_register

!> This subroutine timesteps the filter equations. It takes model output u at the current time step as the input,
!! and returns tidal signal u1 as the output, which is the solution of a set of two ODEs (the filter equations).
subroutine Filt_accum(u, u1, Time, US, CS)
  real, dimension(:,:,:), pointer, intent(out)   :: u1   !< Output of the filter [A]
  type(time_type),                 intent(in)    :: Time !< The current model time
  type(unit_scale_type),           intent(in)    :: US   !< A dimensional unit scaling type
  type(Filter_CS),        target,  intent(inout) :: CS   !< Control structure of the MOM_streaming_filter module
  real, dimension(CS%is:CS%ie,CS%js:CS%je), intent(in) :: u !< Input into the filter [A]

  ! Local variables
  real    :: now, &              !< The current model time [T ~> s]
             dt, &               !< Time step size for the filter equations [T ~> s]
             c1, c2              !< Coefficients for the filter equations [nondim]
  integer :: i, j, k, is, ie, js, je

  now = US%s_to_T * time_type_to_real(Time)
  is = CS%is ; ie = CS%ie ; js = CS%js ; je = CS%je

  ! Initialize u1
  if (CS%old_time < 0.0) then
    CS%old_time = now
    do k=1,CS%nf
      CS%u1(:,:,k) = u(:,:)
    enddo
  endif

  dt = now - CS%old_time
  CS%old_time = now

  ! Timestepping
  do k=1,CS%nf
    c1 = CS%filter_omega(k) * dt
    c2 = 1.0 - CS%filter_alpha(k) * c1

    do j=js,je ; do i=is,ie
      CS%s1(i,j,k) =  c1 *  CS%u1(i,j,k) + CS%s1(i,j,k)
      CS%u1(i,j,k) = -c1 * (CS%s1(i,j,k) - CS%filter_alpha(k) * u(i,j)) + c2 * CS%u1(i,j,k)
    enddo; enddo
  enddo
  u1 => CS%u1

end subroutine Filt_accum

!> \namespace streaming_filter
!!
!! This module detects instantaneous tidal signals in the model output using a set of coupled ODEs (the filter
!! equations), given the target frequency (filter_omega) and the bandwidth parameter (filter_alpha) of the filter.
!! At each timestep, the filter takes model output (u) as the input and returns a time series (u1) consisting of
!! sinusoidal motions near its target frequency. The filtered tidal signals can be used to parameterize
!! frequency-dependent drag, or to detide the model output. See Xu & Zaron (2024) for detail.
!!
!! Reference: Xu, C., & Zaron, E. D. (2024). Detecting instantaneous tidal signals in ocean models utilizing
!! streaming band-pass filters. Journal of Advances in Modeling Earth Systems, 16, e2024MS004319.
!! https://doi.org/10.1029/2024MS004319

end module MOM_streaming_filter

