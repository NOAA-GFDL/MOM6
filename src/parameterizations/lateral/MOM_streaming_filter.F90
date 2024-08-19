!> Streaming band-pass filter for detecting the instantaneous tidal signals in the simulation
module MOM_streaming_filter

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public Filt_register, Filt_accum

#include <MOM_memory.h>

!> The control structure for storing the filter infomation of a particular field
type, public :: Filter_CS ; private
  real                     :: a, &                   !< Parameter that determines the bandwidth [nondim]
                              om, &                  !< Target frequency of the filter [T-1 ~> s-1]
                              old_time = -1.0        !< The time of the previous accumulating step [T ~> s]
  real, allocatable        :: s1(:,:), &             !< Dummy variable [A]
                              u1(:,:)                !< Filtered data [A]
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
end type Filter_CS

contains

!> This subroutine registers each of the fields to be filtered.
subroutine Filt_register(a, om, CS)
  real,                      intent(in)    :: a      !< Parameter that determines the bandwidth [nondim]
  real,                      intent(in)    :: om     !< Target frequency of the filter [T-1 ~> s-1]
  type(Filter_CS),           intent(out)   :: CS     !< Control structure for the current field

  if (a  <= 0.0) call MOM_error(FATAL, "MOM_streaming_filter: bandwidth <= 0")
  if (om <= 0.0) call MOM_error(FATAL, "MOM_streaming_filter: target frequency <= 0")

  CS%a  = a
  CS%om = om

end subroutine Filt_register

!> This subroutine timesteps the filter equations. It takes model output u at the current time step as the input,
!! and returns tidal signal u1 as the output, which is the solution of a set of two ODEs (the filter equations).
subroutine Filt_accum(u, u1, Time, US, CS)
  real, dimension(:,:),          intent(in)  :: u    !< Input into the filter [A]
  real, dimension(:,:), pointer, intent(out) :: u1   !< Output of the filter [A]
  type(time_type),               intent(in)  :: Time !< The current model time
  type(unit_scale_type),         intent(in)  :: US   !< A dimensional unit scaling type
  type(Filter_CS),      target,  intent(inout) :: CS !< Control structure of the MOM_streaming_filter module

  ! Local variables
  real                      :: now, &                !< The current model time [T ~> s]
                               dt, &                 !< Time step size for the filter equations [T ~> s]
                               c1, c2                !< Coefficients for the filter equations [nondim]
  integer                   :: i, j, is, ie, js, je

  now = US%s_to_T * time_type_to_real(Time)

  ! Additional processing at the initial accumulating step
  if (CS%old_time < 0.0) then
    CS%old_time = now

    CS%is = LBOUND(u,1) ; is = CS%is
    CS%ie = UBOUND(u,1) ; ie = CS%ie
    CS%js = LBOUND(u,2) ; js = CS%js
    CS%je = UBOUND(u,2) ; je = CS%je

    allocate(CS%s1(is:ie,js:je), source=0.0)
    allocate(CS%u1(is:ie,js:je), source=0.0)

    do j=js,je ; do i=is,ie
      CS%s1(i,j)  = 0.0
      CS%u1(i,j)  = u(i,j)
    enddo ; enddo
  endif

  dt = now - CS%old_time
  CS%old_time = now

  is = CS%is ; ie = CS%ie ; js = CS%js ; je = CS%je

  ! Timestepping
  c1 = CS%om * dt
  c2 = 1.0 - CS%a * c1
  do j=js,je ; do i=is,ie
     CS%s1(i,j)  =  c1 * CS%u1(i,j) + CS%s1(i,j)
     CS%u1(i,j)  = -c1 * (CS%s1(i,j) - CS%a * u(i,j)) + c2 * CS%u1(i,j)
  enddo; enddo

  u1 => CS%u1

end subroutine Filt_accum

!> \namespace streaming_filter
!!
!! This module detects instantaneous tidal signals in the model output using a set of coupled ODEs (the filter
!! equations), given the target frequency (om) and the bandwidth parameter (a) of the filter. At each timestep,
!! the filter takes model output (u) as the input and returns a time series consisting of sinusoidal motions (u1)
!! near its target frequency. The filtered tidal signals can be used to parameterize frequency-dependent drag, or
!! to detide the model output. See Xu & Zaron (2024) for detail.
!!
!! Reference: Xu, C., & Zaron, E. D. (2024). Detecting instantaneous tidal signals in ocean models utilizing
!! streaming band-pass filters. Journal of Advances in Modeling Earth Systems. Under review.

end module MOM_streaming_filter

