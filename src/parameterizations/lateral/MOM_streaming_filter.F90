!> Streaming band-pass filter for detecting the instantaneous tidal signals in the simulation
module MOM_streaming_filter

use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public Filt_init, Filt_register, Filt_accum

#include <MOM_memory.h>

!> The private control structure for storing the filter infomation of a particular field
type, private :: Filt_type
  character(len=16)        :: key = 'none'           !< Name of the current field
  real                     :: a, &                   !< Parameter that determines the bandwidth [nondim]
                              om, &                  !< Target frequency of the filter [T-1 ~> s-1]
                              old_time = -1.0        !< The time of the previous accumulating step [T ~> s]
  integer                  :: is, ie, js, je         !< Lower and upper bounds of input data
  real, allocatable        :: s1(:,:), &             !< Dummy variable [A]
                              u1(:,:)                !< Filtered data [A]
end type Filt_type

!> A linked list of control structures that store the filter infomation of different fields
type, private :: Filt_node
  type(Filt_type)          :: this                   !< Control structure of the current field in the list
  type(Filt_node), pointer :: next                   !< The list of other fields
end type Filt_node

!> The public control structure of the MOM_streaming_filter module
type, public :: streaming_filter_CS ; private
  integer                  :: length                 !< Number of fields to be filtered
  type(Filt_node), pointer :: list => NULL()         !< A linked list for storing filter info of different fields
end type streaming_filter_CS

contains

!> This subroutine initializes CS%list.
subroutine Filt_init(CS)
  type(streaming_filter_CS), intent(out)   :: CS     !< Control structure of the MOM_streaming_filter module

  ! Local variable
  type(Filt_type) :: ha1                             !< A temporary, null field used for initializing CS%list

  allocate(CS%list)
  CS%list%this  =  ha1
  nullify(CS%list%next)
  CS%length     =  0

end subroutine Filt_init

!> This subroutine registers each of the fields to be filtered.
subroutine Filt_register(key, a, om, CS)
  character(len=*),          intent(in)    :: key    !< Name of the current field
  real,                      intent(in)    :: a      !< Parameter that determines the bandwidth [nondim]
  real,                      intent(in)    :: om     !< Target frequency of the filter [T-1 ~> s-1]
  type(streaming_filter_CS), intent(inout) :: CS     !< Control structure of the MOM_streaming_filter module

  ! Local variables
  type(Filt_type)          :: ha1                    !< Control structure for the current field
  type(Filt_node), pointer :: tmp                    !< A temporary list to hold the current field
  character(len=128)       :: mesg

  if (a <= 0.0) then
    write(mesg, *) "MOM_streaming_filter: bandwidth <= 0, key ", trim(key), " not registered."
    call MOM_error(NOTE, trim(mesg))
    return                                           !< Setting a <= 0.0 turns off the filter
  endif
  if (om <= 0.0) then
    write(mesg, *) "MOM_streaming_filter: target frequency <= 0, key ", trim(key), " not registered."
    call MOM_error(NOTE, trim(mesg))
    return                                           !< Setting om <= 0.0 turns off the filter
  endif

  ha1%key   =  trim(key)
  ha1%a     =  a
  ha1%om    =  om
  allocate(tmp)
  tmp%this  =  ha1
  tmp%next  => CS%list
  CS%list   => tmp
  CS%length =  CS%length + 1

end subroutine Filt_register

!> This subroutine timesteps the filter equations. It takes model output u at the current time step as the input,
!! and returns tidal signal u1 as the output, which is the solution of a set of two ODEs (the filter equations).
subroutine Filt_accum(key, u, u1, Time, US, CS)
  character(len=*),              intent(in)  :: key  !< Name of the current field
  real, dimension(:,:),          intent(in)  :: u    !< Input into the filter [A]
  real, dimension(:,:), pointer, intent(out) :: u1   !< Output of the filter [A]
  type(time_type),               intent(in)  :: Time !< The current model time
  type(unit_scale_type),         intent(in)  :: US   !< A dimensional unit scaling type
  type(streaming_filter_CS),     intent(in)  :: CS   !< Control structure of the MOM_streaming_filter module

  ! Local variables
  type(Filt_type), pointer  :: ha1
  type(Filt_node), pointer  :: tmp
  real                      :: now, &                !< The current model time [T ~> s]
                               dt, &                 !< Time step size for the filter equations [T ~> s]
                               c1, c2                !< Coefficients for the filter equations [nondim]
  real, allocatable, target :: u0(:,:)               !< Output (zeros) when a field is not registered [A]
  integer                   :: i, j, k, is, ie, js, je
  character(len=128)        :: mesg

  ! Exit the accumulator and return zeros if no field is registered
  if (CS%length == 0) then
    is = LBOUND(u,1) ; ie = UBOUND(u,1) ; js = LBOUND(u,2) ; je = UBOUND(u,2)
    allocate(u0(is:ie,js:je), source=0.0) ; u0(:,:) = 0.0 ; u1 => u0 ; return
  endif

  ! Loop through the full list to find the current field
  tmp => CS%list
  do k=1,CS%length
    ha1 => tmp%this
    if (trim(key) == trim(ha1%key)) exit
    tmp => tmp%next

    ! Exit the accumulator and return zeros if the field is not registered
    if (k == CS%length) then
      is = LBOUND(u,1) ; ie = UBOUND(u,1) ; js = LBOUND(u,2) ; je = UBOUND(u,2)
      allocate(u0(is:ie,js:je), source=0.0) ; u0(:,:) = 0.0 ; u1 => u0 ; return
    endif
  enddo

  now = US%s_to_T * time_type_to_real(Time)

  ! Additional processing at the initial accumulating step
  if (ha1%old_time < 0.0) then
    ha1%old_time = now

    write(mesg,*) "MOM_streaming_filter: initializing filter equations, key = ", trim(ha1%key)
    call MOM_error(NOTE, trim(mesg))

    ha1%is = LBOUND(u,1) ; is = ha1%is
    ha1%ie = UBOUND(u,1) ; ie = ha1%ie
    ha1%js = LBOUND(u,2) ; js = ha1%js
    ha1%je = UBOUND(u,2) ; je = ha1%je

    allocate(ha1%s1(is:ie,js:je), source=0.0)
    allocate(ha1%u1(is:ie,js:je), source=0.0)
    do j=js,je ; do i=is,ie
      ha1%s1(i,j)  = 0.0
      ha1%u1(i,j)  = u(i,j)
    end do ; end do
  endif

  dt = now - ha1%old_time
  ha1%old_time = now

  is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  ! Timestepping
  c1 = ha1%om * dt
  c2 = 1.0 - ha1%a * c1
  do j=js,je ; do i=is,ie
     ha1%s1(i,j)  =  c1 * (ha1%u1(i,j) + (ha1%a**2) * u(i,j)) + c2 * ha1%s1(i,j)
     ha1%u1(i,j)  = -c1 * (ha1%s1(i,j) + (-2*ha1%a) * u(i,j)) + c2 * ha1%u1(i,j)
  enddo; enddo

  u1 => ha1%u1

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

