!> Inline harmonic analysis (conventional)
module MOM_HA

! This module computes the harmonic constants which can be used to reconstruct the tidal elevation (or other 
! fields) through SSH = F * x, where F is an nt-by-2*nc matrix (nt is the number of time steps and nc is the 
! number of tidal constituents) containing the cosine/sine functions for each frequency evaluated at each time 
! step, and x is a 2*nc-by-1 vector containing the constant coefficients of the sine/cosine for each constituent 
! (i.e., the harmonic constants). At each grid point, the harmonic constants are computed using least squares, 
!     x = (F' * F)^{-1} * (F' * SSH_in),
! where the prime denotes matrix transpose, and SSH_in is the sea surface height (or other fields) determined by 
! the model. The dot products (F' * F) and (F' * SSH_in) are computed by accumulating the sums as the model is 
! running and stored in the arrays FtF and FtSSH, respectively. The FtF matrix is inverted as needed before 
! computing and writing out the harmonic constants.
!
! The start and end times of harmonic analysis (HA), in number of days after the start of the run segment, are 
! specified by HA_start_time and HA_end_time, respectively, in MOM_input. HA will start on the first time step 
! if HA_start_time <= 0. HA will not be performed if HA_start_time >= HA_end_time, or if HA_end_time <= 0 or is 
! greater than the length of the run segment. This module does not have restarting capability.
!
! Ed Zaron and William Xu (chengzhu.xu@oregonstate.edu), April 2024.

use MOM_time_manager,  only : time_type, time_type_to_real, get_date, operator(-)
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_file_parser,   only : param_file_type, get_param
use MOM_tidal_forcing, only : tidal_forcing_CS
use MOM_io,            only : MOM_infra_file, vardesc, MOM_field, &
                              var_desc, create_MOM_file, SINGLE_FILE, MOM_write_field
use MOM_error_handler, only : MOM_mesg, MOM_error, WARNING

implicit none ; private

public HA_init, HA_register, HA_accum

#include <MOM_memory.h>

integer, parameter              :: HA_strlen = 255
type(ocean_grid_type)           :: HA_G
type(tidal_forcing_CS), pointer :: HA_CS => NULL() !< Control structure for tides

!> The control structure for storing the HA info of a particular field
type, public :: HA_type
  character(len=HA_strlen) :: key = "none"         !< Name of the field on which HA is performed
  logical :: isfirst = .true.                      !< If the current step is the first accumulating step
  logical :: ioready = .true.                      !< Perform HA and write results if true
  integer :: is, ie, js, je                        !< Lower and upper bounds of input data
  real, allocatable :: ref(:,:)                    !< The initial field
  real, allocatable :: FtF(:,:)                    !< Accumulator of (F' * F)
  real, allocatable :: FtSSH(:,:,:)                !< Accumulator of (F' * SSH_in)
end type HA_type

!> A linked list of control structures that store the HA info of different fields
type :: HA_node
  type(HA_type) :: h
  type(HA_node), pointer :: next
end type HA_node

integer :: HA_length                               !< Number of fields of which HA is to be performed
real    :: HA_start_time                           !< Start time of harmonic analysis
real    :: HA_end_time                             !< End time of harmonic analysis
type(HA_node), pointer   :: HA_list                !< A linked list for storing the HA info
character(len=HA_strlen) :: HA_path                !< Path to directory where output will be written

contains

!> This subroutine sets private global parameters for the MOM_HA module.
!! THIS MUST BE CALLED AFTER tidal_forcing_init.
subroutine HA_init(Time, G, US, param_file, CS)
  type(time_type),        intent(in) :: Time       !< The current model time
  type(ocean_grid_type),  intent(in) :: G          !< The ocean's grid structure
  type(unit_scale_type),  intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type),  intent(in) :: param_file !< A structure to parse for run-time parameters
  type(tidal_forcing_CS), pointer, intent(in) :: CS         !< Tidal forcing control structure

  ! Local variables
  real :: now
  type(HA_type) :: ha1
  character(len=40) :: mdl="MOM_HA"                !< This module's name
  character(len=128) :: mesg

  HA_G  =  G                                       !< Make a copy of G passed in from MOM_tidal_forcing
  HA_CS => CS                                      !< Associate CS passed in from MOM_tidal_forcing

  ! Initialize the linked list for storing the HA info
  allocate(HA_list)
  ha1%key   = "none"
  HA_list%h = ha1
  HA_length = 0
  nullify(HA_list%next)

  now = US%s_to_T * time_type_to_real(Time - CS%time_ref)

  ! Determine the start time of harmonic analysis
  call get_param(param_file, mdl, "HA_start_time", HA_start_time, &
                "Start time of harmonic analysis, in number of days after the start of run segment. "//&
                "If HA_start_time <= 0, then harmonic analysis will start on the first time step. "//&
                "If HA_start_time >= HA_end_time, then harmonic analysis will not be performed.", &
                units="days", default=0.0, scale=US%T_to_s, do_not_log=.True., fail_if_missing=.false.)
  HA_start_time = now + HA_start_time * 8.64e4

  ! Determine the end time of harmonic analysis
  call get_param(param_file, mdl, "HA_end_time", HA_end_time, &
                "End time of harmonic analysis, in number of days after the start of run segment. "//&
                "If HA_end_time <= 0, then harmonic analysis will not be performed. "//&
                "Also, HA_end_time needs to be smaller than the length of the run segment, "//&
                "in order for harmonic analysis to finish before the simulation ends.", &
                units="days", default=0.0, scale=US%T_to_s, do_not_log=.True., fail_if_missing=.false.)
  HA_end_time = now + HA_end_time * 8.64e4

  write(mesg,*) "HAmod: run segment starting time = ", now/8.64e4
  call MOM_error(WARNING, trim(mesg))
  write(mesg,*) "HAmod: HA starting time = ", HA_start_time/8.64e4
  call MOM_error(WARNING, trim(mesg))
  write(mesg,*) "HAmod: HA end time = ", HA_end_time/8.64e4
  call MOM_error(WARNING, trim(mesg))

 ! Set path to directory where output will be written
  call get_param(param_file, mdl, "HA_path", HA_path, &
                 "Path to output files for runtime harmonic analysis.", &
                 default="./", fail_if_missing=.false.)

end subroutine HA_init

!> This subroutine registers each of the fields on which HA is to be performed.
subroutine HA_register(key)
  character(len=*), intent(in) :: key

  ! Local variables
  type(HA_type) :: ha1
  type(HA_node), pointer :: tmp
  character(len=128) :: mesg

  write(mesg,*) "HAmod: register key = ", trim(key)
  call MOM_error(WARNING, trim(mesg))

  allocate(tmp)
  ha1%key   =  trim(key)
  tmp%h     =  ha1
  tmp%next  => HA_list
  HA_list   => tmp
  HA_length =  HA_length + 1

end subroutine HA_register

!> This subroutine accumulates the temporal basis functions in FtF and FtSSH and then calls HA_write to compute 
!! harmonic constants and write results. The tidal constituents are those used in MOM_tidal_forcing, plus the 
!! mean (of zero frequency).
subroutine HA_accum(key, data, Time, US)
  character(len=*),       intent(in) :: key
  real, dimension(:,:),   intent(in) :: data
  type(time_type),        intent(in) :: Time       !< The current model time
  type(unit_scale_type),  intent(in) :: US         !< A dimensional unit scaling type

  ! Local variables
  type(HA_type), pointer :: ha1
  type(HA_node), pointer :: tmp
  character(len=128) :: mesg
  integer :: nc, i, j, k, c, icos, isin, cc, iccos, issin, is, ie, js, je
  real :: now, cosomegat, sinomegat, ccosomegat, ssinomegat

  ! Do not perform harmonic analysis in the following cases
  if (HA_end_time <= 0) return
  if (HA_start_time >= HA_end_time) return

  ! Loop through the full list to find the current field
  tmp => HA_list
  do k=1,HA_length
    ha1 => tmp%h
    if (trim(key) == trim(ha1%key)) exit
    tmp => tmp%next
    if (k == HA_length) return                     !< Do not perform HA of a field that is not registered.
  enddo

  nc = HA_CS%nc

  ! Additional processing at the initial accumulating step
  if (ha1%isfirst) then
    ha1%isfirst = .false.

    write(mesg,*) "HAmod: initializing accumulator, key = ", trim(ha1%key)
    call MOM_error(WARNING, trim(mesg))

    ! Get the lower and upper bounds of input data
    ha1%is = LBOUND(data,1) ; is = ha1%is
    ha1%ie = UBOUND(data,1) ; ie = ha1%ie
    ha1%js = LBOUND(data,2) ; js = ha1%js
    ha1%je = UBOUND(data,2) ; je = ha1%je

    allocate(ha1%ref(is:ie,js:je), source=0.0)
    allocate(ha1%FtF(2*nc+1,2*nc+1), source=0.0)
    allocate(ha1%FtSSH(is:ie,js:je,2*nc+1), source=0.0)

    do j=js,je ; do i=is,ie
      ha1%ref(i,j) = data(i,j)
    enddo ; enddo
  else
    is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je 
  endif

  now = US%s_to_T * time_type_to_real(Time - HA_CS%time_ref)

  if (now > HA_start_time) then

    ! Accumulate FtF and FtSSH
    if (now < HA_end_time) then
      ha1%FtF(1,1) = ha1%FtF(1,1) + 1.0            !< For the zero frequency
      do c=1,nc
        icos = 2*c
        isin = 2*c+1  
        cosomegat = cos(HA_CS%freq(c) * now + HA_CS%phase0(c))
        sinomegat = sin(HA_CS%freq(c) * now + HA_CS%phase0(c))
        ha1%FtF(icos,1) = ha1%FtF(icos,1) + cosomegat
        ha1%FtF(isin,1) = ha1%FtF(isin,1) + sinomegat
        ha1%FtF(1,icos) = ha1%FtF(icos,1)
        ha1%FtF(1,isin) = ha1%FtF(isin,1)
        do cc=c,nc
          iccos = 2*cc  
          issin = 2*cc+1  
          ccosomegat = cos(HA_CS%freq(cc) * now + HA_CS%phase0(cc))
          ssinomegat = sin(HA_CS%freq(cc) * now + HA_CS%phase0(cc))
          ha1%FtF(icos,iccos) = ha1%FtF(icos,iccos) + cosomegat * ccosomegat
          ha1%FtF(icos,issin) = ha1%FtF(icos,issin) + cosomegat * ssinomegat
          ha1%FtF(isin,iccos) = ha1%FtF(isin,iccos) + sinomegat * ccosomegat
          ha1%FtF(isin,issin) = ha1%FtF(isin,issin) + sinomegat * ssinomegat
        enddo ! cc=c,nc
        do j=js,je ; do i=is,ie
          ha1%FtSSH(i,j,1)    = ha1%FtSSH(i,j,1)    + (data(i,j) - ha1%ref(i,j))
          ha1%FtSSH(i,j,icos) = ha1%FtSSH(i,j,icos) + (data(i,j) - ha1%ref(i,j)) * cosomegat
          ha1%FtSSH(i,j,isin) = ha1%FtSSH(i,j,isin) + (data(i,j) - ha1%ref(i,j)) * sinomegat
        enddo ; enddo
      enddo ! c=1,nc

    ! Compute harmonic constants and write output
    elseif (ha1%ioready) then
      ha1%ioready = .false.
      call HA_write(ha1, Time)
    endif ! (now < HA_end_time)

  endif ! (now > HA_start_time)

end subroutine HA_accum

subroutine HA_write(ha1, Time)
  type(HA_type), pointer, intent(in)  :: ha1
  type(time_type),        intent(in)  :: Time

  ! Local variables
  real, dimension(:,:,:), allocatable :: FtSSHw
  integer :: year, month, day, hour, minute, second
  integer :: nc, k, is, ie, js, je

  character(len=HA_strlen)     :: filename         !< Output file name
  type(MOM_infra_file)         :: cdf              !< The file handle for output harmonic constants
  type(vardesc),   allocatable :: cdf_vars(:)      !< Output variable names
  type(MOM_field), allocatable :: cdf_fields(:)    !< Field type variables for the output fields

  nc = HA_CS%nc ; is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  allocate(FtSSHw(is:ie,js:je,2*nc+1), source=0.0)

  ! Compute the harmonic coefficients
  call HA_solver(ha1, FtSSHw)

  ! Output file name
  call get_date(Time, year, month, day, hour, minute, second)
  write(filename, '(a,"HA_",a,i0.4,i0.2,i0.2,".nc")') &
      trim(HA_path), trim(ha1%key), year, month, day

  allocate(cdf_vars(2*nc+1))
  allocate(cdf_fields(2*nc+1))

  ! Variable names
  cdf_vars(1) = var_desc("z0", "m" ,"mean value", 'h', '1')
  do k=1,nc
    cdf_vars(2*k  ) = var_desc(trim(HA_CS%const_name(k))//"cos", "m", "cosine coefficient", 'h', '1')
    cdf_vars(2*k+1) = var_desc(trim(HA_CS%const_name(k))//"sin", "m", "sine coefficient",   'h', '1')
  enddo

  ! Create output file
  call create_MOM_file(cdf, trim(filename), cdf_vars, &
                       2*nc+1, cdf_fields, SINGLE_FILE, 86400.0, G=HA_G)

  ! Write data
  call MOM_write_field(cdf, cdf_fields(1), HA_G%domain, FtSSHw(:,:,1), 0.0)
  do k=1,nc
    call MOM_write_field(cdf, cdf_fields(2*k  ), HA_G%domain, FtSSHw(:,:,2*k  ), 0.0)
    call MOM_write_field(cdf, cdf_fields(2*k+1), HA_G%domain, FtSSHw(:,:,2*k+1), 0.0)
  enddo

  call cdf%flush()
  deallocate(cdf_vars)
  deallocate(cdf_fields)

end subroutine HA_write

!> This subroutine computes the harmonic constants (stored in FtSSHw) using the dot products of the temporal 
!! basis functions accumulated in FtF, and the dot products of the SSH (or other fields) with the temporal basis
!! functions accumulated in FtSSH. The system is solved by Cholesky decomposition.
subroutine HA_solver(ha1, FtSSHw)
  type(HA_type), pointer,              intent(in)  :: ha1
  real, dimension(:,:,:), allocatable, intent(out) :: FtSSHw

  ! Local variables
  real, dimension(:,:), allocatable :: tmp
  real, dimension(:,:), allocatable :: FtFw
  integer :: nc, k, l, is, ie, js, je

  nc = HA_CS%nc ; is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  allocate(tmp(is:ie,js:je), source=0.0)
  allocate(FtFw(1:2*nc+1,1:2*nc+1), source=0.0)
  allocate(FtSSHw(is:ie,js:je,2*nc+1), source=0.0)

  FtFw(:,:) = 0.0
  do l=1,2*nc+1
    FtFw(l,l) = sqrt(ha1%FtF(l,l) - dot_product(FtFw(l,1:l-1), FtFw(l,1:l-1)))
    do k=l+1,2*nc+1
      FtFw(k,l) = (ha1%FtF(k,l) - dot_product(FtFw(k,1:l-1), FtFw(l,1:l-1))) / FtFw(l,l)
    enddo
  enddo

  FtSSHw(:,:,:) = ha1%FtSSH(:,:,:)
  do k=1,2*nc+1
    tmp(:,:) = 0.0
    do l=1,k-1
      tmp(:,:) = tmp(:,:) + FtFw(k,l) * FtSSHw(:,:,l)
    enddo
    FtSSHw(:,:,k) = (FtSSHw(:,:,k) - tmp(:,:)) / FtFw(k,k)
  enddo
  do k=2*nc+1,1,-1
    tmp(:,:) = 0.0
    do l=k+1,2*nc+1
      tmp(:,:) = tmp(:,:) + FtSSHw(:,:,l) * FtFw(l,k)
    enddo
    FtSSHw(:,:,k) = (FtSSHw(:,:,k) - tmp(:,:)) / FtFw(k,k)
  enddo

  deallocate(tmp)
  deallocate(FtFw)

end subroutine HA_solver

end module MOM_HA

