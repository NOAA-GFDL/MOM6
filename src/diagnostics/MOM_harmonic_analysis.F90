!> Inline harmonic analysis (conventional)
module MOM_harmonic_analysis

use MOM_time_manager,  only : time_type, real_to_time, time_type_to_real, get_date, increment_date, &
                              operator(+), operator(-), operator(<), operator(>), operator(>=)
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_file_parser,   only : param_file_type, get_param
use MOM_io,            only : file_exists, open_ASCII_file, READONLY_FILE, close_file, &
                              MOM_infra_file, vardesc, MOM_field, &
                              var_desc, create_MOM_file, SINGLE_FILE, MOM_write_field
use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE

implicit none ; private

public HA_init, HA_register, HA_accum_FtF, HA_accum_FtSSH

#include <MOM_memory.h>

integer, parameter :: MAX_CONSTITUENTS = 10  !< The maximum number of tidal constituents

!> The private control structure for storing the HA info of a particular field
type, private :: HA_type
  character(len=16) :: key = "none"          !< Name of the field of which harmonic analysis is to be performed
  real :: old_time = -1.0                    !< The time of the previous accumulating step [T ~> s]
  real, allocatable :: ref(:,:)              !< The initial field in arbitrary units [A]
  real, allocatable :: FtSSH(:,:,:)          !< Accumulator of (F' * SSH_in) in arbitrary units [A]
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
end type HA_type

!> A linked list of control structures that store the HA info of different fields
type, private :: HA_node
  type(HA_type)          :: this             !< Control structure of the current field in the list
  type(HA_node), pointer :: next             !< The list of other fields
end type HA_node

!> The public control structure of the MOM_harmonic_analysis module
type, public :: harmonic_analysis_CS ; private
  logical :: HAready = .true.                !< If true, perform harmonic analysis
  type(time_type) :: &
    time_start, &                            !< Start time of harmonic analysis
    time_end, &                              !< End time of harmonic analysis
    time_ref                                 !< Reference time (t = 0) used to calculate tidal forcing
  real, dimension(MAX_CONSTITUENTS) :: &
    freq, &                                  !< The frequency of a tidal constituent [T-1 ~> s-1]
    phase0                                   !< The phase of a tidal constituent at time 0 [rad]
  real, allocatable :: FtF(:,:)              !< Accumulator of (F' * F) for all fields [nondim]
  integer :: nc                              !< The number of tidal constituents in use
  integer :: length                          !< Number of fields of which harmonic analysis is to be performed
  character(len=16)  :: const_name(MAX_CONSTITUENTS) !< The name of each constituent
  character(len=255) :: path                 !< Path to directory where output will be written
  type(HA_node), pointer :: list => NULL()   !< A linked list for storing the HA info of different fields
end type harmonic_analysis_CS

contains

!> This subroutine sets static variables used by this module and initializes CS%list.
!! THIS MUST BE CALLED AT THE END OF tidal_forcing_init.
subroutine HA_init(Time, US, param_file, time_ref, nc, freq, phase0, const_name, CS)
  type(time_type),       intent(in)  :: Time, &     !< The current model time
                                        time_ref    !< Reference time (t = 0) used to calculate tidal forcing
  type(unit_scale_type), intent(in)  :: US          !< A dimensional unit scaling type
  type(param_file_type), intent(in)  :: param_file  !< A structure to parse for run-time parameters
  real, dimension(MAX_CONSTITUENTS), &
                         intent(in)  :: freq, &     !< The frequency of a tidal constituent [T-1 ~> s-1]
                                        phase0      !< The phase of a tidal constituent at time 0 [rad]
  integer,               intent(in)  :: nc          !< The number of tidal constituents in use
  character(len=16),     intent(in)  :: const_name(MAX_CONSTITUENTS) !< The name of each constituent
  type(harmonic_analysis_CS), intent(out) :: CS     !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  type(HA_type) :: ha1                              !< A temporary, null field used for initializing CS%list
  real :: HA_start_time                             !< Start time of harmonic analysis [T ~> s]
  real :: Time_unit                                 !< The time unit for CS%time_end [T ~> s]
  character(len=40)  :: mdl="MOM_harmonic_analysis" !< This module's name
  character(len=255) :: mesg
  integer :: year, month, day, hour, minute, second

  integer :: unit, io_status
  integer :: date_init(6)=0                         !< The start date of the whole simulation.
  character(len=16) :: calendar = 'julian'          !< The name of the calendar type.
  integer :: years=0, months=0, days=0              !< These may determine the segment run
  integer :: hours=0, minutes=0, seconds=0          !! length, if read from a namelist

  namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds

  ! Determine CS%time_end by end time of the run segment
  if (file_exists('input.nml')) then
    call open_ASCII_file(unit, 'input.nml', action=READONLY_FILE)
    read(unit, ocean_solo_nml, iostat=io_status)
    call close_file(unit)
  endif

  if (years+months+days+hours+minutes+seconds > 0) then
    CS%time_end = increment_date(Time, years, months, days, hours, minutes, seconds)
    call MOM_mesg('HAmod: end time of harmonic analysis determined from ocean_solo_nml.', 2)
  else
    call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
                   "The time unit for DAYMAX, ENERGYSAVEDAYS, and RESTINT.", &
                   units="s", default=86400.0)
    call get_param(param_file, mdl, "DAYMAX", CS%time_end, &
                   "The final time of the whole simulation, in units of "//&
                   "TIMEUNIT seconds.  This also sets the potential end "//&
                   "time of the present run segment if the end time is "//&
                   "not set via ocean_solo_nml in input.nml.", &
                   timeunit=Time_unit, fail_if_missing=.true.)
  endif

  ! Determine CS%time_start
  call get_param(param_file, mdl, "HA_START_TIME", HA_start_time, &
                 "Start time of harmonic analysis (HA), in units of days after "//&
                 "the start of the current run segment. If equal to or greater than "//&
                 "the length of the run segment, HA will not be performed. "//&
                 "If negative, |HA_start_time| determines the length of simulation time "//&
                 "over which HA will be performed. In this case, HA will start |HA_start_time| days "//&
                 "before the end of the run segment, or at the beginning of the run segment, "//&
                 "whichever occurs later. HA will always finish at the end of the run segment.", &
                 units="days", default=0.0, scale=86400.0*US%s_to_T, do_not_log=.True., fail_if_missing=.false.)

  if (HA_start_time >= 0.0) then
    CS%time_start = Time + real_to_time(US%T_to_s * HA_start_time)
    if (CS%time_start >= CS%time_end) then
      call MOM_mesg('HAmod: HA_start_time equal to or greater than the length of run segment, '//&
                    'harmonic analysis will not be performed.')
      CS%HAready = .false. ; return
    endif
  else
    HA_start_time = 0.0 - HA_start_time
    if (time_type_to_real(CS%time_end - Time) >= HA_start_time) then
      CS%time_start = CS%time_end - real_to_time(US%T_to_s * HA_start_time)
    else
      call MOM_mesg('HAmod: harmonic analysis will be performed over the entire run segment.')
      CS%time_start = Time
    endif
  endif

  call get_date(Time, year, month, day, hour, minute, second)
  write(mesg,*) "HAmod: run segment starts on ", year, month, day, hour, minute, second
  call MOM_error(NOTE, trim(mesg))
  call get_date(CS%time_start, year, month, day, hour, minute, second)
  write(mesg,*) "HAmod: harmonic analysis starts on ", year, month, day, hour, minute, second
  call MOM_error(NOTE, trim(mesg))
  call get_date(CS%time_end, year, month, day, hour, minute, second)
  write(mesg,*) "HAmod: harmonic analysis ends on ", year, month, day, hour, minute, second
  call MOM_error(NOTE, trim(mesg))

  ! Set path to directory where output will be written
  call get_param(param_file, mdl, "HA_path", CS%path, &
                 "Path to output files for runtime harmonic analysis.", &
                 default="./", fail_if_missing=.false.)

  ! Populate some parameters of the control structure
  CS%time_ref   =  time_ref
  CS%freq       =  freq
  CS%phase0     =  phase0
  CS%nc         =  nc
  CS%const_name =  const_name
  CS%length     =  0

  allocate(CS%FtF(2*nc+1,2*nc+1), source=0.0)

  ! Initialize CS%list
  allocate(CS%list)
  CS%list%this  =  ha1
  nullify(CS%list%next)

end subroutine HA_init

!> This subroutine registers each of the fields on which HA is to be performed.
subroutine HA_register(key, CS)
  character(len=*),           intent(in)    :: key     !< Name of the current field
  type(harmonic_analysis_CS), intent(inout) :: CS      !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  type(HA_type)          :: ha1                        !< Control structure for the current field
  type(HA_node), pointer :: tmp                        !< A temporary list to hold the current field

  if (.not. CS%HAready) return

  allocate(tmp)
  ha1%key   =  trim(key)
  tmp%this  =  ha1
  tmp%next  => CS%list
  CS%list   => tmp
  CS%length =  CS%length + 1

end subroutine HA_register

!> This subroutine accumulates the temporal basis functions in FtF.
!! The tidal constituents are those used in MOM_tidal_forcing, plus the mean (of zero frequency).
subroutine HA_accum_FtF(Time, US, CS)
  type(time_type),            intent(in)    :: Time    !< The current model time
  type(unit_scale_type),      intent(in)    :: US      !< A dimensional unit scaling type
  type(harmonic_analysis_CS), intent(inout) :: CS      !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  real :: now                                          !< The relative time compared with tidal reference [T ~> s]
  real :: cosomegat, sinomegat, ccosomegat, ssinomegat !< The components of the phase [nondim]
  integer :: nc, c, icos, isin, cc, iccos, issin

  ! Exit the accumulator in the following cases
  if (.not. CS%HAready) return
  if (CS%length == 0) return
  if (Time < CS%time_start) return

  nc  = CS%nc
  now = US%s_to_T * time_type_to_real(Time - CS%time_ref)

  ! Accumulate FtF
  CS%FtF(1,1) = CS%FtF(1,1) + 1.0         !< For the zero frequency
  do c=1,nc
    icos = 2*c
    isin = 2*c+1
    cosomegat = cos(CS%freq(c) * now + CS%phase0(c))
    sinomegat = sin(CS%freq(c) * now + CS%phase0(c))
    CS%FtF(icos,1) = CS%FtF(icos,1) + cosomegat
    CS%FtF(isin,1) = CS%FtF(isin,1) + sinomegat
    CS%FtF(1,icos) = CS%FtF(icos,1)
    CS%FtF(1,isin) = CS%FtF(isin,1)
    do cc=c,nc
      iccos = 2*cc
      issin = 2*cc+1
      ccosomegat = cos(CS%freq(cc) * now + CS%phase0(cc))
      ssinomegat = sin(CS%freq(cc) * now + CS%phase0(cc))
      CS%FtF(icos,iccos) = CS%FtF(icos,iccos) + cosomegat * ccosomegat
      CS%FtF(icos,issin) = CS%FtF(icos,issin) + cosomegat * ssinomegat
      CS%FtF(isin,iccos) = CS%FtF(isin,iccos) + sinomegat * ccosomegat
      CS%FtF(isin,issin) = CS%FtF(isin,issin) + sinomegat * ssinomegat
    enddo ! cc=c,nc
  enddo ! c=1,nc

end subroutine HA_accum_FtF

!> This subroutine accumulates the temporal basis functions in FtSSH and then calls HA_write to compute
!! harmonic constants and write results. The tidal constituents are those used in MOM_tidal_forcing, plus the
!! mean (of zero frequency).
subroutine HA_accum_FtSSH(key, data, Time, G, US, CS)
  character(len=*),           intent(in) :: key  !< Name of the current field
  real, dimension(:,:),       intent(in) :: data !< Input data of which harmonic analysis is to be performed [A]
  type(time_type),            intent(in) :: Time !< The current model time
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(harmonic_analysis_CS), intent(in) :: CS   !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  type(HA_type), pointer :: ha1
  type(HA_node), pointer :: tmp
  real :: now                                    !< The relative time compared with the tidal reference [T ~> s]
  real :: dt                                     !< The current time step size of the accumulator [T ~> s]
  real :: cosomegat, sinomegat                   !< The components of the phase [nondim]
  integer :: nc, i, j, k, c, icos, isin, is, ie, js, je
  character(len=128) :: mesg

  ! Exit the accumulator in the following cases
  if (.not. CS%HAready) return
  if (CS%length == 0) return
  if (Time < CS%time_start) return

  ! Loop through the full list to find the current field
  tmp => CS%list
  do k=1,CS%length
    ha1 => tmp%this
    if (trim(key) == trim(ha1%key)) exit
    tmp => tmp%next
    if (k == CS%length) return              !< Do not perform harmonic analysis of a field that is not registered
  enddo

  nc  = CS%nc
  now = US%s_to_T * time_type_to_real(Time - CS%time_ref)

  ! Additional processing at the initial accumulating step
  if (ha1%old_time < 0.0) then
    ha1%old_time = now

    write(mesg,*) "HAmod: initializing accumulator, key = ", trim(ha1%key)
    call MOM_error(NOTE, trim(mesg))

    ! Get the lower and upper bounds of input data
    ha1%is = LBOUND(data,1) ; is = ha1%is
    ha1%ie = UBOUND(data,1) ; ie = ha1%ie
    ha1%js = LBOUND(data,2) ; js = ha1%js
    ha1%je = UBOUND(data,2) ; je = ha1%je

    allocate(ha1%ref(is:ie,js:je), source=0.0)
    allocate(ha1%FtSSH(is:ie,js:je,2*nc+1), source=0.0)
    ha1%ref(:,:) = data(:,:)
  endif

  dt = now - ha1%old_time
  ha1%old_time = now                        !< Keep track of time so we know when Time approaches CS%time_end

  is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  ! Accumulate FtF and FtSSH
  do c=1,nc
    icos = 2*c
    isin = 2*c+1
    cosomegat = cos(CS%freq(c) * now + CS%phase0(c))
    sinomegat = sin(CS%freq(c) * now + CS%phase0(c))
    do j=js,je ; do i=is,ie
      ha1%FtSSH(i,j,1)    = ha1%FtSSH(i,j,1)    + (data(i,j) - ha1%ref(i,j))
      ha1%FtSSH(i,j,icos) = ha1%FtSSH(i,j,icos) + (data(i,j) - ha1%ref(i,j)) * cosomegat
      ha1%FtSSH(i,j,isin) = ha1%FtSSH(i,j,isin) + (data(i,j) - ha1%ref(i,j)) * sinomegat
    enddo ; enddo
  enddo ! c=1,nc

  ! Compute harmonic constants and write output as Time approaches CS%time_end
  ! This guarantees that HA_write will be called before Time becomes larger than CS%time_end
  if (time_type_to_real(CS%time_end - Time) <= dt) then
    call HA_write(ha1, Time, G, CS)

    write(mesg,*) "HAmod: harmonic analysis done, key = ", trim(ha1%key)
    call MOM_error(NOTE, trim(mesg))

    deallocate(ha1%ref)
    deallocate(ha1%FtSSH)
  endif

end subroutine HA_accum_FtSSH

!> This subroutine computes the harmonic constants and write output for the current field
subroutine HA_write(ha1, Time, G, CS)
  type(HA_type), pointer,     intent(in) :: ha1
  type(time_type),            intent(in) :: Time   !< The current model time
  type(ocean_grid_type),      intent(in) :: G      !< The ocean's grid structure
  type(harmonic_analysis_CS), intent(in) :: CS     !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  real, dimension(:,:,:), allocatable :: FtSSHw    !< An array containing the harmonic constants [A]
  integer :: year, month, day, hour, minute, second
  integer :: nc, k, is, ie, js, je

  character(len=255)           :: filename         !< Output file name
  type(MOM_infra_file)         :: cdf              !< The file handle for output harmonic constants
  type(vardesc),   allocatable :: cdf_vars(:)      !< Output variable names
  type(MOM_field), allocatable :: cdf_fields(:)    !< Field type variables for the output fields

  nc = CS%nc ; is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  allocate(FtSSHw(is:ie,js:je,2*nc+1), source=0.0)

  ! Compute the harmonic coefficients
  call HA_solver(ha1, nc, CS%FtF, FtSSHw)

  ! Output file name
  call get_date(Time, year, month, day, hour, minute, second)
  write(filename, '(a,"HA_",a,i0.4,i0.2,i0.2,".nc")') &
      trim(CS%path), trim(ha1%key), year, month, day

  allocate(cdf_vars(2*nc+1))
  allocate(cdf_fields(2*nc+1))

  ! Variable names
  cdf_vars(1) = var_desc("z0", "m" ,"mean value", 'h', '1')
  do k=1,nc
    cdf_vars(2*k  ) = var_desc(trim(CS%const_name(k))//"cos", "m", "cosine coefficient", 'h', '1')
    cdf_vars(2*k+1) = var_desc(trim(CS%const_name(k))//"sin", "m", "sine coefficient",   'h', '1')
  enddo

  ! Create output file
  call create_MOM_file(cdf, trim(filename), cdf_vars, &
                       2*nc+1, cdf_fields, SINGLE_FILE, 86400.0, G=G)

  ! Write data
  call MOM_write_field(cdf, cdf_fields(1), G%domain, FtSSHw(:,:,1), 0.0)
  do k=1,nc
    call MOM_write_field(cdf, cdf_fields(2*k  ), G%domain, FtSSHw(:,:,2*k  ), 0.0)
    call MOM_write_field(cdf, cdf_fields(2*k+1), G%domain, FtSSHw(:,:,2*k+1), 0.0)
  enddo

  call cdf%flush()
  deallocate(cdf_vars)
  deallocate(cdf_fields)
  deallocate(FtSSHw)

end subroutine HA_write

!> This subroutine computes the harmonic constants (stored in FtSSHw) using the dot products of the temporal
!! basis functions accumulated in FtF, and the dot products of the SSH (or other fields) with the temporal basis
!! functions accumulated in FtSSH. The system is solved by Cholesky decomposition.
subroutine HA_solver(ha1, nc, FtF, FtSSHw)
  type(HA_type), pointer,              intent(in)  :: ha1
  integer,                             intent(in)  :: nc
  real, dimension(:,:),                intent(in)  :: FtF    !< Accumulator of (F' * F) for all fields [nondim]
  real, dimension(:,:,:), allocatable, intent(out) :: FtSSHw !< Work array for Cholesky decomposition [A]

  ! Local variables
  real, dimension(:,:), allocatable :: tmp                   !< Work array for Cholesky decomposition [A]
  real, dimension(:,:), allocatable :: FtFw                  !< Work array for Cholesky decomposition [nondim]
  integer :: k, l, is, ie, js, je

  is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  allocate(tmp(is:ie,js:je), source=0.0)
  allocate(FtFw(1:2*nc+1,1:2*nc+1), source=0.0)
  allocate(FtSSHw(is:ie,js:je,2*nc+1), source=0.0)

  FtFw(:,:) = 0.0
  do l=1,2*nc+1
    FtFw(l,l) = sqrt(FtF(l,l) - dot_product(FtFw(l,1:l-1), FtFw(l,1:l-1)))
    do k=l+1,2*nc+1
      FtFw(k,l) = (FtF(k,l) - dot_product(FtFw(k,1:l-1), FtFw(l,1:l-1))) / FtFw(l,l)
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

!> \namespace harmonic_analysis
!!
!! This module computes the harmonic constants which can be used to reconstruct the tidal elevation (or other
!! fields) through SSH = F * x, where F is an nt-by-2*nc matrix (nt is the number of time steps and nc is the
!! number of tidal constituents) containing the cosine/sine functions for each frequency evaluated at each time
!! step, and x is a 2*nc-by-1 vector containing the constant coefficients of the sine/cosine for each constituent
!! (i.e., the harmonic constants). At each grid point, the harmonic constants are computed using least squares,
!!
!!     x = (F' * F)^{-1} * (F' * SSH_in),
!!
!! where the prime denotes matrix transpose, and SSH_in is the sea surface height (or other fields) determined by
!! the model. The dot products (F' * F) and (F' * SSH_in) are computed by accumulating the sums as the model is
!! running and stored in the arrays FtF and FtSSH, respectively. The FtF matrix is inverted as needed before
!! computing and writing out the harmonic constants.
!!
!! Ed Zaron and William Xu (chengzhu.xu@oregonstate.edu), April 2024.

end module MOM_harmonic_analysis

