!> Provides diagnostics of work due to a given diffusivity
module MOM_diagnose_kdwork

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public diagnoseKdWork

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains
!> Diagnose the implied "work", or buoyancy forcing & its integral, due to a given diffusivity and column state.
  subroutine diagnoseKdWork(G, GV, US, N2, Kd, Bdif_flx, dz, Bdif_flx_dz)
  type(ocean_grid_type),   intent(in) :: G     !< Grid type
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
       intent(in)                     :: N2    !< Buoyancy frequency [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                         Kd    !< Diffusivity [H2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
       intent(out)                    :: Bdif_flx !< Buoyancy flux [H2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in), optional           :: dz    !< Grid spacing [H ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(out), optional          :: Bdif_flx_dz !< Buoyancy flux over layer [H3 T-3 ~> m3 s-3]

  integer :: i, j, k

  !$OMP parallel do default(shared)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    Bdif_flx(i,j,1) = 0.0
    Bdif_flx(i,j,GV%ke+1) = 0.0
    do K=2,GV%ke
      Bdif_flx(i,j,K) = - N2(i,j,K) * Kd(i,j,K)
    enddo
  enddo; enddo

  if (present(Bdif_flx_dz)) then
    !$OMP parallel do default(shared)
    if (.not. present(dz)) &
      call MOM_error(FATAL,'diagnoseKdWork called requesting integrated output but not passing dz')
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      do K=1,GV%ke
        Bdif_flx_dz(i,j,k) = (Bdif_flx(i,j,K)+Bdif_flx(i,j,K+1))*dz(i,j,k)
      enddo
    enddo; enddo
  endif

end subroutine diagnoseKdWork

!> \namespace mom_diagnose_kdwork
!!
!!    The subroutine diagnoseKdWork diagnoses the energetics associated with various vertical diffusivities
!!    inside MOM6 diabatic routines.
!!

end module MOM_diagnose_kdwork
